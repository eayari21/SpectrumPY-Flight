#!/usr/bin/env python3
"""
drive_idex_packet.py — fast, parallel orchestrator for idex_packet.py

Examples:
  # Basic: scan ./Data, write to ./HDF5, use all CPUs
  python drive_idex_packet.py --inputs Data --out HDF5

  # Limit to 8 concurrent jobs and 2 BLAS threads per job
  python drive_idex_packet.py --inputs Data --out HDF5 --max-procs 8 --threads-per-proc 2

  # Multiple input roots
  python drive_idex_packet.py --inputs Pre_Env Post_Env --out HDF5

  # Dry-run (show what would run)
  python drive_idex_packet.py --inputs Data --out HDF5 --dry-run
"""
from __future__ import annotations
import argparse, os, sys, subprocess, time, shutil, importlib.util
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Tuple

def find_source_files(roots: List[Path]) -> List[Path]:
    """Return all files with **no extension** under the given roots."""
    files: List[Path] = []
    for root in roots:
        if not root.exists():
            continue
        for p in root.rglob("*"):
            if not p.is_file():
                continue
            if p.suffix != "":  # Skip anything that already has an extension
                continue
            if p.name.startswith('.'):
                # Ignore metadata files like .DS_Store that lack extensions
                continue
            files.append(p)
    return files

def needs_conversion(src: Path, out_dir: Path) -> Tuple[bool, Path]:
    dst = out_dir / (src.name + ".h5")
    return (not dst.exists()), dst

def build_env(threads_per_proc: int | None) -> dict:
    env = os.environ.copy()
    if threads_per_proc:
        t = str(threads_per_proc)
        # Cap common BLAS/OpenMP envs so parallel jobs don’t oversubscribe cores
        env["OMP_NUM_THREADS"] = t
        env["OPENBLAS_NUM_THREADS"] = t
        env["MKL_NUM_THREADS"] = t
        env["NUMEXPR_NUM_THREADS"] = t
        env["VECLIB_MAXIMUM_THREADS"] = t  # Accelerate (macOS)
    return env

def run_one(idx: int, py: str, entrypoint: tuple[str, str | Path], src: Path, env: dict,
            log_dir: Path, nice: int | None) -> Tuple[Path, int, float]:
    """Run the idex_packet converter and return ``(src, rc, seconds)``."""
    log_dir.mkdir(parents=True, exist_ok=True)
    out_log = log_dir / f"job_{idx:04d}.out"
    err_log = log_dir / f"job_{idx:04d}.err"

    kind, target = entrypoint
    cmd = [py]
    if kind == "module":
        cmd.extend(["-m", str(target)])
    else:
        cmd.append(str(target))
    cmd.extend(["-f", str(src)])
    start = time.time()

    # Optional niceness (POSIX); ignore on platforms without `nice`
    pre_cmd = []
    if nice is not None and shutil.which("nice"):
        pre_cmd = ["nice", f"-n{nice}"]

    with open(out_log, "wb") as out_f, open(err_log, "wb") as err_f:
        out_f.write(f"$ {' '.join(pre_cmd + cmd)}\n".encode())
        out_f.flush()
        proc = subprocess.Popen(pre_cmd + cmd, stdout=out_f, stderr=err_f, env=env)
        rc = proc.wait()
    dur = time.time() - start
    return (src, rc, dur)

def locate_idex_entrypoint(explicit: str | None) -> tuple[str, str | Path]:
    """Locate the best way to invoke ``idex_packet``.

    Returns a tuple of ("module"|"script", module_name or Path).
    """

    if explicit:
        candidate = Path(explicit)
        if candidate.is_file():
            return ("script", candidate.resolve())
        # Allow passing a file relative to this module
        relative = Path(__file__).resolve().parent / explicit
        if relative.is_file():
            return ("script", relative)
        raise FileNotFoundError(f"Cannot locate idex script at '{explicit}'")

    # Prefer the installed module when available; works for pip installs.
    try:
        spec = importlib.util.find_spec("spectrumpy_flight.idex_packet")
    except ModuleNotFoundError:
        spec = None
    if spec is not None:
        return ("module", "spectrumpy_flight.idex_packet")

    # Fallback to a sibling file, which covers running from a source checkout/zip.
    sibling = Path(__file__).with_name("idex_packet.py")
    if sibling.is_file():
        return ("script", sibling.resolve())

    raise FileNotFoundError("Unable to find spectrumpy_flight.idex_packet entrypoint")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", nargs="+", required=True, help="Input directories to scan.")
    ap.add_argument("--out", required=True, help="Output HDF5 directory (where .h5 are created).")
    ap.add_argument("--idex-script", default=None, help="Path to idex_packet.py (optional).")
    ap.add_argument("--python", default=sys.executable, help="Python interpreter to run idex_packet.py")
    ap.add_argument("--max-procs", type=int, default=os.cpu_count() or 4, help="Max concurrent conversions.")
    ap.add_argument("--threads-per-proc", type=int, default=None, help="Thread cap inside each conversion.")
    ap.add_argument("--log-dir", default="idex_logs", help="Directory for per-job logs.")
    ap.add_argument("--nice", type=int, default=None, help="POSIX niceness (e.g., 10).")
    ap.add_argument("--dry-run", action="store_true", help="Print planned work only.")
    args = ap.parse_args()

    inputs = [Path(p).resolve() for p in args.inputs]
    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    entrypoint = locate_idex_entrypoint(args.idex_script)
    log_dir = Path(args.log_dir).resolve()

    # 1) Discover work
    candidates = find_source_files(inputs)
    to_do: List[Tuple[int, Path]] = []
    for i, src in enumerate(sorted(candidates)):
        needed, dst = needs_conversion(src, out_dir)
        if needed:
            to_do.append((i, src))

    print(f"Discovered {len(candidates)} source files; {len(to_do)} need conversion.")
    if args.dry_run:
        for _, src in to_do[:20]:
            print(f"[DRY] would process: {src}")
        if len(to_do) > 20:
            print(f"... and {len(to_do) - 20} more")
        return

    if not to_do:
        print("Nothing to do. All .h5 outputs are present.")
        return

    # 2) Run in parallel
    env = build_env(args.threads_per_proc)
    max_workers = max(1, args.max_procs)
    target_desc = entrypoint[1] if entrypoint[0] == "script" else f"module {entrypoint[1]}"
    print(f"Running up to {max_workers} concurrent jobs "
          f"(threads-per-proc={args.threads_per_proc or 'default'}), using {target_desc}")

    ok, fail = 0, 0
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futs = [
            ex.submit(run_one, idx, args.python, entrypoint, src, env, log_dir, args.nice)
            for idx, src in to_do
        ]
        for fut in as_completed(futs):
            src, rc, dur = fut.result()
            if rc == 0:
                ok += 1
                print(f"[OK]  {src.name} in {dur:.1f}s")
            else:
                fail += 1
                print(f"[ERR] {src.name} exit={rc} (see logs)")
    dt = time.time() - t0
    print(f"\nDone in {dt:.1f}s — success: {ok}, failed: {fail}, total: {ok+fail}")

if __name__ == "__main__":
    main()
