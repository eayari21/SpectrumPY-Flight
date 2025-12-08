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
from typing import List, Sequence, Tuple

def find_source_files(roots: List[Path]) -> List[Path]:
    """Return all candidate input files under the given roots.

    Historically we only looked for extensionless files, but IDEX raw packet
    captures often arrive as ``.pkts`` or ``.bin``.  The extension does not
    matter to ``idex_packet.py`` so include everything except hidden files and
    existing ``.h5`` outputs.
    """
    files: List[Path] = []
    for root in roots:
        if not root.exists():
            continue
        for p in root.rglob("*"):
            if not p.is_file():
                continue
            if p.name.startswith('.'):
                # Ignore metadata files like .DS_Store that lack extensions
                continue
            if p.suffix.lower() == ".h5":
                # Skip already-produced outputs
                continue
            files.append(p)
    return files

def needs_conversion(src: Path, out_dir: Path) -> Tuple[bool, Path]:
    stem = src.stem if src.suffix else src.name
    dst = out_dir / f"{stem}.h5"
    return (not dst.exists()), dst

def _maybe_prepend_pythonpath(env: dict) -> None:
    """Ensure the spectrumpy_flight source tree is importable.

    When running from a source checkout, the package lives under ``src`` and
    is not necessarily installed.  Add that directory to ``PYTHONPATH`` so
    ``python -m spectrumpy_flight.idex_packet`` works in child processes.
    """

    here = Path(__file__).resolve()
    # ``.../src/spectrumpy_flight/drive_idex_packet.py`` -> ``.../src``
    src_root = here.parents[1]
    if not (src_root / "spectrumpy_flight").exists():
        return

    existing = env.get("PYTHONPATH")
    src_str = str(src_root)
    if existing:
        if src_str in existing.split(os.pathsep):
            return
        env["PYTHONPATH"] = os.pathsep.join([src_str, existing])
    else:
        env["PYTHONPATH"] = src_str


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
    _maybe_prepend_pythonpath(env)
    return env

def run_one(
    idx: int,
    py: str,
    entrypoint: tuple[str, str | Path],
    src: Path,
    dst: Path,
    env: dict,
    log_dir: Path,
    nice: int | None,
) -> Tuple[Path, int, float]:
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
    if rc == 0:
        produced = src.with_suffix(".h5")
        dst.parent.mkdir(parents=True, exist_ok=True)
        try:
            if produced.exists() and produced.resolve() != dst.resolve():
                shutil.move(str(produced), str(dst))
            elif produced.exists():
                # Already in the desired location
                pass
            elif not dst.exists():
                rc = 1
        except Exception:
            rc = 1
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
    module_name = "spectrumpy_flight.idex_packet"
    try:
        spec = importlib.util.find_spec(module_name)
    except ModuleNotFoundError:
        spec = None
    if spec is not None:
        return ("module", module_name)

    # Running from a source checkout: ensure the package directory exists and
    # prefer invoking it as a module so relative imports succeed.
    package_dir = Path(__file__).resolve().parent
    if (package_dir / "idex_packet.py").is_file():
        return ("module", module_name)

    raise FileNotFoundError("Unable to find spectrumpy_flight.idex_packet entrypoint")


def main(argv: Sequence[str] | None = None):
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
    args = ap.parse_args(argv)

    inputs = [Path(p).resolve() for p in args.inputs]
    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    entrypoint = locate_idex_entrypoint(args.idex_script)
    log_dir = Path(args.log_dir).resolve()

    # 1) Discover work
    candidates = find_source_files(inputs)
    to_do: List[Tuple[int, Path, Path]] = []
    for i, src in enumerate(sorted(candidates)):
        needed, dst = needs_conversion(src, out_dir)
        if needed:
            to_do.append((i, src, dst))

    print(f"Discovered {len(candidates)} source files; {len(to_do)} need conversion.")
    if args.dry_run:
        for _, src, dst in to_do[:20]:
            print(f"[DRY] would process: {src} -> {dst}")
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
            ex.submit(
                run_one,
                idx,
                args.python,
                entrypoint,
                src,
                dst,
                env,
                log_dir,
                args.nice,
            )
            for idx, src, dst in to_do
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
