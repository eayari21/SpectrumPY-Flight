"""Parallel packet processing entrypoint used by ``process_packets.sh``."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable, List, Sequence

from .resource_utils import recommended_worker_count, system_resources


def _resolve_processor(explicit: str | None) -> Path:
    candidates: List[Path] = []
    if explicit is not None:
        candidates.append(Path(explicit))
    candidates.extend([Path("lmfit_idex_packet.py"), Path("idex_packet.py")])
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        "Unable to locate an IDEX packet decoder. Checked: "
        + ", ".join(str(c) for c in candidates)
    )


def _iter_raw_files(paths: Sequence[Path]) -> Iterable[Path]:
    for base in paths:
        if base.is_file():
            yield base
            continue
        for candidate in sorted(base.rglob("*")):
            if candidate.is_file() and candidate.suffix == "":
                yield candidate


def _should_skip(output_dir: Path, raw_file: Path) -> bool:
    output_path = output_dir / f"{raw_file.name}.h5"
    return output_path.exists()


def _process_single(processor: Path, raw_file: Path, extra_args: Sequence[str]) -> tuple[Path, int, str]:
    cmd = [sys.executable, str(processor), "-f", str(raw_file), *extra_args]
    env = os.environ.copy()
    env.setdefault("NUMEXPR_MAX_THREADS", "auto")
    result = subprocess.run(cmd, capture_output=True, text=True, env=env)
    message = result.stdout + result.stderr
    return raw_file, result.returncode, message


def _log_system_profile() -> None:
    resources = system_resources()
    print(
        "[process_packets] Host profile: "
        f"CPUs={resources.logical_cpus}, "
        f"AvailableMemory={resources.available_memory_gb:.1f} GiB, "
        f"GPUs={resources.gpu_count}",
        flush=True,
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Parallel packet processing driver")
    parser.add_argument("paths", nargs="*", default=["Data"], help="Input directories or files")
    parser.add_argument("--processor", help="Explicit decoder to invoke per file")
    parser.add_argument("--output", default="HDF5", help="Directory storing the resulting HDF5 files")
    parser.add_argument(
        "--max-workers",
        type=int,
        default=None,
        help="Upper bound for worker processes (auto-detected when omitted)",
    )
    parser.add_argument(
        "--memory-per-job",
        type=int,
        default=2 * 1024 ** 3,
        help="Approximate memory budget per decoding job in bytes",
    )
    parser.add_argument(
        "--extra-arg",
        dest="extra_args",
        action="append",
        default=[],
        help="Additional arguments passed to the decoder (repeatable)",
    )

    args = parser.parse_args(argv)

    processor = _resolve_processor(args.processor)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    inputs = [Path(p).resolve() for p in args.paths]
    raw_files = [p for p in _iter_raw_files(inputs) if not _should_skip(output_dir, p)]

    if not raw_files:
        print("[process_packets] No new files detected — nothing to do.")
        return 0

    _log_system_profile()

    workers = recommended_worker_count(args.memory_per_job, max_workers=args.max_workers)
    print(f"[process_packets] Launching {workers} worker(s) using {processor}")

    errors: list[str] = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        future_map = {
            executor.submit(
                _process_single,
                processor,
                raw_file,
                args.extra_args,
            ): raw_file
            for raw_file in raw_files
        }

        for future in as_completed(future_map):
            raw_file, returncode, message = future.result()
            if returncode == 0:
                print(f"[process_packets] ✓ {raw_file.name} processed successfully")
            else:
                print(f"[process_packets] ✗ {raw_file.name} failed with exit code {returncode}")
                errors.append(message)

    if errors:
        print("\n".join(errors))
        return 1

    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entrypoint
    raise SystemExit(main())
