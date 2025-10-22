"""Resource discovery utilities for SpectrumPY parallel workflows."""

from __future__ import annotations

import os
import shutil
from dataclasses import dataclass
from typing import Optional

try:
    import psutil  # type: ignore
except Exception:  # pragma: no cover - optional dependency
    psutil = None


@dataclass(frozen=True)
class SystemResources:
    """Snapshot of the host system resources relevant for parallel jobs."""

    logical_cpus: int
    available_memory: int
    gpu_count: int

    @property
    def available_memory_gb(self) -> float:
        """Return the available memory expressed in gigabytes."""

        return self.available_memory / (1024 ** 3)


def _detect_cpu_count() -> int:
    if psutil is not None:
        try:
            affinity = psutil.Process().cpu_affinity()  # type: ignore[attr-defined]
            if affinity:
                return len(affinity)
        except Exception:
            pass

    return os.cpu_count() or 1


def _detect_available_memory() -> int:
    if psutil is not None:
        try:
            return int(psutil.virtual_memory().available)
        except Exception:
            pass

    if hasattr(os, "sysconf"):
        try:
            page_size = os.sysconf("SC_PAGE_SIZE")
            available_pages = os.sysconf("SC_AVPHYS_PAGES")
            return int(page_size * available_pages)
        except (OSError, ValueError):
            pass

    return 0


def _detect_gpu_count() -> int:
    """Detect NVIDIA GPUs via ``nvidia-smi`` when available."""

    nvidia_smi = shutil.which("nvidia-smi")
    if not nvidia_smi:
        return 0

    try:
        completed = os.popen(f"{nvidia_smi} -L").read().strip()
    except Exception:
        return 0
    if not completed:
        return 0
    return completed.count("GPU ")


def system_resources() -> SystemResources:
    """Collect a snapshot of the host system capabilities."""

    return SystemResources(
        logical_cpus=_detect_cpu_count(),
        available_memory=_detect_available_memory(),
        gpu_count=_detect_gpu_count(),
    )


def recommended_worker_count(memory_per_job: Optional[int] = None, *, max_workers: Optional[int] = None) -> int:
    """Select a conservative worker count based on CPU and memory limits."""

    resources = system_resources()
    workers = resources.logical_cpus

    if max_workers is not None:
        workers = min(workers, max_workers)

    if memory_per_job and resources.available_memory:
        mem_limited = max(1, resources.available_memory // memory_per_job)
        workers = min(workers, mem_limited)

    return max(1, workers)
