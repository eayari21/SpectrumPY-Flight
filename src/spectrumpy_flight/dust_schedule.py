"""Utilities for loading dust accelerator test schedules."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from . import package_path

__all__ = ["DustScheduleEntry", "load_dust_schedule"]


@dataclass
class DustScheduleEntry:
    """Describe a dust accelerator campaign window."""

    __slots__ = ("label", "start", "end", "instrument_model", "material", "count")

    label: str
    start: datetime
    end: datetime
    instrument_model: str
    material: str
    count: int

    def contains(self, timestamp: datetime) -> bool:
        return self.start <= timestamp <= self.end


def _parse_day(text: str) -> int:
    digits = "".join(ch for ch in text if ch.isdigit())
    if not digits:
        raise ValueError(f"Could not parse day from '{text}'.")
    return int(digits)


def _parse_date_range(text: str) -> tuple[datetime, datetime]:
    cleaned = text.strip().replace(",", "")
    parts = cleaned.split()
    if not parts:
        raise ValueError("Empty date string in dust schedule.")

    if len(parts) == 2:
        month_name, year_text = parts
        month = datetime.strptime(month_name, "%B").month
        year = int(year_text)
        start = datetime(year, month, 1, tzinfo=timezone.utc)
        if month == 12:
            end_month = 1
            end_year = year + 1
        else:
            end_month = month + 1
            end_year = year
        end = datetime(end_year, end_month, 1, tzinfo=timezone.utc) - timedelta(seconds=1)
        return start, end

    if len(parts) == 3 and "-" in parts[1]:
        month_name, day_range, year_text = parts
        start_day_text, end_day_text = day_range.split("-", 1)
        month = datetime.strptime(month_name, "%B").month
        year = int(year_text)
        start_day = _parse_day(start_day_text)
        end_day = _parse_day(end_day_text)
        start = datetime(year, month, start_day, tzinfo=timezone.utc)
        end = datetime(year, month, end_day, 23, 59, 59, tzinfo=timezone.utc)
        return start, end

    raise ValueError(f"Unsupported dust schedule date format: '{text}'.")


def load_dust_schedule(path: Path | None = None) -> list[DustScheduleEntry]:
    """Load the dust testing schedule from a CSV file."""

    if path is None:
        path = package_path("lookup", "IDEX_Dust_Testing.csv")

    entries: list[DustScheduleEntry] = []
    with path.open("r", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            date_label = (row.get("Date") or "").strip()
            instrument_model = (row.get("Instrument Model") or "").strip()
            material = (row.get("Material") or "").strip()
            count_text = (row.get("Count") or "0").strip()
            if not date_label:
                continue
            start, end = _parse_date_range(date_label)
            try:
                count = int(count_text)
            except Exception:
                count = 0
            entries.append(
                DustScheduleEntry(
                    label=date_label,
                    start=start,
                    end=end,
                    instrument_model=instrument_model,
                    material=material,
                    count=count,
                )
            )
    entries.sort(key=lambda entry: entry.start)
    return entries
