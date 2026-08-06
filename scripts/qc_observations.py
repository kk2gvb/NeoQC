"""Extract normalized QC observations from NeoQC TSV artifacts."""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Callable


class ObservationError(ValueError):
    """Raised when an input artifact cannot produce a valid observation."""


NumericRows = list[dict[str, float]]
Extractor = Callable[[Path, str], dict[str, float]]


SOURCE_PREFIXES = {
    "per_base_quality": "per_cycle",
    "per_sequence_quality": "per_sequence_quality",
    "per_base_sequence_content": "per_base_sequence_content",
    "per_sequence_gc_content": "per_sequence_gc_content",
    "per_base_n_content": "per_base_n_content",
    "sequence_length_distribution": "sequence_length_distribution",
    "sequence_duplication_levels": "sequence_duplication_levels",
    "adapter_content": "adapter_content",
}


def source_path(input_dir: Path, metric_id: str, read: str) -> Path:
    try:
        prefix = SOURCE_PREFIXES[metric_id]
    except KeyError as error:
        raise ObservationError(f"unsupported metric: {metric_id}") from error
    return input_dir / f"{prefix}_{read}.tsv"


def _read_numeric(
    path: Path,
    required: tuple[str, ...],
    *,
    allow_extra: bool = False,
) -> tuple[NumericRows, tuple[str, ...]]:
    try:
        handle = path.open("r", encoding="utf-8", newline="")
    except OSError as error:
        raise ObservationError(f"cannot read {path.name}: {error}") from error
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = tuple(reader.fieldnames or ())
        if not columns or any(not column.strip() for column in columns):
            raise ObservationError(f"{path.name}: invalid header")
        if len(set(columns)) != len(columns):
            raise ObservationError(f"{path.name}: duplicate column names")
        missing = [column for column in required if column not in columns]
        if missing:
            raise ObservationError(
                f"{path.name}: missing required column(s): {', '.join(missing)}"
            )
        if not allow_extra:
            unexpected = [column for column in columns if column not in required]
            if unexpected:
                raise ObservationError(
                    f"{path.name}: unexpected column(s): {', '.join(unexpected)}"
                )
        rows: NumericRows = []
        for line_number, raw in enumerate(reader, start=2):
            if None in raw:
                raise ObservationError(f"{path.name}:{line_number}: too many fields")
            row: dict[str, float] = {}
            for column in columns:
                value = (raw.get(column) or "").strip()
                try:
                    number = float(value)
                except ValueError as error:
                    raise ObservationError(
                        f"{path.name}:{line_number}: {column} is not numeric"
                    ) from error
                if not math.isfinite(number):
                    raise ObservationError(
                        f"{path.name}:{line_number}: {column} is not finite"
                    )
                if number < 0:
                    raise ObservationError(
                        f"{path.name}:{line_number}: {column} must not be negative"
                    )
                row[column] = number
            rows.append(row)
    if not rows:
        raise ObservationError(f"{path.name}: contains no data rows")
    return rows, columns


def _range(values: list[float], label: str, lower: float, upper: float) -> None:
    if any(value < lower or value > upper for value in values):
        raise ObservationError(f"{label} must be between {lower:g} and {upper:g}")


def _per_base_quality(path: Path, _read: str) -> dict[str, float]:
    rows, _ = _read_numeric(
        path, ("cycle", "mean_quality", "lower_quartile", "median")
    )
    return {
        "minimum_lower_quartile": min(row["lower_quartile"] for row in rows),
        "minimum_median": min(row["median"] for row in rows),
    }


def _per_sequence_quality(path: Path, _read: str) -> dict[str, float]:
    rows, _ = _read_numeric(path, ("mean_quality", "read_count"))
    total = sum(row["read_count"] for row in rows)
    if total <= 0:
        raise ObservationError("per-sequence quality contains no observations")
    mode = max(rows, key=lambda row: row["read_count"])["mean_quality"]
    return {"modal_mean_quality": mode}


def _base_content(path: Path, _read: str) -> dict[str, float]:
    rows, _ = _read_numeric(path, ("position", "A", "C", "G", "T", "N"))
    percentages = [row[base] for row in rows for base in ("A", "C", "G", "T", "N")]
    _range(percentages, "base content", 0, 100)
    difference = max(
        max(abs(row["A"] - row["T"]), abs(row["G"] - row["C"]))
        for row in rows
    )
    return {"maximum_base_difference_percent": difference}


def fastqc_theoretical_gc(
    x: list[float], counts: list[float]
) -> tuple[list[int], list[float], float, float]:
    """Return FastQC-style GC model coordinates, centre and deviation."""
    observed = [0.0] * 101
    for percentage, count in zip(x, counts):
        index = int(round(percentage))
        if not math.isclose(percentage, index, abs_tol=1e-6) or index < 0 or index > 100:
            raise ObservationError("GC percentage values must be whole numbers from 0 to 100")
        observed[index] += count
    total = sum(observed)
    if total <= 1:
        raise ObservationError("GC distribution needs at least two reads")
    first_mode = max(range(101), key=observed.__getitem__)
    threshold = observed[first_mode] * 0.9
    modal_bins = [first_mode]
    for index in range(first_mode + 1, 101):
        if observed[index] <= threshold:
            break
        modal_bins.append(index)
    for index in range(first_mode - 1, -1, -1):
        if observed[index] <= threshold:
            break
        modal_bins.append(index)
    touches_edge = min(modal_bins) == 0 or max(modal_bins) == 100
    centre = float(first_mode) if touches_edge else sum(modal_bins) / len(modal_bins)
    variance = sum(((index - centre) ** 2) * count for index, count in enumerate(observed))
    stdev = math.sqrt(variance / (total - 1))
    if stdev <= 0:
        theoretical = [0.0] * 101
        theoretical[int(round(centre))] = total
    else:
        scale = total / (math.sqrt(2.0 * math.pi) * stdev)
        theoretical = [
            scale * math.exp(-((index - centre) ** 2) / (2.0 * stdev * stdev))
            for index in range(101)
        ]
    deviation = (
        sum(abs(model - actual) for model, actual in zip(theoretical, observed))
        / total
        * 100.0
    )
    return list(range(101)), theoretical, centre, deviation


def modeled_gc_deviation(x: list[float], counts: list[float]) -> float:
    return fastqc_theoretical_gc(x, counts)[3]


def _gc_content(path: Path, _read: str) -> dict[str, float]:
    rows, _ = _read_numeric(path, ("gc_percent", "reads"))
    deviation = modeled_gc_deviation(
        [row["gc_percent"] for row in rows],
        [row["reads"] for row in rows],
    )
    return {"modeled_deviation_percent": deviation}


def _n_content(path: Path, _read: str) -> dict[str, float]:
    rows, _ = _read_numeric(path, ("position", "N_percent"))
    values = [row["N_percent"] for row in rows]
    _range(values, "N percentage", 0, 100)
    return {"maximum_n_percent": max(values)}


def _length_distribution(path: Path, _read: str) -> dict[str, float]:
    rows, _ = _read_numeric(path, ("length", "reads"))
    observed = sorted({row["length"] for row in rows if row["reads"] > 0})
    if not observed:
        raise ObservationError("length distribution contains no observations")
    return {
        "observed_length_count": float(len(observed)),
        "minimum_length": min(observed),
    }


def _duplication(path: Path, _read: str) -> dict[str, float]:
    try:
        handle = path.open("r", encoding="utf-8", newline="")
    except OSError as error:
        raise ObservationError(f"cannot read {path.name}: {error}") from error
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected = (
            "duplication_level",
            "total_sequences_percent",
            "deduplicated_sequences_percent",
        )
        if tuple(reader.fieldnames or ()) != expected:
            raise ObservationError(f"{path.name}: expected columns: {', '.join(expected)}")
        level_one: tuple[float, float] | None = None
        seen: set[str] = set()
        for line_number, raw in enumerate(reader, start=2):
            label = (raw.get("duplication_level") or "").strip()
            if not label or label in seen:
                raise ObservationError(f"{path.name}:{line_number}: invalid duplication level")
            seen.add(label)
            try:
                total = float(raw["total_sequences_percent"] or "")
                deduplicated = float(raw["deduplicated_sequences_percent"] or "")
            except (TypeError, ValueError) as error:
                raise ObservationError(f"{path.name}:{line_number}: percentages must be numeric") from error
            if not math.isfinite(total) or not math.isfinite(deduplicated):
                raise ObservationError(f"{path.name}:{line_number}: percentages must be finite")
            _range([total, deduplicated], "duplication percentage", 0, 100)
            if label in {"1", "1.0"}:
                level_one = (total, deduplicated)
    if level_one is None:
        raise ObservationError("duplication level 1 is required for the headline statistic")
    total_one, deduplicated_one = level_one
    if deduplicated_one <= 0:
        raise ObservationError("deduplicated percentage at level 1 must be positive")
    remaining = 100.0 * total_one / deduplicated_one
    if remaining > 100.0 + 1e-6:
        raise ObservationError("duplication level percentages are inconsistent")
    remaining = min(100.0, remaining)
    return {
        "deduplicated_remaining_percent": remaining,
        "deduplicated_loss_percent": 100.0 - remaining,
    }


def _adapter_content(path: Path, _read: str) -> dict[str, float]:
    rows, columns = _read_numeric(path, ("pos",), allow_extra=True)
    series = [column for column in columns if column != "pos"]
    if not series:
        raise ObservationError("adapter content contains no adapter series")
    values = [row[column] for row in rows for column in series]
    _range(values, "adapter percentage", 0, 100)
    return {"maximum_adapter_percent": max(values)}


EXTRACTORS: dict[str, Extractor] = {
    "per_base_quality": _per_base_quality,
    "per_sequence_quality": _per_sequence_quality,
    "per_base_sequence_content": _base_content,
    "per_sequence_gc_content": _gc_content,
    "per_base_n_content": _n_content,
    "sequence_length_distribution": _length_distribution,
    "sequence_duplication_levels": _duplication,
    "adapter_content": _adapter_content,
}


def extract_observations(input_dir: Path, metric_id: str, read: str) -> dict[str, float]:
    try:
        extractor = EXTRACTORS[metric_id]
    except KeyError as error:
        raise ObservationError(f"unsupported metric: {metric_id}") from error
    return extractor(source_path(input_dir, metric_id, read), read)
