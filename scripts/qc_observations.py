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
    optional: tuple[str, ...] = (),
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
            allowed = required + optional
            unexpected = [column for column in columns if column not in allowed]
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
    rows, columns = _read_numeric(
        path,
        ("mean_quality", "read_count"),
        optional=("read_count_truncate",),
    )
    count_column = (
        "read_count_truncate" if "read_count_truncate" in columns else "read_count"
    )
    total = sum(row[count_column] for row in rows)
    if total <= 0:
        raise ObservationError("per-sequence quality contains no observations")
    mode = max(rows, key=lambda row: row[count_column])["mean_quality"]
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


def _native_duplication_summary(path: Path, read: str) -> float | None:
    marker = path.parent / f"sequence_duplication_{read}.incomplete"
    if marker.exists():
        raise ObservationError(
            f"{marker.name}: native duplication artifacts are incomplete"
        )

    summary = path.parent / f"sequence_duplication_summary_{read}.tsv"
    if not summary.exists():
        return None

    exact_columns = (
        "source_kind",
        "algorithm",
        "source_fastq",
        "prefix_length",
        "total_reads",
        "unique_sequences",
        "deduplicated_remaining_percent",
    )
    bounded_columns = (
        "source_kind",
        "algorithm",
        "source_fastq",
        "prefix_length",
        "max_tracked_unique",
        "total_reads",
        "tracked_unique_sequences",
        "count_at_unique_limit",
        "sampling_limited",
        "deduplicated_remaining_percent",
    )
    try:
        handle = summary.open("r", encoding="utf-8", newline="")
    except OSError as error:
        raise ObservationError(f"cannot read {summary.name}: {error}") from error
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = tuple(reader.fieldnames or ())
        if columns not in {exact_columns, bounded_columns}:
            raise ObservationError(
                f"{summary.name}: unsupported summary schema"
            )
        rows = list(reader)
    if len(rows) != 1 or None in rows[0]:
        raise ObservationError(f"{summary.name}: expected exactly one complete data row")
    row = {key: (value or "").strip() for key, value in rows[0].items()}
    if row["source_kind"] != "native_fastq":
        raise ObservationError(f"{summary.name}: unsupported source_kind")
    if not row["source_fastq"] or Path(row["source_fastq"]).name != row["source_fastq"]:
        raise ObservationError(f"{summary.name}: source_fastq must be a filename")

    integer_fields = tuple(
        column
        for column in columns
        if column
        in {
            "prefix_length",
            "max_tracked_unique",
            "total_reads",
            "unique_sequences",
            "tracked_unique_sequences",
            "count_at_unique_limit",
        }
    )
    try:
        integers = {field: int(row[field]) for field in integer_fields}
        remaining = float(row["deduplicated_remaining_percent"])
    except ValueError as error:
        raise ObservationError(f"{summary.name}: invalid numeric value") from error
    if any(value < 0 for value in integers.values()) or not math.isfinite(remaining):
        raise ObservationError(f"{summary.name}: invalid numeric range")
    if integers["prefix_length"] != 50:
        raise ObservationError(f"{summary.name}: unsupported prefix_length")
    if integers["total_reads"] <= 0:
        raise ObservationError(f"{summary.name}: total_reads must be positive")
    _range([remaining], "deduplicated remaining percentage", 0, 100)

    if columns == exact_columns:
        if row["algorithm"] != "neoqc-exact-prefix-v1":
            raise ObservationError(f"{summary.name}: unsupported exact algorithm")
        unique = integers["unique_sequences"]
        total = integers["total_reads"]
        if not 0 < unique <= total:
            raise ObservationError(f"{summary.name}: unique_sequences is inconsistent")
        calculated = 100.0 * unique / total
        if not math.isclose(remaining, calculated, abs_tol=1e-8):
            raise ObservationError(f"{summary.name}: exact headline is inconsistent")
        return remaining

    # Compatibility with artifacts produced by the earlier bounded prototype.
    if row["algorithm"] != "fastqc-compatible-bounded-v1":
        raise ObservationError(f"{summary.name}: unsupported bounded algorithm")
    if integers["max_tracked_unique"] != 100_000:
        raise ObservationError(f"{summary.name}: unsupported bounded-analysis parameters")
    if not 0 < integers["tracked_unique_sequences"] <= integers["max_tracked_unique"]:
        raise ObservationError(f"{summary.name}: tracked_unique_sequences is inconsistent")
    if not 0 < integers["count_at_unique_limit"] <= integers["total_reads"]:
        raise ObservationError(f"{summary.name}: count_at_unique_limit is inconsistent")
    if row["sampling_limited"] not in {"true", "false"}:
        raise ObservationError(f"{summary.name}: sampling_limited must be true or false")
    sampling_limited = row["sampling_limited"] == "true"
    if sampling_limited != (
        integers["tracked_unique_sequences"] == integers["max_tracked_unique"]
    ):
        raise ObservationError(f"{summary.name}: sampling_limited is inconsistent")
    if not sampling_limited and (
        integers["count_at_unique_limit"] != integers["total_reads"]
    ):
        raise ObservationError(f"{summary.name}: unsampled counts are inconsistent")
    return remaining


def _duplication(path: Path, read: str) -> dict[str, float]:
    native_remaining = _native_duplication_summary(path, read)
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
        totals: list[float] = []
        deduplicated_values: list[float] = []
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
            totals.append(total)
            deduplicated_values.append(deduplicated)
            if label in {"1", "1.0"}:
                level_one = (total, deduplicated)

    if native_remaining is not None:
        if not math.isclose(sum(totals), 100.0, abs_tol=1e-5):
            raise ObservationError("native total-sequence percentages must sum to 100")
        if not math.isclose(sum(deduplicated_values), 100.0, abs_tol=1e-5):
            raise ObservationError("native deduplicated percentages must sum to 100")
        if level_one is not None and level_one[1] > 0:
            level_one_remaining = 100.0 * level_one[0] / level_one[1]
            if not math.isclose(native_remaining, level_one_remaining, abs_tol=1e-5):
                raise ObservationError(
                    "native duplication summary and level percentages are inconsistent"
                )
        remaining = native_remaining
    else:
        # Compatibility path for previously imported FastQC two-line profiles.
        if level_one is None:
            raise ObservationError(
                "duplication level 1 or a native summary is required for the headline statistic"
            )
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
