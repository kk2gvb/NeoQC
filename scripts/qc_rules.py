"""Versioned, explainable PASS / WARNING / FAIL evaluation for NeoQC."""

from __future__ import annotations

import hashlib
import json
import math
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

from qc_observations import EXTRACTORS, ObservationError, extract_observations, source_path


DEFAULT_RULESET = (
    Path(__file__).resolve().parents[1]
    / "config"
    / "qc_rules"
    / "fastqc-compatible-v1.json"
)
EVALUATION_FILENAME = "qc_evaluation.json"
STATUSES = ("pass", "warning", "fail", "not_evaluated")
SEVERITY = {"not_evaluated": -1, "pass": 0, "warning": 1, "fail": 2}
OPERATORS = {
    "<": lambda observed, threshold: observed < threshold,
    "<=": lambda observed, threshold: observed <= threshold,
    ">": lambda observed, threshold: observed > threshold,
    ">=": lambda observed, threshold: observed >= threshold,
}


class QcRuleError(ValueError):
    """Raised when rules or evaluation data violate the QC contract."""


@dataclass(frozen=True)
class Threshold:
    operator: str
    value: float

    @classmethod
    def from_value(cls, raw: object, path: str) -> "Threshold | None":
        if raw is None:
            return None
        if not isinstance(raw, Mapping):
            raise QcRuleError(f"{path} must be an object")
        operator = raw.get("operator")
        value = raw.get("value")
        if operator not in OPERATORS:
            raise QcRuleError(f"{path}.operator is invalid")
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise QcRuleError(f"{path}.value must be numeric")
        number = float(value)
        if not math.isfinite(number):
            raise QcRuleError(f"{path}.value must be finite")
        return cls(str(operator), number)

    def matches(self, observed: float) -> bool:
        return OPERATORS[self.operator](observed, self.value)

    def as_dict(self) -> dict[str, object]:
        return {"operator": self.operator, "value": self.value}


@dataclass(frozen=True)
class Check:
    observation: str
    label: str
    unit: str
    warning: Threshold | None
    fail: Threshold | None


@dataclass(frozen=True)
class MetricRule:
    metric_id: str
    title: str
    checks: tuple[Check, ...]


@dataclass(frozen=True)
class Ruleset:
    ruleset_id: str
    version: str
    library_type: str
    description: str
    sha256: str
    rules: tuple[MetricRule, ...]


def _required_text(raw: Mapping[str, object], key: str, path: str) -> str:
    value = raw.get(key)
    if not isinstance(value, str) or not value.strip():
        raise QcRuleError(f"{path}.{key} is required")
    return value.strip()


def load_ruleset(path: Path = DEFAULT_RULESET) -> Ruleset:
    try:
        payload = path.read_bytes()
        raw = json.loads(payload)
    except OSError as error:
        raise QcRuleError(f"cannot read ruleset {path}: {error}") from error
    except json.JSONDecodeError as error:
        raise QcRuleError(f"invalid JSON in ruleset {path}: {error}") from error
    if not isinstance(raw, Mapping) or raw.get("schema_version") != 1:
        raise QcRuleError("unsupported ruleset schema")
    raw_rules = raw.get("rules")
    if not isinstance(raw_rules, Mapping) or not raw_rules:
        raise QcRuleError("ruleset.rules must be a non-empty object")
    rules: list[MetricRule] = []
    for metric_id, raw_rule in raw_rules.items():
        path_prefix = f"ruleset.rules.{metric_id}"
        if metric_id not in EXTRACTORS or not isinstance(raw_rule, Mapping):
            raise QcRuleError(f"{path_prefix} is invalid")
        raw_checks = raw_rule.get("checks")
        if not isinstance(raw_checks, list) or not raw_checks:
            raise QcRuleError(f"{path_prefix}.checks must be a non-empty array")
        checks: list[Check] = []
        for index, raw_check in enumerate(raw_checks):
            check_path = f"{path_prefix}.checks[{index}]"
            if not isinstance(raw_check, Mapping):
                raise QcRuleError(f"{check_path} must be an object")
            warning = Threshold.from_value(raw_check.get("warning"), f"{check_path}.warning")
            fail = Threshold.from_value(raw_check.get("fail"), f"{check_path}.fail")
            if warning is None and fail is None:
                raise QcRuleError(f"{check_path} needs a warning or fail threshold")
            checks.append(
                Check(
                    observation=_required_text(raw_check, "observation", check_path),
                    label=_required_text(raw_check, "label", check_path),
                    unit=str(raw_check.get("unit") or ""),
                    warning=warning,
                    fail=fail,
                )
            )
        rules.append(
            MetricRule(
                metric_id=str(metric_id),
                title=_required_text(raw_rule, "title", path_prefix),
                checks=tuple(checks),
            )
        )
    return Ruleset(
        ruleset_id=_required_text(raw, "id", "ruleset"),
        version=_required_text(raw, "version", "ruleset"),
        library_type=_required_text(raw, "library_type", "ruleset"),
        description=str(raw.get("description") or ""),
        sha256=hashlib.sha256(payload).hexdigest(),
        rules=tuple(rules),
    )


def _number(value: float) -> str:
    return f"{value:.4g}"


def _evaluate_rule(rule: MetricRule, observations: Mapping[str, float]) -> tuple[str, list[dict[str, object]]]:
    status = "pass"
    reasons: list[dict[str, object]] = []
    for check in rule.checks:
        if check.observation not in observations:
            raise QcRuleError(
                f"{rule.metric_id}: observation {check.observation} is missing"
            )
        observed = observations[check.observation]
        check_status = "pass"
        threshold: Threshold | None = None
        if check.fail is not None and check.fail.matches(observed):
            check_status, threshold = "fail", check.fail
        elif check.warning is not None and check.warning.matches(observed):
            check_status, threshold = "warning", check.warning
        if SEVERITY[check_status] > SEVERITY[status]:
            status = check_status
        if threshold is not None:
            unit = f" {check.unit}" if check.unit else ""
            reasons.append(
                {
                    "code": f"{rule.metric_id}.{check.observation}.{check_status}",
                    "message": (
                        f"{check.label}: {_number(observed)}{unit}; "
                        f"{check_status.upper()} threshold {threshold.operator} "
                        f"{_number(threshold.value)}{unit}."
                    ),
                    "observation": check.observation,
                    "observed": observed,
                    "threshold": threshold.value,
                    "operator": threshold.operator,
                    "unit": check.unit,
                }
            )
    if not reasons:
        reasons.append(
            {
                "code": f"{rule.metric_id}.within_thresholds",
                "message": "All evaluated observations are within configured thresholds.",
            }
        )
    return status, reasons


def _not_evaluated(rule: MetricRule, read: str, code: str, message: str) -> dict[str, object]:
    return {
        "metric_id": rule.metric_id,
        "read": read,
        "title": rule.title,
        "qc_status": "not_evaluated",
        "observations": {},
        "checks": [],
        "reasons": [{"code": code, "message": message}],
    }


def _checks(rule: MetricRule) -> list[dict[str, object]]:
    result: list[dict[str, object]] = []
    for check in rule.checks:
        item: dict[str, object] = {
            "observation": check.observation,
            "label": check.label,
            "unit": check.unit,
        }
        if check.warning is not None:
            item["warning"] = check.warning.as_dict()
        if check.fail is not None:
            item["fail"] = check.fail.as_dict()
        result.append(item)
    return result


def _active_reads(input_dir: Path, ruleset: Ruleset) -> tuple[str, ...]:
    return tuple(
        read
        for read in ("R1", "R2")
        if any(source_path(input_dir, rule.metric_id, read).is_file() for rule in ruleset.rules)
    )


def evaluate_directory(input_dir: Path, ruleset_path: Path = DEFAULT_RULESET) -> dict[str, object]:
    input_dir = input_dir.resolve()
    ruleset = load_ruleset(ruleset_path.resolve())
    reads = _active_reads(input_dir, ruleset)
    if not reads:
        raise QcRuleError(f"no recognized NeoQC TSV files found in {input_dir}")
    evaluations: list[dict[str, object]] = []
    for read in reads:
        for rule in ruleset.rules:
            path = source_path(input_dir, rule.metric_id, read)
            if not path.is_file():
                evaluations.append(
                    _not_evaluated(
                        rule,
                        read,
                        "evaluation.source_not_found",
                        f"Source data was not produced: {path.name}.",
                    )
                )
                continue
            try:
                observations = extract_observations(input_dir, rule.metric_id, read)
                status, reasons = _evaluate_rule(rule, observations)
                evaluations.append(
                    {
                        "metric_id": rule.metric_id,
                        "read": read,
                        "title": rule.title,
                        "qc_status": status,
                        "observations": observations,
                        "checks": _checks(rule),
                        "reasons": reasons,
                    }
                )
            except (ObservationError, QcRuleError) as error:
                evaluations.append(
                    _not_evaluated(
                        rule,
                        read,
                        "evaluation.data_invalid",
                        f"QC evaluation unavailable: {error}",
                    )
                )
    counts = {status: 0 for status in STATUSES}
    for evaluation in evaluations:
        counts[str(evaluation["qc_status"])] += 1
    evaluated = [status for status in ("pass", "warning", "fail") if counts[status]]
    overall = max(evaluated, key=SEVERITY.__getitem__) if evaluated else "not_evaluated"
    return {
        "schema_version": 1,
        "ruleset": {
            "id": ruleset.ruleset_id,
            "version": ruleset.version,
            "sha256": ruleset.sha256,
            "library_type": ruleset.library_type,
            "description": ruleset.description,
        },
        "reads": list(reads),
        "summary": {**counts, "overall_status": overall},
        "evaluations": evaluations,
    }


def write_evaluation(
    input_dir: Path,
    output_path: Path | None = None,
    ruleset_path: Path = DEFAULT_RULESET,
) -> Path:
    input_dir = input_dir.resolve()
    destination = (output_path or input_dir / EVALUATION_FILENAME).resolve()
    result = evaluate_directory(input_dir, ruleset_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=destination.parent,
            prefix=f".{destination.name}.",
            suffix=".tmp",
            delete=False,
        ) as temporary:
            temporary_path = temporary.name
            json.dump(result, temporary, indent=2, ensure_ascii=False)
            temporary.write("\n")
            temporary.flush()
            os.fsync(temporary.fileno())
        os.replace(temporary_path, destination)
    except Exception:
        if temporary_path is not None:
            Path(temporary_path).unlink(missing_ok=True)
        raise
    return destination
