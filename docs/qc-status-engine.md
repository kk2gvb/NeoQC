# NeoQC PASS / WARNING / FAIL decision engine

## Purpose

NeoQC needs an explainable, versioned evaluation layer between calculated QC
data and report presentation. A chart must never decide its own status, and a
plot-rendering failure must never be reported as biological QC failure.

The implementation therefore keeps two independent state dimensions:

- `artifact_status`: `generated`, `skipped` or `error`; this remains in
  `plots_manifest.json` and describes report artifact availability;
- `qc_status`: `pass`, `warning`, `fail` or `not_evaluated`; this is written to
  `qc_evaluation.json` and describes the result of a versioned QC rule.

PASS / WARNING / FAIL are technical QC flags, not clinical conclusions. Their
interpretation must take the library type and the downstream workflow into
account.

## Processing flow

```text
NeoQC TSV / summary data
        |
        +--> validated observations --> versioned rule engine
        |                                  |
        |                                  +--> qc_evaluation.json
        |
        +--> plot renderer -------------> plots_manifest.json + SVG / PNG
                                               |
qc_evaluation.json + plots_manifest.json ------+--> NeoQC HTML report
```

`plot_results.py` orchestrates evaluation and plotting, while the evaluator
remains usable as a standalone Python API and CLI. Report renderers
join plot and evaluation records by `(metric_id, read)`.

## Python modules

- `scripts/qc_observations.py`: strict TSV readers and extraction of normalized
  observed values used by both plotting and evaluation;
- `scripts/qc_rules.py`: immutable status enum, rule models, comparisons,
  aggregation and reason generation;
- `scripts/evaluate_qc.py`: standalone CLI which writes `qc_evaluation.json`;
- `config/qc_rules/fastqc-compatible-v1.json`: reproducible baseline rules;
- future `config/qc_rules/neoqc-rna-v1.json`: RNA-seq profile, introduced only
  after its deviations from the baseline have been reviewed and documented.

No plotting function should contain a QC threshold. Graph shading may display
thresholds supplied by an evaluation result, but does not calculate them.

## Evaluation contract

The initial `qc_evaluation.json` contract is:

```json
{
  "schema_version": 1,
  "ruleset": {
    "id": "fastqc-compatible-v1",
    "version": "1.0.0",
    "sha256": "...",
    "library_type": "rna-seq"
  },
  "summary": {
    "pass": 9,
    "warning": 4,
    "fail": 1,
    "not_evaluated": 2,
    "overall_status": "fail"
  },
  "evaluations": [
    {
      "metric_id": "sequence_duplication_levels",
      "read": "R1",
      "qc_status": "warning",
      "observations": {
        "deduplicated_remaining_percent": 68.6,
        "deduplicated_loss_percent": 31.4
      },
      "checks": [{
        "observation": "deduplicated_loss_percent",
        "label": "Sequences removed by deduplication",
        "warning": {"operator": ">", "value": 20.0},
        "fail": {"operator": ">", "value": 50.0},
        "unit": "%"
      }],
      "reasons": [
        {
          "code": "duplication.non_unique_above_warning",
          "message": "Non-unique sequences exceed the warning threshold.",
          "observed": 31.4,
          "threshold": 20.0
        }
      ]
    }
  ]
}
```

Every decision includes its observed value, exact comparison, threshold and a
stable machine-readable reason code. Human-readable messages are escaped by
the renderer and must not be the integration contract.

## Baseline rules and data readiness

The FastQC-compatible profile is a useful reproducible baseline. It does not
automatically define acceptance criteria for an RNA-seq or clinical workflow.

| NeoQC metric | Baseline decision input | Current data readiness |
| --- | --- | --- |
| `per_base_quality` | minimum lower quartile and minimum median by position | Ready; C++ writes both values with the mean |
| `per_sequence_quality` | mode of truncated mean read quality | Ready; legacy two-column TSV files fall back to the rounded distribution |
| `per_base_sequence_content` | maximum `abs(A-T)` or `abs(G-C)` by position | Ready |
| `per_sequence_gc_content` | total deviation from modeled normal distribution | Ready |
| `per_base_n_content` | maximum N percentage at any position | Ready |
| `sequence_length_distribution` | variable lengths and presence of zero-length reads | Ready |
| `sequence_duplication_levels` | percentage lost after exact 50-nt-prefix deduplication | Ready; native runs use exact total/unique counts, legacy imported profiles fall back to the level-1 ratio |
| `adapter_content` | maximum cumulative percentage for any adapter | Ready |

FastQC-compatible thresholds encoded in the versioned configuration are:

| Metric | WARNING | FAIL |
| --- | --- | --- |
| Per-base quality | any lower quartile `< 10` or median `< 25` | any lower quartile `< 5` or median `< 20` |
| Per-sequence quality | modal mean quality `< 27` | modal mean quality `< 20` |
| Base content | `abs(A-T)` or `abs(G-C) > 10%` | difference `> 20%` |
| GC content | modeled deviation `> 15%` | deviation `> 30%` |
| N content | any position `> 5%` | any position `> 20%` |
| Length distribution | more than one observed length | any zero-length sequence |
| Duplication | non-unique sequences `> 20%` | non-unique sequences `> 50%` |
| Adapter content | any adapter `> 5%` | any adapter `> 10%` |

Threshold boundary semantics (`>`, `<`, `>=`, `<=`) are part of the ruleset and
must have dedicated tests.

## RNA-seq interpretation

The report must display the active ruleset. A FastQC-compatible status assumes
a random, diverse library; RNA-seq libraries can have expected initial base
composition bias and biological duplication. The first release should show the
baseline flag with an explanatory RNA-seq note. A separate `neoqc-rna-v1`
profile must not be activated until its thresholds and rationale are approved.

Changing a ruleset never changes historical output silently. The identifier,
semantic version and content hash are stored in every evaluation result.

## Aggregation

Status severity is ordered as `fail > warning > pass`. `not_evaluated` is not a
pass and is excluded from the evaluated denominator.

- a metric module combines R1 and R2 using the worst evaluated status;
- the summary counts each read-level evaluation explicitly;
- overall technical QC is the worst evaluated status;
- missing required data can make a module `not_evaluated`, but never `pass`;
- artifact `error` is displayed independently and does not become QC `fail`;
- technical QC overall must not overwrite processing status or clinical case
  status.

## Представление в HTML-отчёте

HTML-отчёт использует модель оценки так:

1. Блок «Обзор QC» показывает итоговый статус, счётчики `PASS / WARNING / FAIL /
   NOT EVALUATED` и пропорциональную полосу.
2. Матрица «модуль × рид» позволяет сравнить R1 и R2 одним взглядом.
3. В каждой дорожке R1/R2 показаны статус рида, наблюдаемые значения и причина
   статуса с порогом.
4. Статус модуля в заголовке панели и в навигации — худший из статусов ридов.
5. Статус никогда не передаётся только цветом: у плашек и ячеек есть символ
   (✓ ▲ ✕ —) и текст или подсказка.
6. Пороговые зоны рисуются на графике только тогда, когда порог честно
   соответствует оси графика; иначе пояснение остаётся под графиком.
7. Печать и PDF сохраняют плашки, матрицу и причины и не зависят от наведения
   мыши.
8. В режиме «Приоритет» модули `FAIL`, `WARNING` и `NOT EVALUATED` идут первыми,
   модули `PASS` — в конце в свёрнутом виде.

Цвета статусов: зелёный — `PASS`, янтарный — `WARNING`, красный — `FAIL`,
серый — `NOT EVALUATED`; значения для светлой и тёмной темы заданы в
`scripts/qc_report.py`. Интерфейс, режимы и контракт отчёта описаны в
[`html-report.md`](html-report.md).

## Backward compatibility

Reports generated from an older `plots_manifest.json` without
`qc_evaluation.json` remain valid. They show artifact state and
`NOT EVALUATED`, never an inferred PASS. `plots_manifest.json` stays on schema
version 1 because QC decisions are kept in a separate contract.

## Test strategy

- unit tests for every rule and exact threshold boundary;
- malformed, missing, NaN, infinite and out-of-range observation tests;
- R1-only, paired R1/R2 and asymmetric status tests;
- aggregation tests including `not_evaluated` and artifact errors;
- ruleset version/hash reproducibility tests;
- JSON contract and escaping tests;
- HTML report tests (summary counts, matrix, badges, reason text, themes,
  languages, ordering) — see [`html-report.md`](html-report.md);
- comparison fixtures against known FastQC examples before calling the profile
  FastQC-compatible.

## Implementation status

Implemented now:

- immutable rules, validated observations and atomic `qc_evaluation.json`;
- the FastQC-compatible versioned configuration and content hash;
- per-base lower quartile/median output and exact duplication headline;
- exact native prefix-duplication calculation, provenance summary and incomplete-run guard;
- automatic evaluation from `plot_results.py` and standalone evaluation CLI;
- backward-compatible `NOT EVALUATED` behaviour;
- summary strip, R1/R2 matrix, accessible badges, observations and reasons in
  the NeoQC report (see [`html-report.md`](html-report.md));
- boundary, malformed-data, integration, HTML and C++ output contract tests.

The RNA-specific profile remains intentionally pending domain approval and must
not silently replace the baseline ruleset.
