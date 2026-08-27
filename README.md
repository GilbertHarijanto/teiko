# Miraclib immune-cell analysis

Python pipeline and interactive dashboard for immune-cell counts from the miraclib clinical trial (`cell-count.csv`).

## Run (GitHub Codespaces)

The repo is set up to be reproduced with Make. From the repository root:

```bash
make setup      # install dependencies into .venv
make pipeline   # load SQLite, print Part 2–4 tables, write boxplots
make dashboard  # start the Streamlit app on port 8501
```

`make pipeline` is headless (matplotlib uses the Agg backend). It:

1. Creates `cell_count.db` from `cell-count.csv` (Part 1)
2. Prints the relative-frequency summary table (Part 2; full table is the `population_frequencies` view)
3. Prints the responder vs non-responder statistics and writes `figures/response_boxplots.png` (Part 3)
4. Prints the baseline melanoma / PBMC / miraclib subset counts (Part 4)

Individual scripts can also be run with `.venv/bin/python load_data.py` (and the other modules) after `make setup`.

### Dashboard

[http://localhost:8501](http://localhost:8501)

After `make dashboard`, open that URL locally, or in Codespaces use **Ports → 8501 → Open in Browser**.

## Database schema

The CSV is wide and denormalized: subject attributes are repeated on every sample row, and the five populations are columns. The SQLite schema is third-normal-form so each fact is stored once.

| Relation | Grain | Contents |
|---|---|---|
| `subjects` | one row per patient | `project`, `condition`, `age`, `sex`, `treatment`, `response` |
| `samples` | one row per biological sample | `sample_type`, `time_from_treatment_start`, FK → `subjects` |
| `cell_counts` | one row per (sample, population) | `population`, `count`, FK → `samples` |
| `population_frequencies` | view | Part 2 table: `sample`, `total_count`, `population`, `count`, `percentage` |

Empty `response` values (healthy / untreated) are stored as `NULL`. Foreign keys, check constraints, and indexes on `samples(subject)`, `samples(sample_type, time_from_treatment_start)`, and `subjects(condition, treatment)` match the Part 3–4 filters.

**Why this shape.** Subject-level fields do not change across timepoints, so they belong on `subjects`, not copied onto every sample. Counts are stored long rather than wide so adding a population is a new row, not a schema change, and relative-frequency queries are a `SUM` plus a join (the view) instead of five column-specific expressions.

**How this scales.** Hundreds of projects: introduce a `projects` table (`project` PK plus study metadata) and keep `subjects.project` as a foreign key. Thousands to millions of samples: the current indexes already serve cohort filters; if several analysts query concurrently, move the same three tables to PostgreSQL without changing the logical model. Other analytics (new assays, more timepoints, per-marker intensities) extend `cell_counts` with an `assay` or `feature` column, or add a `measurements` table with the same `(sample, feature, value)` grain. Derived metrics stay in views or materialized tables so the raw counts remain the source of truth.

## Code structure

| File | Role |
|---|---|
| `load_data.py` | Part 1: schema + CSV load; stdlib only |
| `frequency.py` | Part 2: query `population_frequencies` and print the summary table |
| `response_analysis.py` | Part 3: melanoma / miraclib / PBMC comparison, Mann–Whitney U, boxplots |
| `subset_analysis.py` | Part 4: SQL subset of baseline melanoma PBMC miraclib samples |
| `dashboard.py` | Streamlit UI over the same queries and statistics |
| `requirements.txt` | pandas, scipy, matplotlib, streamlit, plotly |
| `Makefile` | `setup`, `pipeline`, `dashboard` |

Each part is a script that can run on its own (`python frequency.py`, etc.). The dashboard imports those functions instead of reimplementing SQL, so CLI output and the UI cannot drift. `load_data.py` and `frequency.py` stay on the standard library so the database and frequency table do not depend on the analysis stack.
