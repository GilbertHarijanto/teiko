#!/usr/bin/env python3
"""Subset melanoma PBMC baseline samples from miraclib-treated patients.

Queries cell_count.db (created by load_data.py) for:

  condition = melanoma
  treatment = miraclib
  sample_type = PBMC
  time_from_treatment_start = 0

Then reports sample counts by project and subject counts by response and sex.
"""

from __future__ import annotations

import argparse
import sqlite3
import sys

from frequency import connect

COHORT_SQL = """
SELECT
    sm.sample,
    s.subject,
    s.project,
    s.response,
    s.sex
FROM samples AS sm
JOIN subjects AS s ON s.subject = sm.subject
WHERE s.condition = 'melanoma'
  AND s.treatment = 'miraclib'
  AND sm.sample_type = 'PBMC'
  AND sm.time_from_treatment_start = 0
ORDER BY sm.sample
"""

SAMPLES_BY_PROJECT_SQL = """
SELECT
    p.project,
    COUNT(c.sample) AS n_samples
FROM (SELECT DISTINCT project FROM subjects) AS p
LEFT JOIN (
    SELECT s.project, sm.sample
    FROM samples AS sm
    JOIN subjects AS s ON s.subject = sm.subject
    WHERE s.condition = 'melanoma'
      AND s.treatment = 'miraclib'
      AND sm.sample_type = 'PBMC'
      AND sm.time_from_treatment_start = 0
) AS c ON c.project = p.project
GROUP BY p.project
ORDER BY p.project
"""

SUBJECTS_BY_RESPONSE_SQL = """
SELECT
    s.response,
    COUNT(DISTINCT s.subject) AS n_subjects
FROM samples AS sm
JOIN subjects AS s ON s.subject = sm.subject
WHERE s.condition = 'melanoma'
  AND s.treatment = 'miraclib'
  AND sm.sample_type = 'PBMC'
  AND sm.time_from_treatment_start = 0
GROUP BY s.response
ORDER BY s.response
"""

SUBJECTS_BY_SEX_SQL = """
SELECT
    s.sex,
    COUNT(DISTINCT s.subject) AS n_subjects
FROM samples AS sm
JOIN subjects AS s ON s.subject = sm.subject
WHERE s.condition = 'melanoma'
  AND s.treatment = 'miraclib'
  AND sm.sample_type = 'PBMC'
  AND sm.time_from_treatment_start = 0
GROUP BY s.sex
ORDER BY s.sex
"""

PREVIEW_ROWS = 8
RESPONSE_LABELS = {"yes": "responder", "no": "non-responder"}
SEX_LABELS = {"M": "male", "F": "female"}


def fetch_cohort(conn: sqlite3.Connection) -> list[sqlite3.Row]:
    return conn.execute(COHORT_SQL).fetchall()


def fetch_counts(
    conn: sqlite3.Connection,
) -> tuple[list[sqlite3.Row], list[sqlite3.Row], list[sqlite3.Row]]:
    return (
        conn.execute(SAMPLES_BY_PROJECT_SQL).fetchall(),
        conn.execute(SUBJECTS_BY_RESPONSE_SQL).fetchall(),
        conn.execute(SUBJECTS_BY_SEX_SQL).fetchall(),
    )


def print_section(title: str, rows: list[tuple[str, int]], total_label: str) -> None:
    print(title)
    print("-" * 40)
    total = 0
    width = max(len(label) for label, _ in rows)
    for label, count in rows:
        print(f"  {label:<{width}}  {count:>5}")
        total += count
    print(f"  {'total':<{width}}  {total:>5}  {total_label}")
    print()


def display(
    cohort: list[sqlite3.Row],
    by_project: list[sqlite3.Row],
    by_response: list[sqlite3.Row],
    by_sex: list[sqlite3.Row],
    show_all: bool,
) -> None:
    n_samples = len(cohort)
    n_subjects = len({row["subject"] for row in cohort})

    print("Part 4: Baseline melanoma PBMC miraclib subset")
    print("=" * 50)
    print("Filters:  condition=melanoma, treatment=miraclib,")
    print("          sample_type=PBMC, time_from_treatment_start=0")
    print(f"Samples:  {n_samples}  ({n_subjects} subjects; one baseline sample each)")
    print()

    print("Identified samples")
    print("-" * 40)
    print(f"  {'sample':<14} {'subject':<10} {'project':<8} {'response':<10} {'sex':<4}")
    visible = cohort if show_all else cohort[:PREVIEW_ROWS]
    for row in visible:
        print(
            f"  {row['sample']:<14} {row['subject']:<10} {row['project']:<8} "
            f"{row['response']:<10} {row['sex']:<4}"
        )
    if not show_all and n_samples > PREVIEW_ROWS:
        print(f"  ... {n_samples - PREVIEW_ROWS} more. Re-run with --all to list every sample.")
    print()

    print_section(
        "Samples by project",
        [(row["project"], int(row["n_samples"])) for row in by_project],
        "samples",
    )
    print_section(
        "Subjects by response",
        [
            (f"{row['response']} ({RESPONSE_LABELS.get(row['response'], row['response'])})", int(row["n_subjects"]))
            for row in by_response
        ],
        "subjects",
    )
    print_section(
        "Subjects by sex",
        [
            (f"{row['sex']} ({SEX_LABELS.get(row['sex'], row['sex'])})", int(row["n_subjects"]))
            for row in by_sex
        ],
        "subjects",
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Query baseline melanoma PBMC samples from miraclib-treated patients."
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="List every sample in the subset instead of a preview.",
    )
    args = parser.parse_args()

    try:
        conn = connect()
    except FileNotFoundError as exc:
        print(exc, file=sys.stderr)
        return 1

    try:
        cohort = fetch_cohort(conn)
        by_project, by_response, by_sex = fetch_counts(conn)
    finally:
        conn.close()

    if not cohort:
        print("Subset query returned no rows. Re-run `python load_data.py`.", file=sys.stderr)
        return 1

    display(cohort, by_project, by_response, by_sex, show_all=args.all)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
