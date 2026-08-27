#!/usr/bin/env python3
"""Initialize a SQLite database from cell-count.csv.

Schema (3NF):
  subjects                 — one row per patient; attributes that do not vary by sample
  samples                  — one row per biological sample (CSV grain)
  cell_counts              — one row per (sample, population); long form of the five count columns
  population_frequencies   — view: relative frequency of each population within each sample

Empty `response` values (healthy / untreated subjects) are stored as NULL.
Re-running this script replaces the database.
"""

from __future__ import annotations

import csv
import sqlite3
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "cell-count.csv"
DB_PATH = ROOT / "cell_count.db"

POPULATIONS = ("b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte")

SCHEMA_SQL = """
PRAGMA foreign_keys = ON;

CREATE TABLE subjects (
    subject   TEXT PRIMARY KEY,
    project   TEXT NOT NULL,
    condition TEXT NOT NULL CHECK (condition IN ('melanoma', 'carcinoma', 'healthy')),
    age       INTEGER NOT NULL CHECK (age >= 0),
    sex       TEXT NOT NULL CHECK (sex IN ('M', 'F')),
    treatment TEXT NOT NULL CHECK (treatment IN ('miraclib', 'phauximab', 'none')),
    response  TEXT CHECK (response IN ('yes', 'no'))
);

CREATE TABLE samples (
    sample                     TEXT PRIMARY KEY,
    subject                    TEXT NOT NULL REFERENCES subjects (subject),
    sample_type                TEXT NOT NULL CHECK (sample_type IN ('PBMC', 'WB')),
    time_from_treatment_start  INTEGER NOT NULL CHECK (time_from_treatment_start IN (0, 7, 14)),
    UNIQUE (subject, time_from_treatment_start)
);

CREATE TABLE cell_counts (
    sample     TEXT NOT NULL REFERENCES samples (sample),
    population TEXT NOT NULL CHECK (
        population IN ('b_cell', 'cd8_t_cell', 'cd4_t_cell', 'nk_cell', 'monocyte')
    ),
    count      INTEGER NOT NULL CHECK (count >= 0),
    PRIMARY KEY (sample, population)
);

CREATE INDEX idx_samples_subject ON samples (subject);
CREATE INDEX idx_samples_type_time ON samples (sample_type, time_from_treatment_start);
CREATE INDEX idx_subjects_condition_treatment ON subjects (condition, treatment);

CREATE VIEW population_frequencies AS
SELECT
    c.sample AS sample,
    totals.total_count AS total_count,
    c.population AS population,
    c.count AS count,
    100.0 * c.count / totals.total_count AS percentage
FROM cell_counts AS c
JOIN (
    SELECT sample, SUM(count) AS total_count
    FROM cell_counts
    GROUP BY sample
) AS totals ON totals.sample = c.sample;
"""


def _response(value: str) -> str | None:
    value = value.strip()
    return value if value else None


def read_csv(path: Path) -> tuple[list[tuple], list[tuple], list[tuple]]:
    subjects: dict[str, tuple] = {}
    samples: list[tuple] = []
    cell_counts: list[tuple] = []

    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            subject_id = row["subject"]
            attrs = (
                subject_id,
                row["project"],
                row["condition"],
                int(row["age"]),
                row["sex"],
                row["treatment"],
                _response(row["response"]),
            )
            existing = subjects.get(subject_id)
            if existing is None:
                subjects[subject_id] = attrs
            elif existing != attrs:
                raise ValueError(f"Inconsistent subject attributes for {subject_id}")

            sample_id = row["sample"]
            samples.append(
                (
                    sample_id,
                    subject_id,
                    row["sample_type"],
                    int(row["time_from_treatment_start"]),
                )
            )
            for population in POPULATIONS:
                cell_counts.append((sample_id, population, int(row[population])))

    return list(subjects.values()), samples, cell_counts


def create_database(db_path: Path) -> sqlite3.Connection:
    if db_path.exists():
        db_path.unlink()
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys = ON")
    conn.executescript(SCHEMA_SQL)
    return conn


def load(conn: sqlite3.Connection, subjects: list[tuple], samples: list[tuple], cell_counts: list[tuple]) -> None:
    conn.executemany(
        """
        INSERT INTO subjects (subject, project, condition, age, sex, treatment, response)
        VALUES (?, ?, ?, ?, ?, ?, ?)
        """,
        subjects,
    )
    conn.executemany(
        """
        INSERT INTO samples (sample, subject, sample_type, time_from_treatment_start)
        VALUES (?, ?, ?, ?)
        """,
        samples,
    )
    conn.executemany(
        """
        INSERT INTO cell_counts (sample, population, count)
        VALUES (?, ?, ?)
        """,
        cell_counts,
    )
    conn.commit()


def verify(conn: sqlite3.Connection, n_subjects: int, n_samples: int, n_counts: int) -> None:
    counts = {
        table: conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
        for table in ("subjects", "samples", "cell_counts")
    }
    expected = {"subjects": n_subjects, "samples": n_samples, "cell_counts": n_counts}
    if counts != expected:
        raise RuntimeError(f"Row-count mismatch: got {counts}, expected {expected}")

    orphan_samples = conn.execute(
        """
        SELECT COUNT(*) FROM samples s
        LEFT JOIN subjects u ON u.subject = s.subject
        WHERE u.subject IS NULL
        """
    ).fetchone()[0]
    orphan_counts = conn.execute(
        """
        SELECT COUNT(*) FROM cell_counts c
        LEFT JOIN samples s ON s.sample = c.sample
        WHERE s.sample IS NULL
        """
    ).fetchone()[0]
    if orphan_samples or orphan_counts:
        raise RuntimeError(
            f"Broken foreign keys: {orphan_samples} samples, {orphan_counts} cell_counts"
        )

    n_frequencies = conn.execute("SELECT COUNT(*) FROM population_frequencies").fetchone()[0]
    if n_frequencies != n_counts:
        raise RuntimeError(
            f"population_frequencies has {n_frequencies} rows, expected {n_counts}"
        )


def main() -> int:
    if not CSV_PATH.exists():
        print(f"CSV not found: {CSV_PATH}", file=sys.stderr)
        return 1

    subjects, samples, cell_counts = read_csv(CSV_PATH)
    conn = create_database(DB_PATH)
    try:
        load(conn, subjects, samples, cell_counts)
        verify(conn, len(subjects), len(samples), len(cell_counts))
    finally:
        conn.close()

    print(
        f"Wrote {DB_PATH.name}: "
        f"{len(subjects)} subjects, {len(samples)} samples, {len(cell_counts)} cell counts"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
