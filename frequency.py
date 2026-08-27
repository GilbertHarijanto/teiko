#!/usr/bin/env python3
"""Display relative frequency of each cell population in each sample.

Reads from the `population_frequencies` view in cell_count.db
(created by load_data.py). Each row is one population from one sample:

    sample, total_count, population, count, percentage

Run `python load_data.py` first if the database does not exist.
"""

from __future__ import annotations

import argparse
import sqlite3
import sys

from load_data import DB_PATH, POPULATIONS

COLUMNS = ("sample", "total_count", "population", "count", "percentage")

FREQUENCY_QUERY = """
SELECT sample, total_count, population, count, percentage
FROM population_frequencies
ORDER BY sample,
    CASE population
        WHEN 'b_cell' THEN 1
        WHEN 'cd8_t_cell' THEN 2
        WHEN 'cd4_t_cell' THEN 3
        WHEN 'nk_cell' THEN 4
        WHEN 'monocyte' THEN 5
    END
"""

PREVIEW_ROWS = 15


def connect() -> sqlite3.Connection:
    if not DB_PATH.exists():
        raise FileNotFoundError(
            f"Database not found: {DB_PATH.name}. Run `python load_data.py` first."
        )
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    return conn


def fetch_frequencies(conn: sqlite3.Connection) -> list[sqlite3.Row]:
    return conn.execute(FREQUENCY_QUERY).fetchall()


def format_row(row: sqlite3.Row) -> str:
    return (
        f"{row['sample']:<14} "
        f"{row['total_count']:>11} "
        f"{row['population']:<12} "
        f"{row['count']:>8} "
        f"{row['percentage']:>11.4f}"
    )


def format_header() -> str:
    return (
        f"{'sample':<14} "
        f"{'total_count':>11} "
        f"{'population':<12} "
        f"{'count':>8} "
        f"{'percentage':>11}"
    )


def display(rows: list[sqlite3.Row], show_all: bool) -> None:
    n_samples = len(rows) // len(POPULATIONS)
    print(
        f"Relative frequency of each cell population in each sample "
        f"({len(rows)} rows = {n_samples} samples × {len(POPULATIONS)} populations)\n"
    )
    print(format_header())
    print("-" * 60)

    visible = rows if show_all else rows[:PREVIEW_ROWS]
    for row in visible:
        print(format_row(row))

    if not show_all and len(rows) > PREVIEW_ROWS:
        remaining = len(rows) - PREVIEW_ROWS
        print(
            f"... {remaining} more rows. "
            f"Re-run with --all to print the full table, or query "
            f"population_frequencies in {DB_PATH.name}."
        )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Display per-sample cell-population relative frequencies."
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Print every row instead of a preview.",
    )
    args = parser.parse_args()

    try:
        conn = connect()
    except FileNotFoundError as exc:
        print(exc, file=sys.stderr)
        return 1

    try:
        rows = fetch_frequencies(conn)
    finally:
        conn.close()

    if not rows:
        print("population_frequencies is empty. Re-run `python load_data.py`.", file=sys.stderr)
        return 1

    display(rows, show_all=args.all)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
