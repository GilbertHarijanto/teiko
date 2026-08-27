#!/usr/bin/env python3
"""Compare immune-cell relative frequencies in miraclib responders vs non-responders.

Cohort (Part 3): melanoma patients treated with miraclib; PBMC samples only.
Outcome: response yes vs no.

Each subject has three PBMC samples (days 0, 7, and 14). Treating those as
independent would inflate sample size and understate p-values. The primary
analysis therefore uses the per-subject mean relative frequency, then a
two-sided Mann–Whitney U test per population, with Bonferroni correction
across the five populations.

Run `python load_data.py` first if the database does not exist.
"""

from __future__ import annotations

import argparse
import sqlite3
import sys
from pathlib import Path

from frequency import connect
from load_data import POPULATIONS, ROOT

try:
    import matplotlib.pyplot as plt
    import pandas as pd
    from scipy import stats
except ImportError:
    print(
        "Missing analysis packages. From the repo root run:\n"
        "  python3 -m venv .venv && source .venv/bin/activate\n"
        "  pip install -r requirements.txt",
        file=sys.stderr,
    )
    raise SystemExit(1)

ALPHA = 0.05
FIGURE_PATH = ROOT / "figures" / "response_boxplots.png"

POPULATION_LABELS = {
    "b_cell": "B cell",
    "cd8_t_cell": "CD8 T cell",
    "cd4_t_cell": "CD4 T cell",
    "nk_cell": "NK cell",
    "monocyte": "Monocyte",
}

COHORT_SQL = """
SELECT
    s.subject,
    s.response,
    sm.time_from_treatment_start,
    f.sample,
    f.population,
    f.percentage
FROM population_frequencies AS f
JOIN samples AS sm ON sm.sample = f.sample
JOIN subjects AS s ON s.subject = sm.subject
WHERE s.condition = 'melanoma'
  AND s.treatment = 'miraclib'
  AND sm.sample_type = 'PBMC'
  AND s.response IN ('yes', 'no')
"""


def fetch_cohort(conn: sqlite3.Connection) -> pd.DataFrame:
    df = pd.read_sql_query(COHORT_SQL, conn)
    if df.empty:
        raise RuntimeError("Cohort query returned no rows. Re-run `python load_data.py`.")
    df["response"] = pd.Categorical(df["response"], categories=["no", "yes"], ordered=True)
    df["population"] = pd.Categorical(df["population"], categories=list(POPULATIONS), ordered=True)
    return df


def subject_means(cohort: pd.DataFrame) -> pd.DataFrame:
    """One independent observation per subject per population."""
    return cohort.groupby(
        ["subject", "response", "population"], observed=True, as_index=False
    )["percentage"].mean()


def _rank_biserial(u: float, n_yes: int, n_no: int) -> float:
    """Positive values mean responders tend to have higher relative frequency."""
    return (2.0 * u) / (n_yes * n_no) - 1.0


def test_populations(means: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for population in POPULATIONS:
        sub = means.loc[means["population"] == population]
        yes = sub.loc[sub["response"] == "yes", "percentage"]
        no = sub.loc[sub["response"] == "no", "percentage"]
        u, p_raw = stats.mannwhitneyu(yes, no, alternative="two-sided")
        rows.append(
            {
                "population": population,
                "n_responder": int(yes.size),
                "n_nonresponder": int(no.size),
                "median_responder": float(yes.median()),
                "median_nonresponder": float(no.median()),
                "iqr_responder": float(yes.quantile(0.75) - yes.quantile(0.25)),
                "iqr_nonresponder": float(no.quantile(0.75) - no.quantile(0.25)),
                "U": float(u),
                "p_raw": float(p_raw),
                "rank_biserial": _rank_biserial(float(u), int(yes.size), int(no.size)),
            }
        )

    results = pd.DataFrame(rows)
    results["p_bonferroni"] = (results["p_raw"] * len(POPULATIONS)).clip(upper=1.0)
    results["significant"] = results["p_bonferroni"] < ALPHA
    return results


def plot_boxplots(means: pd.DataFrame, results: pd.DataFrame, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    n_yes = int((means.drop_duplicates("subject")["response"] == "yes").sum())
    n_no = int((means.drop_duplicates("subject")["response"] == "no").sum())

    fig, axes = plt.subplots(1, len(POPULATIONS), figsize=(16, 4.8), layout="constrained")
    fig.suptitle(
        "Melanoma · miraclib · PBMC\n"
        "Subject-mean relative frequency, responders vs non-responders",
        fontsize=13,
    )

    box_colors = ("#E76F51", "#2A9D8F")
    for ax, population in zip(axes, POPULATIONS):
        sub = means.loc[means["population"] == population]
        no_vals = sub.loc[sub["response"] == "no", "percentage"].to_numpy()
        yes_vals = sub.loc[sub["response"] == "yes", "percentage"].to_numpy()
        bp = ax.boxplot(
            [no_vals, yes_vals],
            tick_labels=[f"Non-resp.\n(n={n_no})", f"Resp.\n(n={n_yes})"],
            patch_artist=True,
            widths=0.6,
            medianprops={"color": "#1d3557", "linewidth": 2},
            whiskerprops={"color": "#1d3557"},
            capprops={"color": "#1d3557"},
            flierprops={"marker": "o", "markersize": 3, "alpha": 0.35, "color": "#1d3557"},
            boxprops={"linewidth": 1.2, "edgecolor": "#1d3557"},
        )
        for patch, color in zip(bp["boxes"], box_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.8)

        p_adj = float(results.loc[results["population"] == population, "p_bonferroni"].iloc[0])
        p_raw = float(results.loc[results["population"] == population, "p_raw"].iloc[0])
        ax.set_title(f"{POPULATION_LABELS[population]}\n$p_{{adj}}$={p_adj:.3f}  ($p$={p_raw:.3f})")
        ax.set_ylim(bottom=0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[0].set_ylabel("Relative frequency (%)")
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def _fmt_p(value: float) -> str:
    if value < 1e-4:
        return f"{value:.2e}"
    return f"{value:.4f}"


def print_report(cohort: pd.DataFrame, results: pd.DataFrame, figure_path: Path) -> None:
    n_subjects = cohort["subject"].nunique()
    n_samples = cohort["sample"].nunique()
    n_yes = cohort.loc[cohort["response"] == "yes", "subject"].nunique()
    n_no = cohort.loc[cohort["response"] == "no", "subject"].nunique()

    print("Part 3: miraclib response vs immune-cell relative frequency")
    print("=" * 78)
    print("Cohort:      melanoma, treatment=miraclib, sample_type=PBMC")
    print(f"Subjects:    {n_subjects}  ({n_yes} responders, {n_no} non-responders)")
    print(f"Samples:     {n_samples}  (days 0, 7, 14 for each subject)")
    print()
    print("Methods")
    print("-------")
    print("Frequencies come from the population_frequencies summary table (Part 2).")
    print("Because each subject contributes three samples, tests use the per-subject")
    print("mean relative frequency (one independent observation per person).")
    print("Test: two-sided Mann–Whitney U (Wilcoxon rank-sum), α=0.05.")
    print(f"Multiple testing: Bonferroni correction across {len(POPULATIONS)} populations")
    print(f"(threshold p < {ALPHA / len(POPULATIONS):.3f} unadjusted, or p_adj < {ALPHA}).")
    print("Effect size: rank-biserial correlation r (positive ⇒ higher in responders).")
    print()
    print(
        f"{'population':<12} {'n_yes':>6} {'n_no':>6} "
        f"{'med_yes':>8} {'med_no':>8} {'U':>9} "
        f"{'p':>8} {'p_adj':>8} {'r_rb':>7} {'sig':>4}"
    )
    print("-" * 78)
    for row in results.itertuples(index=False):
        sig = "yes" if row.significant else "no"
        print(
            f"{row.population:<12} {row.n_responder:>6} {row.n_nonresponder:>6} "
            f"{row.median_responder:>8.3f} {row.median_nonresponder:>8.3f} {row.U:>9.1f} "
            f"{_fmt_p(row.p_raw):>8} {_fmt_p(row.p_bonferroni):>8} "
            f"{row.rank_biserial:>7.3f} {sig:>4}"
        )
    print()

    significant = results.loc[results["significant"], "population"].tolist()
    print("Conclusion")
    print("----------")
    if significant:
        labels = ", ".join(POPULATION_LABELS[p] for p in significant)
        print(
            f"After Bonferroni correction, {labels} "
            f"{'has' if len(significant) == 1 else 'have'} a significant difference "
            "in relative frequency between responders and non-responders."
        )
    else:
        closest = results.loc[results["p_raw"].idxmin()]
        print(
            "After Bonferroni correction, no cell population has a statistically"
        )
        print("significant difference in relative frequency between responders and")
        print("non-responders. There is not enough evidence in this cohort to treat")
        print("any of these five frequencies as a predictor of miraclib response.")
        print()
        print(
            f"Closest result: {POPULATION_LABELS[closest.population]} "
            f"(raw p={_fmt_p(closest.p_raw)}, adjusted p={_fmt_p(closest.p_bonferroni)}, "
            f"r={closest.rank_biserial:.3f})."
        )
        print("That effect is small and does not survive multiple-testing correction.")
    print()
    print(f"Boxplots written to {figure_path}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare cell-population frequencies in miraclib responders vs non-responders."
    )
    parser.add_argument(
        "--figure",
        type=Path,
        default=FIGURE_PATH,
        help=f"Where to save the boxplot figure (default: {FIGURE_PATH})",
    )
    args = parser.parse_args()

    try:
        conn = connect()
    except FileNotFoundError as exc:
        print(exc, file=sys.stderr)
        return 1

    try:
        cohort = fetch_cohort(conn)
    finally:
        conn.close()

    means = subject_means(cohort)
    results = test_populations(means)
    figure_path = plot_boxplots(means, results, args.figure)
    print_report(cohort, results, figure_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
