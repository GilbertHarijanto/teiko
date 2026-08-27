#!/usr/bin/env python3
"""Interactive dashboard for the miraclib immune-cell analysis.

    streamlit run dashboard.py
"""

from __future__ import annotations

import math
import sqlite3

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from frequency import FREQUENCY_QUERY, connect
from load_data import DB_PATH, POPULATIONS
from response_analysis import (
    ALPHA,
    POPULATION_LABELS,
    fetch_cohort,
    subject_means,
    test_populations,
)
from subset_analysis import (
    RESPONSE_LABELS,
    SEX_LABELS,
    fetch_cohort as fetch_baseline,
    fetch_counts,
)

RESPONSE_COLORS = {"Non-responder": "#E76F51", "Responder": "#2A9D8F"}


def _ensure_database() -> None:
    if DB_PATH.exists():
        return
    st.warning("SQLite database not found. Load `cell-count.csv` first.")
    if st.button("Create database now"):
        from load_data import main as load_main

        with st.spinner("Loading cell-count.csv into SQLite…"):
            if load_main() == 0:
                st.success(f"Created {DB_PATH.name}")
                st.rerun()
            st.error("Database load failed.")
    st.stop()


@st.cache_data(show_spinner="Loading frequencies…")
def load_frequencies() -> pd.DataFrame:
    conn = sqlite3.connect(DB_PATH)
    try:
        df = pd.read_sql_query(FREQUENCY_QUERY, conn)
    finally:
        conn.close()
    return df


@st.cache_data(show_spinner="Running response analysis…")
def load_response_analysis() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    conn = connect()
    try:
        cohort = fetch_cohort(conn)
    finally:
        conn.close()
    means = subject_means(cohort)
    results = test_populations(means)
    return cohort, means, results


@st.cache_data(show_spinner="Querying baseline subset…")
def load_baseline_subset() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    conn = connect()
    try:
        cohort_rows = fetch_baseline(conn)
        by_project, by_response, by_sex = fetch_counts(conn)
    finally:
        conn.close()
    cohort = pd.DataFrame([dict(row) for row in cohort_rows])
    return (
        cohort,
        pd.DataFrame([dict(row) for row in by_project]),
        pd.DataFrame([dict(row) for row in by_response]),
        pd.DataFrame([dict(row) for row in by_sex]),
    )


def _fmt_p(value: float) -> str:
    if value < 1e-4:
        return f"{value:.1e}"
    return f"{value:.4g}"


def _fmt_sig(value: float, sig: int = 4) -> str:
    if abs(value) >= 0.9995 and abs(value) < 1.0005:
        return "1"
    if abs(value) < 1e-3:
        return f"{value:.1e}"
    return f"{value:.{sig}g}"


def _sort_frequencies(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["population"] = pd.Categorical(out["population"], categories=list(POPULATIONS), ordered=True)
    return out.sort_values(["sample", "population"], kind="mergesort").reset_index(drop=True)


def _response_legend() -> None:
    st.markdown(
        """
        <div style="display:flex;gap:1.5rem;align-items:center;margin:0.25rem 0 0.75rem 0;">
          <span style="display:flex;align-items:center;gap:0.45rem;font-size:0.95rem;">
            <span style="width:14px;height:14px;background:#E76F51;border-radius:3px;display:inline-block;"></span>
            Non-responder
          </span>
          <span style="display:flex;align-items:center;gap:0.45rem;font-size:0.95rem;">
            <span style="width:14px;height:14px;background:#2A9D8F;border-radius:3px;display:inline-block;"></span>
            Responder
          </span>
        </div>
        """,
        unsafe_allow_html=True,
    )


def _project_bar_figure(by_project: pd.DataFrame) -> go.Figure:
    ymax = max(int(by_project["n_samples"].max()), 1)
    stub = ymax * 0.04
    heights = [stub if int(n) == 0 else int(n) for n in by_project["n_samples"]]
    colors = ["#C5D4CF" if int(n) == 0 else "#2A9D8F" for n in by_project["n_samples"]]
    outlines = ["#8A9A94" if int(n) == 0 else "#1F7A6E" for n in by_project["n_samples"]]
    patterns = ["/" if int(n) == 0 else "" for n in by_project["n_samples"]]
    texts = ["0 · WB only" if int(n) == 0 else f"{int(n)}" for n in by_project["n_samples"]]
    fig = go.Figure(
        go.Bar(
            x=by_project["project"],
            y=heights,
            text=texts,
            textposition="outside",
            cliponaxis=False,
            marker=dict(
                color=colors,
                line=dict(color=outlines, width=1.2),
                pattern=dict(shape=patterns, fgcolor="#6B7C76", size=6, solidity=0.35),
            ),
            customdata=by_project["n_samples"],
            hovertemplate="%{x}: %{customdata} samples<extra></extra>",
        )
    )
    fig.update_layout(
        title="Samples by project",
        yaxis_title="Samples",
        xaxis_title="",
        yaxis_range=[0, ymax * 1.25],
        margin=dict(t=50, b=40),
        height=320,
        showlegend=False,
    )
    return fig


def render_overview(n_freq: int, cohort: pd.DataFrame, results: pd.DataFrame, baseline: pd.DataFrame) -> None:
    n_samples = n_freq // len(POPULATIONS)
    n_subjects = cohort["subject"].nunique()
    n_yes = cohort.loc[cohort["response"] == "yes", "subject"].nunique()
    n_no = cohort.loc[cohort["response"] == "no", "subject"].nunique()

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Samples", f"{n_samples:,}")
    c2.metric("Frequency rows", f"{n_freq:,}")
    c3.metric("Response cohort", f"{n_subjects:,} subjects")
    c4.metric("Baseline subset", f"{len(baseline):,} samples")

    st.subheader("Headline findings")
    significant = results.loc[results["significant"], "population"].tolist()
    if significant:
        labels = ", ".join(POPULATION_LABELS[p] for p in significant)
        st.success(
            f"After Bonferroni correction, {labels} "
            f"{'differs' if len(significant) == 1 else 'differ'} significantly "
            "between miraclib responders and non-responders."
        )
    else:
        closest = results.loc[results["p_raw"].idxmin()]
        st.info(
            "No cell population differs significantly between miraclib responders "
            "and non-responders after Bonferroni correction. Closest result: "
            f"**{POPULATION_LABELS[closest.population]}** "
            f"(raw *p* = {_fmt_p(float(closest.p_raw))}, "
            f"adjusted *p* = {_fmt_p(float(closest.p_bonferroni))})."
        )

    st.caption(
        f"Response analysis cohort: {n_yes} responders and {n_no} non-responders "
        "(melanoma, miraclib, PBMC). Baseline subset is the same cohort at day 0."
    )


def render_frequencies(freq: pd.DataFrame) -> None:
    st.subheader("Relative frequency of each cell type in each sample")
    st.caption(
        "For each sample, total cell count is the sum of the five populations. "
        "Percentage is that population’s share of the sample total."
    )

    left, right = st.columns([2, 3])
    with left:
        sample_query = st.text_input("Filter by sample ID", placeholder="e.g. sample00000")
    with right:
        populations = st.multiselect(
            "Populations",
            options=list(POPULATIONS),
            default=list(POPULATIONS),
        )

    view = freq
    if sample_query.strip():
        view = view[view["sample"].str.contains(sample_query.strip(), case=False, regex=False)]
    if populations:
        view = view[view["population"].isin(populations)]
    else:
        st.warning("Select at least one population.")
        return

    view = _sort_frequencies(view)
    filter_key = (sample_query.strip(), tuple(populations))
    if st.session_state.get("_freq_filter") != filter_key:
        st.session_state["_freq_filter"] = filter_key
        st.session_state["freq_page"] = 1

    pager, sizer = st.columns([2, 1])
    with sizer:
        page_size = st.selectbox("Rows per page", [25, 50], index=0)
    n_pages = max(1, math.ceil(len(view) / page_size))
    if st.session_state.get("freq_page", 1) > n_pages:
        st.session_state["freq_page"] = 1
    with pager:
        page = st.number_input("Page", min_value=1, max_value=n_pages, step=1, key="freq_page")

    start = (int(page) - 1) * int(page_size)
    page_df = view.iloc[start : start + int(page_size)]

    st.dataframe(
        page_df,
        width="stretch",
        hide_index=True,
        column_config={
            "sample": st.column_config.TextColumn("sample"),
            "total_count": st.column_config.NumberColumn("total_count", format="%d"),
            "population": st.column_config.TextColumn("population"),
            "count": st.column_config.NumberColumn("count", format="%d"),
            "percentage": st.column_config.NumberColumn(
                "percentage",
                format="%.2f%%",
                help="Rounded to 2 decimals in the table. Download CSV for full precision.",
            ),
        },
    )
    st.caption(
        f"Rows {start + 1:,}–{start + len(page_df):,} of {len(view):,} matching "
        f"({len(freq):,} total). Sorted by sample, then population. "
        "Percentages shown to 2 decimals; CSV keeps full precision."
    )
    st.download_button(
        "Download matching rows (CSV)",
        view.to_csv(index=False),
        file_name="population_frequencies.csv",
        mime="text/csv",
        help="Full-precision percentages for the current filters (all matching rows, not just this page).",
    )


def render_response(cohort: pd.DataFrame, means: pd.DataFrame, results: pd.DataFrame) -> None:
    n_yes = int((means.drop_duplicates("subject")["response"] == "yes").sum())
    n_no = int((means.drop_duplicates("subject")["response"] == "no").sum())

    st.subheader("Responders vs non-responders")
    st.caption(
        "Melanoma patients treated with miraclib; PBMC samples only. "
        "Each subject has samples at days 0, 7, and 14. Tests use the per-subject "
        "mean relative frequency so observations are independent."
    )

    m1, m2, m3 = st.columns(3)
    m1.metric("Subjects", f"{n_yes + n_no}")
    m2.metric("Responders", f"{n_yes}")
    m3.metric("Non-responders", f"{n_no}")

    with st.expander("Statistical methods"):
        st.markdown(
            f"""
- Source: `population_frequencies` summary table (Part 2)
- Test: two-sided Mann–Whitney U (Wilcoxon rank-sum)
- Multiple testing: Bonferroni across {len(POPULATIONS)} populations
  (α = {ALPHA}, significant if adjusted *p* < {ALPHA})
- Effect size: rank-biserial *r* (positive ⇒ higher frequency in responders)
            """
        )

    show_points = st.checkbox("Show every subject on the boxplots", value=False)
    _response_legend()
    plot_df = means.copy()
    plot_df["Response"] = plot_df["response"].map({"yes": "Responder", "no": "Non-responder"})
    plot_df["Population"] = pd.Categorical(
        plot_df["population"].map(POPULATION_LABELS),
        categories=[POPULATION_LABELS[p] for p in POPULATIONS],
        ordered=True,
    )
    fig = px.box(
        plot_df,
        x="Response",
        y="percentage",
        color="Response",
        facet_col="Population",
        points="all" if show_points else "outliers",
        color_discrete_map=RESPONSE_COLORS,
        category_orders={"Response": ["Non-responder", "Responder"]},
        labels={"percentage": "Relative frequency (%)", "Response": ""},
    )
    fig.update_yaxes(matches=None, rangemode="tozero")
    fig.update_xaxes(title_text="", showticklabels=False)
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    fig.update_layout(
        showlegend=False,
        margin=dict(l=40, r=20, t=40, b=20),
        height=400,
    )
    st.plotly_chart(fig, width="stretch")

    display = results.copy()
    display["population"] = display["population"].map(POPULATION_LABELS)
    display["significant"] = display["significant"].map({True: "yes", False: "no"})
    display = display.rename(
        columns={
            "n_responder": "n_yes",
            "n_nonresponder": "n_no",
            "median_responder": "median_yes",
            "median_nonresponder": "median_no",
            "p_raw": "p",
            "p_bonferroni": "p_adj",
            "rank_biserial": "r_rb",
        }
    )
    show_cols = [
        "population",
        "n_yes",
        "n_no",
        "median_yes",
        "median_no",
        "U",
        "p",
        "p_adj",
        "r_rb",
        "significant",
    ]
    table = display[show_cols].copy()
    table["p"] = table["p"].map(_fmt_sig)
    table["p_adj"] = table["p_adj"].map(_fmt_sig)
    table["r_rb"] = table["r_rb"].map(lambda v: f"{v:.3f}")
    st.dataframe(
        table,
        width="stretch",
        hide_index=True,
        column_config={
            "median_yes": st.column_config.NumberColumn(format="%.2f"),
            "median_no": st.column_config.NumberColumn(format="%.2f"),
            "U": st.column_config.NumberColumn(format="%.0f"),
            "p": st.column_config.TextColumn("p", help="Four significant figures. CSV has full precision."),
            "p_adj": st.column_config.TextColumn("p_adj", help="Four significant figures. CSV has full precision."),
            "r_rb": st.column_config.TextColumn("r_rb", help="Three decimal places. CSV has full precision."),
        },
    )
    st.download_button(
        "Download statistics (CSV)",
        display[show_cols].to_csv(index=False),
        file_name="response_statistics.csv",
        mime="text/csv",
        help="Full-precision p-values and effect sizes.",
    )

    significant = results.loc[results["significant"], "population"].tolist()
    if significant:
        labels = ", ".join(POPULATION_LABELS[p] for p in significant)
        st.success(
            f"Significant after Bonferroni correction: {labels}."
        )
    else:
        closest = results.loc[results["p_raw"].idxmin()]
        st.info(
            "No population is significant after Bonferroni correction. "
            f"Closest: {POPULATION_LABELS[closest.population]} "
            f"(raw *p* = {_fmt_p(float(closest.p_raw))}, "
            f"adjusted *p* = {_fmt_p(float(closest.p_bonferroni))}, "
            f"*r* = {float(closest.rank_biserial):.3f}). "
            "There is not enough evidence to treat any of these frequencies "
            "as a predictor of miraclib response."
        )


def render_subset(
    cohort: pd.DataFrame,
    by_project: pd.DataFrame,
    by_response: pd.DataFrame,
    by_sex: pd.DataFrame,
) -> None:
    st.subheader("Baseline melanoma PBMC samples (miraclib)")
    st.caption(
        "Filters: condition = melanoma, treatment = miraclib, "
        "sample_type = PBMC, time_from_treatment_start = 0."
    )

    k1, k2, k3 = st.columns(3)
    k1.metric("Samples", f"{len(cohort):,}")
    k2.metric("Subjects", f"{cohort['subject'].nunique():,}")
    k3.metric("Projects with PBMC", int((by_project["n_samples"] > 0).sum()))

    c1, c2, c3 = st.columns(3)
    with c1:
        fig = _project_bar_figure(by_project)
        st.plotly_chart(fig, width="stretch")
        if (by_project["n_samples"] == 0).any():
            missing = ", ".join(by_project.loc[by_project["n_samples"] == 0, "project"])
            st.caption(f"{missing} is whole blood only, so it contributes 0 PBMC samples.")
    with c2:
        resp = by_response.copy()
        resp["label"] = resp["response"].map(lambda v: RESPONSE_LABELS.get(v, v))
        fig = px.bar(
            resp,
            x="label",
            y="n_subjects",
            text="n_subjects",
            title="Subjects by response",
            color="label",
            color_discrete_map={"non-responder": "#E76F51", "responder": "#2A9D8F"},
        )
        fig.update_traces(textposition="outside", cliponaxis=False, showlegend=False)
        fig.update_layout(
            yaxis_title="Subjects",
            xaxis_title="",
            yaxis_range=[0, int(resp["n_subjects"].max()) * 1.2],
            margin=dict(t=50, b=40),
            height=320,
        )
        st.plotly_chart(fig, width="stretch")
    with c3:
        sex = by_sex.copy()
        sex["label"] = sex["sex"].map(lambda v: SEX_LABELS.get(v, v))
        fig = px.bar(
            sex,
            x="label",
            y="n_subjects",
            text="n_subjects",
            title="Subjects by sex",
            color_discrete_sequence=["#457B9D"],
        )
        fig.update_traces(textposition="outside", cliponaxis=False)
        fig.update_layout(
            yaxis_title="Subjects",
            xaxis_title="",
            yaxis_range=[0, int(sex["n_subjects"].max()) * 1.2],
            margin=dict(t=50, b=40),
            height=320,
        )
        st.plotly_chart(fig, width="stretch")

    sample_query = st.text_input("Search identified samples", placeholder="sample, subject, or project")
    view = cohort
    if sample_query.strip():
        q = sample_query.strip()
        mask = (
            view["sample"].str.contains(q, case=False, regex=False)
            | view["subject"].str.contains(q, case=False, regex=False)
            | view["project"].str.contains(q, case=False, regex=False)
        )
        view = view[mask]
    st.dataframe(view, width="stretch", hide_index=True)
    st.caption(f"{len(view):,} of {len(cohort):,} baseline samples shown.")


def main() -> None:
    st.set_page_config(
        page_title="Miraclib immune-cell analysis",
        page_icon="🔬",
        layout="wide",
    )
    _ensure_database()

    st.title("Miraclib immune-cell analysis")
    st.caption(
        "Interactive view of Bob Loblaw’s trial analysis, cell-population "
        "frequencies, responder comparison, and the baseline subset."
    )

    freq = load_frequencies()
    cohort, means, results = load_response_analysis()
    baseline, by_project, by_response, by_sex = load_baseline_subset()

    overview, tab_freq, tab_resp, tab_subset = st.tabs(
        [
            "Overview",
            "Part 2 · Frequencies",
            "Part 3 · Response",
            "Part 4 · Baseline",
        ]
    )
    with overview:
        render_overview(len(freq), cohort, results, baseline)
    with tab_freq:
        render_frequencies(freq)
    with tab_resp:
        render_response(cohort, means, results)
    with tab_subset:
        render_subset(baseline, by_project, by_response, by_sex)


if __name__ == "__main__":
    main()
