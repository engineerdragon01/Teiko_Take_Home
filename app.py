"""
app.py — Interactive Streamlit dashboard for Teiko Take-Home.

Usage:  streamlit run app.py
"""

import sqlite3
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from analysis import (
    get_baseline_melanoma_miraclib,
    get_cell_frequencies,
    get_responder_data,
    make_boxplot,
    run_statistical_tests,
)

DB_PATH = Path(__file__).parent / "teiko.db"

st.set_page_config(
    page_title="Teiko: Immune Cell Analysis",
    page_icon="🧬",
    layout="wide",
)


@st.cache_resource
def get_connection():
    return sqlite3.connect(DB_PATH, check_same_thread=False)


@st.cache_data
def load_frequencies():
    return get_cell_frequencies(get_connection())


@st.cache_data
def load_responder_data():
    return get_responder_data(get_connection())


@st.cache_data
def load_stats():
    return run_statistical_tests(load_responder_data())


@st.cache_data
def load_subset():
    return get_baseline_melanoma_miraclib(get_connection())


def check_db():
    if not DB_PATH.exists():
        st.error(
            f"Database not found at `{DB_PATH}`. "
            "Please run `python load_data.py` first (or `make pipeline`)."
        )
        st.stop()


def main():
    check_db()

    st.title("Immune Cell Population Analysis")
    st.caption("Clinical trial data — Bob Loblaw Bio drug candidate")

    tab1, tab2, tab3 = st.tabs([
        "📊 Data Overview",
        "🔬 Statistical Analysis",
        "🔍 Subset Analysis",
    ])

    # ------------------------------------------------------------------
    # Tab 1: Cell Frequency Table
    # ------------------------------------------------------------------
    with tab1:
        st.header("Cell Population Frequencies by Sample")
        st.markdown(
            "Relative frequency of each immune cell population as a percentage "
            "of the total cell count per sample."
        )

        freq_df = load_frequencies()

        # Sidebar-style filters inside the tab
        col_f1, col_f2 = st.columns([1, 3])
        with col_f1:
            pop_filter = st.multiselect(
                "Filter by population",
                options=sorted(freq_df["population"].unique()),
                default=[],
                placeholder="All populations",
            )
        with col_f2:
            sample_search = st.text_input("Search sample ID", placeholder="e.g. sample00042")

        display_df = freq_df.copy()
        if pop_filter:
            display_df = display_df[display_df["population"].isin(pop_filter)]
        if sample_search:
            display_df = display_df[display_df["sample"].str.contains(sample_search, case=False)]

        fmt_df = display_df.copy()
        fmt_df["percentage"] = fmt_df["percentage"].map("{:.2f}%".format)
        fmt_df["count"] = fmt_df["count"].map("{:,}".format)
        fmt_df["total_count"] = fmt_df["total_count"].map("{:,}".format)
        st.dataframe(fmt_df, use_container_width=True, height=440)
        st.caption(f"{len(display_df):,} rows shown")

        st.subheader("Population Average Frequencies (All Samples)")
        avg = (
            freq_df.groupby("population")["percentage"]
            .agg(["mean", "median", "std"])
            .round(3)
            .rename(columns={"mean": "Mean %", "median": "Median %", "std": "Std Dev"})
            .reset_index()
        )
        st.dataframe(avg, use_container_width=True, hide_index=True)

    # ------------------------------------------------------------------
    # Tab 2: Statistical Analysis
    # ------------------------------------------------------------------
    with tab2:
        st.header("Responders vs Non-Responders")
        st.markdown(
            "Melanoma patients treated with **miraclib**, PBMC samples only. "
            "Comparing relative frequencies between responders and non-responders."
        )

        resp_df = load_responder_data()
        stats_df = load_stats()

        # Boxplot
        fig = make_boxplot(resp_df)
        st.plotly_chart(fig, use_container_width=True)

        # Stats table with highlighting
        st.subheader("Mann-Whitney U Test Results")
        st.markdown(
            "Non-parametric test appropriate for clinical data that may not follow "
            "a normal distribution. Significance threshold: **p < 0.05**."
        )

        display_stats = stats_df[[
            "population_label", "n_responders", "n_non_responders",
            "median_responders", "median_non_responders", "p_value", "significant"
        ]].rename(columns={
            "population_label": "Population",
            "n_responders": "N (Responders)",
            "n_non_responders": "N (Non-Responders)",
            "median_responders": "Median % (R)",
            "median_non_responders": "Median % (NR)",
            "p_value": "p-value",
            "significant": "Significant (p<0.05)",
        })

        def highlight_sig(row):
            if row["Significant (p<0.05)"]:
                return ["background-color: #e65100; color: #ffffff"] * len(row)
            return [""] * len(row)

        st.dataframe(
            display_stats.style.apply(highlight_sig, axis=1).format({"p-value": "{:.6f}"}),
            use_container_width=True,
            hide_index=True,
        )

        sig_pops = stats_df.loc[stats_df["significant"], "population_label"].tolist()
        if sig_pops:
            st.success(
                f"**Significant differences found** (p < 0.05): {', '.join(sig_pops)}"
            )
        else:
            st.info("No populations reached statistical significance (p < 0.05).")

    # ------------------------------------------------------------------
    # Tab 3: Subset Analysis
    # ------------------------------------------------------------------
    with tab3:
        st.header("Baseline Melanoma Subset — Miraclib Treated")
        st.markdown(
            "Melanoma PBMC samples at **baseline** (time from treatment start = 0) "
            "from patients treated with **miraclib**."
        )

        subset = load_subset()
        n_samples = len(subset["filtered_samples"])
        n_subjects = subset["filtered_samples"]["subject"].nunique()

        m1, m2 = st.columns(2)
        m1.metric("Total Samples", n_samples)
        m2.metric("Unique Subjects", n_subjects)

        col1, col2, col3 = st.columns(3)

        with col1:
            st.subheader("Samples per Project")
            st.dataframe(subset["samples_per_project"], use_container_width=True, hide_index=True)

        with col2:
            st.subheader("Subjects by Response")
            st.dataframe(subset["response_breakdown"], use_container_width=True, hide_index=True)
            fig_r = go.Figure(go.Pie(
                labels=subset["response_breakdown"]["response"],
                values=subset["response_breakdown"]["subject_count"],
                hole=0.4,
                marker_colors=["#2196F3", "#F44336"],
            ))
            fig_r.update_layout(
                title="Response Distribution",
                height=280,
                margin=dict(t=40, b=0, l=0, r=0),
                showlegend=True,
            )
            st.plotly_chart(fig_r, use_container_width=True)

        with col3:
            st.subheader("Subjects by Sex")
            st.dataframe(subset["sex_breakdown"], use_container_width=True, hide_index=True)
            fig_s = go.Figure(go.Pie(
                labels=subset["sex_breakdown"]["sex"],
                values=subset["sex_breakdown"]["subject_count"],
                hole=0.4,
                marker_colors=["#9C27B0", "#FF9800"],
            ))
            fig_s.update_layout(
                title="Sex Distribution",
                height=280,
                margin=dict(t=40, b=0, l=0, r=0),
                showlegend=True,
            )
            st.plotly_chart(fig_s, use_container_width=True)

        with st.expander("View filtered sample list"):
            st.dataframe(subset["filtered_samples"], use_container_width=True, hide_index=True)


if __name__ == "__main__":
    main()
