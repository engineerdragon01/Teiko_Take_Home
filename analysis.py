"""
analysis.py — Parts 2-4 pipeline.

Run directly:  python analysis.py
Outputs written to ./outputs/
"""

import sqlite3
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy import stats

DB_PATH = Path(__file__).parent / "teiko.db"
OUT_DIR = Path(__file__).parent / "outputs"

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
POPULATION_LABELS = {
    "b_cell": "B Cell",
    "cd8_t_cell": "CD8 T Cell",
    "cd4_t_cell": "CD4 T Cell",
    "nk_cell": "NK Cell",
    "monocyte": "Monocyte",
}


# ---------------------------------------------------------------------------
# Part 2: Cell frequency table
# ---------------------------------------------------------------------------

def get_cell_frequencies(conn: sqlite3.Connection) -> pd.DataFrame:
    """Return per-sample relative frequency for each cell population."""
    query = """
        SELECT
            s.sample,
            cc.b_cell, cc.cd8_t_cell, cc.cd4_t_cell, cc.nk_cell, cc.monocyte
        FROM samples s
        JOIN cell_counts cc ON s.sample = cc.sample
    """
    df = pd.read_sql_query(query, conn)

    df["total_count"] = df[POPULATIONS].sum(axis=1)

    long = df.melt(
        id_vars=["sample", "total_count"],
        value_vars=POPULATIONS,
        var_name="population",
        value_name="count",
    )
    long["percentage"] = (long["count"] / long["total_count"] * 100).round(4)

    return long[["sample", "total_count", "population", "count", "percentage"]].reset_index(drop=True)


# ---------------------------------------------------------------------------
# Part 3: Responder vs non-responder comparison
# ---------------------------------------------------------------------------

def get_responder_data(conn: sqlite3.Connection) -> pd.DataFrame:
    """Fetch melanoma + miraclib + PBMC frequency data for responder comparison."""
    query = """
        SELECT
            s.sample,
            sub.response,
            cc.b_cell, cc.cd8_t_cell, cc.cd4_t_cell, cc.nk_cell, cc.monocyte
        FROM samples s
        JOIN subjects sub ON s.subject = sub.subject
        JOIN cell_counts cc ON s.sample = cc.sample
        WHERE sub.condition = 'melanoma'
          AND sub.treatment = 'miraclib'
          AND s.sample_type = 'PBMC'
          AND sub.response IN ('yes', 'no')
    """
    df = pd.read_sql_query(query, conn)
    df["total_count"] = df[POPULATIONS].sum(axis=1)
    for pop in POPULATIONS:
        df[f"{pop}_pct"] = df[pop] / df["total_count"] * 100
    return df


def run_statistical_tests(responder_df: pd.DataFrame) -> pd.DataFrame:
    """Mann-Whitney U test for each population, responders vs non-responders."""
    rows = []
    for pop in POPULATIONS:
        col = f"{pop}_pct"
        yes = responder_df.loc[responder_df["response"] == "yes", col].dropna()
        no = responder_df.loc[responder_df["response"] == "no", col].dropna()
        stat, p = stats.mannwhitneyu(yes, no, alternative="two-sided")
        rows.append({
            "population": pop,
            "population_label": POPULATION_LABELS[pop],
            "n_responders": len(yes),
            "n_non_responders": len(no),
            "median_responders": round(yes.median(), 3),
            "median_non_responders": round(no.median(), 3),
            "mann_whitney_u": round(stat, 2),
            "p_value": round(p, 6),
            "significant": p < 0.05,
        })
    return pd.DataFrame(rows)


def make_boxplot(responder_df: pd.DataFrame) -> go.Figure:
    """Faceted boxplot: one subplot per population, responders vs non-responders."""
    fig = make_subplots(
        rows=1, cols=5,
        subplot_titles=[POPULATION_LABELS[p] for p in POPULATIONS],
        shared_yaxes=False,
    )

    colors = {"yes": "#2196F3", "no": "#F44336"}
    labels = {"yes": "Responders", "no": "Non-Responders"}
    shown = set()

    for col_idx, pop in enumerate(POPULATIONS, start=1):
        pct_col = f"{pop}_pct"
        for resp in ["yes", "no"]:
            subset = responder_df.loc[responder_df["response"] == resp, pct_col]
            show_legend = resp not in shown
            shown.add(resp)
            fig.add_trace(
                go.Box(
                    y=subset,
                    name=labels[resp],
                    marker_color=colors[resp],
                    legendgroup=resp,
                    showlegend=show_legend,
                    boxpoints="all",
                    jitter=0.3,
                    pointpos=0,
                    marker=dict(size=4, opacity=0.6),
                ),
                row=1, col=col_idx,
            )

    fig.update_layout(
        title="Cell Population Frequencies: Responders vs Non-Responders<br>"
              "<sup>Melanoma patients treated with miraclib (PBMC samples)</sup>",
        yaxis_title="Relative Frequency (%)",
        legend_title="Response",
        height=520,
        template="plotly_white",
        font=dict(size=12),
    )
    for i in range(1, 6):
        fig.update_yaxes(title_text="%" if i == 1 else "", row=1, col=i)

    return fig


# ---------------------------------------------------------------------------
# Part 4: Subset analysis
# ---------------------------------------------------------------------------

def get_baseline_melanoma_miraclib(conn: sqlite3.Connection) -> dict:
    """
    Baseline melanoma PBMC samples from miraclib-treated patients.
    Returns dict with three summary DataFrames.
    """
    base_query = """
        SELECT
            s.sample,
            sub.subject,
            sub.project,
            sub.response,
            sub.sex
        FROM samples s
        JOIN subjects sub ON s.subject = sub.subject
        WHERE sub.condition    = 'melanoma'
          AND s.sample_type    = 'PBMC'
          AND s.time_from_treatment_start = 0
          AND sub.treatment    = 'miraclib'
    """
    df = pd.read_sql_query(base_query, conn)

    samples_per_project = (
        df.groupby("project")["sample"]
        .count()
        .reset_index()
        .rename(columns={"sample": "sample_count"})
    )

    response_breakdown = (
        df.drop_duplicates("subject")
        .groupby("response")["subject"]
        .count()
        .reset_index()
        .rename(columns={"subject": "subject_count"})
    )

    sex_breakdown = (
        df.drop_duplicates("subject")
        .groupby("sex")["subject"]
        .count()
        .reset_index()
        .rename(columns={"subject": "subject_count"})
    )

    return {
        "filtered_samples": df,
        "samples_per_project": samples_per_project,
        "response_breakdown": response_breakdown,
        "sex_breakdown": sex_breakdown,
    }


# ---------------------------------------------------------------------------
# Pipeline entry point
# ---------------------------------------------------------------------------

def main() -> None:
    if not DB_PATH.exists():
        raise FileNotFoundError(
            f"Database not found at {DB_PATH}. Run `python load_data.py` first."
        )

    OUT_DIR.mkdir(exist_ok=True)
    conn = sqlite3.connect(DB_PATH)

    try:
        # Part 2
        print("=" * 60)
        print("Part 2: Cell Frequency Table")
        print("=" * 60)
        freq_df = get_cell_frequencies(conn)
        freq_df.to_csv(OUT_DIR / "cell_frequencies.csv", index=False)
        print(freq_df.head(10).to_string(index=False))
        print(f"\n→ Saved {len(freq_df)} rows to outputs/cell_frequencies.csv")

        # Part 3
        print("\n" + "=" * 60)
        print("Part 3: Responder vs Non-Responder Statistical Analysis")
        print("=" * 60)
        resp_df = get_responder_data(conn)
        stats_df = run_statistical_tests(resp_df)
        stats_df.to_csv(OUT_DIR / "responder_stats.csv", index=False)
        print(stats_df[["population_label", "n_responders", "n_non_responders",
                         "median_responders", "median_non_responders",
                         "p_value", "significant"]].to_string(index=False))

        sig = stats_df.loc[stats_df["significant"], "population_label"].tolist()
        print(f"\nSignificant populations (p < 0.05): {sig if sig else 'None'}")
        print("→ Saved to outputs/responder_stats.csv")

        fig = make_boxplot(resp_df)
        fig.write_html(OUT_DIR / "responder_boxplot.html")
        print("→ Saved outputs/responder_boxplot.html")
        try:
            fig.write_image(OUT_DIR / "responder_boxplot.png", scale=2)
            print("→ Saved outputs/responder_boxplot.png")
        except Exception:
            pass  # kaleido may not be available in all environments

        # Part 4
        print("\n" + "=" * 60)
        print("Part 4: Baseline Melanoma + Miraclib Subset Analysis")
        print("=" * 60)
        subset = get_baseline_melanoma_miraclib(conn)

        print(f"\nTotal filtered samples: {len(subset['filtered_samples'])}")

        print("\nSamples per project:")
        print(subset["samples_per_project"].to_string(index=False))

        print("\nSubjects by response:")
        print(subset["response_breakdown"].to_string(index=False))

        print("\nSubjects by sex:")
        print(subset["sex_breakdown"].to_string(index=False))

        subset["samples_per_project"].to_csv(OUT_DIR / "subset_samples_per_project.csv", index=False)
        subset["response_breakdown"].to_csv(OUT_DIR / "subset_response_breakdown.csv", index=False)
        subset["sex_breakdown"].to_csv(OUT_DIR / "subset_sex_breakdown.csv", index=False)
        subset["filtered_samples"].to_csv(OUT_DIR / "subset_filtered_samples.csv", index=False)
        print("\n→ Saved subset CSVs to outputs/")

    finally:
        conn.close()


if __name__ == "__main__":
    main()
