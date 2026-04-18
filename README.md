# Teiko Take-Home: Immune Cell Population Analysis

Interactive analysis of a clinical trial evaluating how **miraclib** affects immune cell populations in melanoma patients.

---

## Quick Start

```bash
make setup      # Install dependencies
make pipeline   # Load data and generate all outputs
make dashboard  # Launch the Streamlit dashboard
```

The dashboard will be available at `http://localhost:8501`.

> **Note for GitHub Codespaces:** After running `make dashboard`, Codespaces will forward port 8501 automatically. Click the "Open in Browser" popup or find it in the **Ports** tab.

---

## Prerequisites

- Python 3.9+
- `cell-count.csv` in the repository root (already included)

---

## Database Schema

### Design

Three normalized tables:

```
subjects
├── subject     TEXT  PRIMARY KEY
├── project     TEXT
├── condition   TEXT  (melanoma, carcinoma, healthy)
├── age         INTEGER
├── sex         TEXT  (M, F)
├── treatment   TEXT  (miraclib, phauximab, none)
└── response    TEXT  (yes, no, NULL for untreated)

samples
├── sample                    TEXT  PRIMARY KEY
├── subject                   TEXT  → subjects.subject
├── sample_type               TEXT  (PBMC, WB)
└── time_from_treatment_start INTEGER

cell_counts
├── sample     TEXT  PRIMARY KEY → samples.sample
├── b_cell     INTEGER
├── cd8_t_cell INTEGER
├── cd4_t_cell INTEGER
├── nk_cell    INTEGER
└── monocyte   INTEGER
```

### Rationale

**Why three tables?**

Subject-level attributes (demographics, treatment assignment, response) are stored once in `subjects` rather than repeated across every sample row. This eliminates update anomalies — if a response annotation changes, it is updated in one place. The `samples` table captures the repeated-measures structure of the study (each subject has multiple samples at different timepoints and sample types). `cell_counts` is separated from `samples` so the schema cleanly distinguishes *what was collected* from *what was measured*, making it easy to add new assay types later.

**Scaling to hundreds of projects, thousands of samples, and varied analytics**

| Concern | Solution |
|---------|----------|
| Many projects | `project` is a column on `subjects`; add a `projects` lookup table with metadata (PI, indication, protocol) and join as needed |
| Thousands of subjects/samples | Index `subjects(project, condition, treatment)` and `samples(subject, sample_type, time_from_treatment_start)` to keep filtering fast |
| New cell populations or assays | Add a `cell_types` lookup table + a long-format `measurements(sample, cell_type_id, count)` table; the current wide format is convenient for 5 populations but doesn't scale to dozens |
| Varied analytics (survival, multi-omics) | Add dedicated tables (e.g., `clinical_outcomes`, `genomics_results`) that FK to `subjects`; the relational model makes cross-domain joins straightforward |
| Multi-user read load | Migrate to PostgreSQL — the SQLite schema is fully compatible with minor dialect changes |

---

## Code Structure

```
.
├── load_data.py         # Part 1 — Schema init + CSV loading
├── analysis.py          # Parts 2-4 — Analysis functions + pipeline runner
├── app.py               # Streamlit dashboard
├── requirements.txt     # Python dependencies
├── Makefile             # setup / pipeline / dashboard targets
├── cell-count.csv       # Input data
├── teiko.db             # Generated SQLite database (after make pipeline)
└── outputs/             # Generated outputs (after make pipeline)
    ├── cell_frequencies.csv
    ├── responder_stats.csv
    ├── responder_boxplot.html
    ├── responder_boxplot.png
    ├── subset_samples_per_project.csv
    ├── subset_response_breakdown.csv
    ├── subset_sex_breakdown.csv
    └── subset_filtered_samples.csv
```

### Design Decisions

**`load_data.py` is self-contained.** The spec requires `python load_data.py` to work without arguments. All schema definitions live in this file so it can be run independently of the rest of the codebase.

**`analysis.py` is dual-purpose.** It exports importable functions (`get_cell_frequencies`, `get_responder_data`, etc.) used by the dashboard, and has a `main()` entry point for the pipeline. This avoids duplicating logic between the batch pipeline and the interactive UI.

**Streamlit for the dashboard.** Streamlit requires minimal boilerplate compared to Dash or Flask, ships with a good table renderer (`st.dataframe`), and integrates natively with Plotly. The `@st.cache_data` decorators ensure database queries only run once per session.

**Mann-Whitney U test for statistics.** Clinical immunology data is rarely normally distributed. Mann-Whitney U is a non-parametric test that makes no distributional assumptions, is robust to outliers, and is standard practice for comparing cell population frequencies between patient groups.

**SQLite for the database.** Appropriate for a single-node analytical workload and zero-configuration setup. The normalized schema is fully portable to PostgreSQL for production use.

---

## Analysis Results Summary

### Part 3 — Responder vs Non-Responder (Melanoma + Miraclib + PBMC)

| Population | Median % (Responders) | Median % (Non-Responders) | p-value | Significant |
|---|---|---|---|---|
| B Cell | 9.43 | 9.79 | 0.056 | No |
| CD8 T Cell | 24.73 | 24.60 | 0.639 | No |
| **CD4 T Cell** | **30.22** | **29.66** | **0.013** | **Yes** |
| NK Cell | 14.51 | 14.80 | 0.121 | No |
| Monocyte | 19.61 | 19.94 | 0.163 | No |

**CD4 T Cell** frequency is significantly higher in responders (p = 0.013), suggesting it may contribute to predicting miraclib response.

### Part 4 — Baseline Melanoma PBMC Subset

Filtered to melanoma patients, miraclib treatment, PBMC samples, at baseline (day 0):

- **656 total samples** from **656 subjects** across 2 projects
- **prj1:** 384 samples · **prj3:** 272 samples
- **Responders:** 331 · **Non-Responders:** 325
- **Male:** 344 · **Female:** 312

---

## Dashboard

**Live dashboard:** *(deploy to Streamlit Cloud and add URL here)*

The dashboard has three tabs:
- **Data Overview** — Searchable/filterable frequency table for all samples
- **Statistical Analysis** — Boxplot comparing responders vs non-responders, Mann-Whitney U results table
- **Subset Analysis** — Summary tables and charts for the baseline miraclib/melanoma/PBMC subset
