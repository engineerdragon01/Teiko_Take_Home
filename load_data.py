import sqlite3
import pandas as pd
from pathlib import Path

DB_PATH = Path(__file__).parent / "teiko.db"
CSV_PATH = Path(__file__).parent / "cell-count.csv"

SCHEMA = """
CREATE TABLE subjects (
    subject   TEXT PRIMARY KEY,
    project   TEXT NOT NULL,
    condition TEXT NOT NULL,
    age       INTEGER,
    sex       TEXT,
    treatment TEXT,
    response  TEXT
);

CREATE TABLE samples (
    sample                    TEXT PRIMARY KEY,
    subject                   TEXT NOT NULL REFERENCES subjects(subject),
    sample_type               TEXT NOT NULL,
    time_from_treatment_start INTEGER NOT NULL
);

CREATE TABLE cell_counts (
    sample     TEXT PRIMARY KEY REFERENCES samples(sample),
    b_cell     INTEGER NOT NULL,
    cd8_t_cell INTEGER NOT NULL,
    cd4_t_cell INTEGER NOT NULL,
    nk_cell    INTEGER NOT NULL,
    monocyte   INTEGER NOT NULL
);
"""


def init_db(conn: sqlite3.Connection) -> None:
    conn.executescript(SCHEMA)
    conn.commit()


def load_csv(conn: sqlite3.Connection, csv_path: Path) -> None:
    df = pd.read_csv(csv_path)

    # Normalize response: empty string → NULL
    df["response"] = df["response"].replace("", None)

    subjects = (
        df[["subject", "project", "condition", "age", "sex", "treatment", "response"]]
        .drop_duplicates(subset=["subject"])
    )
    subjects.to_sql("subjects", conn, if_exists="append", index=False)
    print(f"  subjects: {len(subjects)} rows")

    samples = df[["sample", "subject", "sample_type", "time_from_treatment_start"]]
    samples.to_sql("samples", conn, if_exists="append", index=False)
    print(f"  samples:  {len(samples)} rows")

    cell_counts = df[["sample", "b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]]
    cell_counts.to_sql("cell_counts", conn, if_exists="append", index=False)
    print(f"  cell_counts: {len(cell_counts)} rows")


def main() -> None:
    if DB_PATH.exists():
        DB_PATH.unlink()

    conn = sqlite3.connect(DB_PATH)
    try:
        print("Initializing schema...")
        init_db(conn)
        print(f"Loading {CSV_PATH.name}...")
        load_csv(conn, CSV_PATH)
        print(f"Done. Database saved to {DB_PATH}")
    finally:
        conn.close()


if __name__ == "__main__":
    main()
