import pandas as pd

# Inputs → Outputs
PAIRS = [
    ("a1-main-data.csv",          "a3-middle-ear.csv"),
    ("a1b-main-no-aggregate.csv", "a3-me-no-agg.csv"),
]

# Leading whitespace (ASCII + common Unicode space characters)
LEADING_WS = r'^[\u00A0\u202F\u2007\u2009\u2000-\u200B\u2060\s]*'

def _normalize(series: pd.Series) -> pd.Series:
    """Trim weird leading spaces and upper-case for robust matching/sorting."""
    return (
        series.astype(str)
              .str.replace(LEADING_WS, "", regex=True)
              .str.upper()
    )

def filter_and_sort(in_csv: str, out_csv: str, prefix: str = "ME") -> None:
    df = pd.read_csv(in_csv)
    first_col = df.columns[0]

    # Filter: first column starts with 'ME' (after normalization)
    norm = _normalize(df[first_col])
    keep_mask = norm.str.startswith(prefix)
    kept = df.loc[keep_mask].copy()
    removed_n = len(df) - len(kept)

    # Stable sort by a normalized key, but keep original values
    kept = kept.sort_values(
        by=first_col,
        key=_normalize,
        ascending=True,
        kind="mergesort"  # stable
    )

    kept.to_csv(out_csv, index=False)
    print(f"Wrote {out_csv}. Removed {removed_n} rows not starting with '{prefix}' in {first_col} from {in_csv}.")

def main():
    for in_csv, out_csv in PAIRS:
        filter_and_sort(in_csv, out_csv, prefix="ME")

if __name__ == "__main__":
    main()
