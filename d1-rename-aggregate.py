import pandas as pd

IN_CSV = "c8-filtered-data-before-renaming.csv"
OUT_CSV = "c8-aggregated-no-renaming.csv"
TOP_K = 31            # keep the 31 most-abundant feature columns
OTHER_NAME = "Other"  # name of the aggregated remainder column

def main():
    df = pd.read_csv(IN_CSV)
    df.columns = df.columns.astype(str).str.strip()

    first_col = df.columns[0]
    feat_cols = list(df.columns[1:])

    # Coerce features to numeric (non-numeric -> 0)
    feats = df[feat_cols].apply(pd.to_numeric, errors="coerce").fillna(0)

    # Compute column-wise sums and sort descending; tie-break by column name
    sums = feats.sum(axis=0)
    order = (
        pd.DataFrame({"col": sums.index, "sum": sums.values})
        .sort_values(["sum", "col"], ascending=[False, True], kind="mergesort")
    )

    top_cols = order["col"].head(TOP_K).tolist()
    other_cols = [c for c in feat_cols if c not in top_cols]

    # Build output: first column + top-K features
    out = pd.concat([df[[first_col]], feats[top_cols]], axis=1)

    # If any remainder, add an "Other" column as row-wise sum of the rest
    if other_cols:
        out[OTHER_NAME] = feats[other_cols].sum(axis=1)

    out.to_csv(OUT_CSV, index=False)

    # Console summary
    print(f"Kept top {len(top_cols)} features. Aggregated {len(other_cols)} into '{OTHER_NAME}'.")
    if other_cols:
        print("Example aggregated (up to 20):", ", ".join(other_cols[:20]))

if __name__ == "__main__":
    main()
