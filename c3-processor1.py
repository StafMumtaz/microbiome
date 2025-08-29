import pandas as pd

def union_c1_c2(c1_csv="c1-new-data.csv", c2_csv="c2-new-data.csv", out_csv="c4-combined-dataset.csv"):
    # Load
    df1 = pd.read_csv(c1_csv).copy()
    df2 = pd.read_csv(c2_csv).copy()

    # Normalize header whitespace
    df1.columns = df1.columns.astype(str).str.strip()
    df2.columns = df2.columns.astype(str).str.strip()

    # First column name from c1 defines the output schema/order
    first_col = df1.columns[0]
    feature_cols = list(df1.columns[1:])

    # Make c2's first column name match c1 to allow vertical stacking
    if df2.columns[0] != first_col:
        df2 = df2.rename(columns={df2.columns[0]: first_col})

    # Ensure every c1 feature column exists in c2 (fill with zeros if missing)
    for col in feature_cols:
        if col not in df2.columns:
            df2[col] = 0

    # Coerce numeric for feature columns; non-numeric -> 0
    for col in feature_cols:
        df1[col] = pd.to_numeric(df1[col], errors="coerce").fillna(0)
        df2[col] = pd.to_numeric(df2[col], errors="coerce").fillna(0)

    # --- NEW STEP: normalize c2 feature columns to row-wise proportions (sum to 1) ---
    if feature_cols:
        row_sums = df2[feature_cols].sum(axis=1)
        # Avoid divide-by-zero by replacing 0 with 1; those rows remain all zeros
        safe_den = row_sums.replace(0, 1)
        df2[feature_cols] = df2[feature_cols].div(safe_den, axis=0)

    # Align orders and drop extras from c2
    df1_aligned = df1[[first_col] + feature_cols]
    df2_aligned = df2[[first_col] + feature_cols]

    # Stack c2 below c1
    out = pd.concat([df1_aligned, df2_aligned], axis=0, ignore_index=True)

    out.to_csv(out_csv, index=False)

if __name__ == "__main__":
    union_c1_c2()
