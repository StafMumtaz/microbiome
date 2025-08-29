import pandas as pd

IN_CSV = "c6-cleaned-dataset.csv"
OUT_CSV = "c8-filtered-data.csv"
MAX_THRESH = 0.01  # drop features whose max across samples is < 1%

def main():
    df = pd.read_csv(IN_CSV)
    df.columns = df.columns.astype(str).str.strip()

    first_col = df.columns[0]
    feat_cols = list(df.columns[1:])

    # Coerce features to numeric; non-numeric -> 0
    feats = df[feat_cols].apply(pd.to_numeric, errors="coerce").fillna(0)

    # Keep columns whose max >= threshold
    keep_mask = (feats.max(axis=0) >= MAX_THRESH)
    kept_feat_cols = list(pd.Index(feat_cols)[keep_mask.values])

    # Build pruned dataframe
    pruned = pd.concat([df[[first_col]], feats.loc[:, keep_mask]], axis=1)

    pruned.to_csv(OUT_CSV, index=False)

    # Console summary
    dropped = len(feat_cols) - len(kept_feat_cols)
    print(f"Dropped {dropped} feature columns with max < {MAX_THRESH}.")
    print(f"Kept {len(kept_feat_cols)} feature columns.")
    if dropped:
        dropped_cols = [c for c in feat_cols if c not in kept_feat_cols]
        print("Example dropped columns (up to 20):", ", ".join(dropped_cols[:20]))

if __name__ == "__main__":
    main()
