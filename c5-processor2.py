import pandas as pd

IN_CSV = "c4-combined-dataset.csv"
OUT_CSV = "c5-cleaned-deduped.csv"

def main():
    df = pd.read_csv(IN_CSV)

    # Normalize headers just in case
    df.columns = df.columns.astype(str).str.strip()
    first_col = df.columns[0]

    # 1) Drop columns (2..n) that are all zeros across rows
    zero_cols = []
    for col in df.columns[1:]:
        s = pd.to_numeric(df[col], errors="coerce").fillna(0)
        if (s == 0).all():
            zero_cols.append(col)
    if zero_cols:
        df = df.drop(columns=zero_cols)

    # 2) Remove duplicate rows by first 30 chars of the first column (keep first)
    seen = set()
    drop_idx = []
    dropped_values = []

    # Treat NaNs as empty strings for the key
    first_vals = df[first_col].fillna("").astype(str)

    for idx, val in first_vals.items():
        key = val[:30]  # works fine even if len(val) < 30
        if key in seen:
            drop_idx.append(idx)
            dropped_values.append(val)
        else:
            seen.add(key)

    if drop_idx:
        df = df.drop(index=drop_idx)

    # Reindex after drops
    df = df.reset_index(drop=True)

    # Save
    df.to_csv(OUT_CSV, index=False)

    # Terminal summary
    print(f"Dropped {len(zero_cols)} all-zero columns.")
    if zero_cols:
        print("Columns dropped:", ", ".join(zero_cols))
    print(f"Removed {len(drop_idx)} duplicate rows (by first 30 chars of '{first_col}').")
    if dropped_values:
        print("First-column values removed:")
        for v in dropped_values:
            print(f" - {v}")

if __name__ == "__main__":
    main()
