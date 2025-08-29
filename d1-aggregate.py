import pandas as pd

IN_CSV = "c8-filtered-data-before-renaming.csv"
DICT_CSV = "d2-dictionary.csv"
OUT_CSV = "d3-renamed-aggregated.csv"
TOP_K = 31

def before_dot(s: str) -> str:
    """Substring before first period; trims whitespace; robust to NaN."""
    if pd.isna(s):
        return ""
    s = str(s).strip()
    return s.split(".", 1)[0]

def safe_other_name(existing_cols, base="Other"):
    """Pick an 'Other' name that won't collide with existing columns."""
    name = base
    while name in existing_cols:
        name += "_"
    return name

def main():
    # --- Load main data ---
    df = pd.read_csv(IN_CSV)
    df.columns = df.columns.astype(str).str.strip()

    first_col = df.columns[0]
    feat_cols = list(df.columns[1:])

    # Convert features to numeric; non-numeric -> 0
    feats = df[feat_cols].apply(pd.to_numeric, errors="coerce").fillna(0)

    # --- Keep top-K features by column-sum; aggregate the rest into 'Other' ---
    sums = feats.sum(axis=0)
    order = (
        pd.DataFrame({"col": sums.index, "sum": sums.values})
        .sort_values(["sum", "col"], ascending=[False, True], kind="mergesort")
    )
    top_cols = order["col"].head(TOP_K).tolist()
    other_cols = [c for c in feat_cols if c not in top_cols]

    out = pd.concat([df[[first_col]], feats[top_cols]], axis=1)

    if other_cols:
        other_name = safe_other_name(out.columns, base="Other")
        out[other_name] = feats[other_cols].sum(axis=1)
    else:
        other_name = None

    # --- Load dictionary (expects at least two columns) ---
    d2 = pd.read_csv(DICT_CSV, dtype=str)
    d2.columns = d2.columns.astype(str).str.strip()
    if d2.shape[1] < 2:
        raise ValueError("d2-dictionary.csv must have at least two columns.")

    d2_col1 = d2.iloc[:, 0].fillna("").astype(str).str.strip()  # replacement value
    d2_col2 = d2.iloc[:, 1].fillna("").astype(str).str.strip()  # match source

    # ========================
    # PASS 1: match by prefix before the first period
    # ========================
    d2_keys1 = d2_col2.apply(before_dot)
    map1 = (
        pd.DataFrame({"key": d2_keys1, "val": d2_col1})
        .drop_duplicates(subset="key", keep="first")
    )
    mapping1 = dict(zip(map1["key"], map1["val"]))

    orig_vals = out[first_col].fillna("").astype(str)
    keys1 = orig_vals.apply(before_dot)
    mask1 = keys1.isin(mapping1)
    out.loc[mask1, first_col] = keys1[mask1].map(mapping1).values
    renamed1 = int(mask1.sum())

    # ========================
    # PASS 2: match by first 30 characters (rows not renamed in pass 1)
    # ========================
    d2_keys2 = d2_col2.str[:30]
    map2 = (
        pd.DataFrame({"key": d2_keys2, "val": d2_col1})
        .drop_duplicates(subset="key", keep="first")
    )
    mapping2 = dict(zip(map2["key"], map2["val"]))

    curr_vals = out[first_col].fillna("").astype(str)  # after pass 1
    keys2 = curr_vals.str[:30]
    mask2 = (~mask1) & keys2.isin(mapping2)
    out.loc[mask2, first_col] = keys2[mask2].map(mapping2).values
    renamed2 = int(mask2.sum())

    # --- Save ---
    out.to_csv(OUT_CSV, index=False)

    # --- Console summary ---
    print(f"Kept top {len(top_cols)} features; aggregated {len(other_cols)} into '{other_name or 'N/A'}'.")
    print(f"Renamed {renamed1} rows via 'before-dot' match, plus {renamed2} via 30-char match.")
    print(f"Wrote: {OUT_CSV}")

if __name__ == "__main__":
    main()
