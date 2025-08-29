import pandas as pd
import re
import os

# I/O
IN_CSV   = "c8-filtered-data-before-renaming.csv"
DICT_CSV = "d2-dictionary.csv"
OUT_CSV  = "f2-non-aggregate.csv"

# --- Helpers ---

def before_dot(s: str) -> str:
    """Substring before first period; trims whitespace; robust to NaN."""
    if pd.isna(s):
        return ""
    s = str(s).strip()
    return s.split(".", 1)[0]

# Matches SDP tokens like "SDP 1", "SDP-1", "SDP1", "SDP 1.00"
_SDP_RE = re.compile(r"^SDP\D*(\d+)(?:\..*)?$")

def normalize_sdp_token(s: str) -> str:
    """
    Normalize SDP-style IDs to SDPXX / SDP##.
    Only normalizes if the token looks like SDP; otherwise returns as-is.
    """
    if s is None:
        return ""
    raw = str(s).strip()
    up  = raw.upper()
    m = _SDP_RE.match(up)
    if m:
        num = int(m.group(1))
        return f"SDP{num:02d}" if num < 10 else f"SDP{num}"
    return raw  # preserve original (including case) if not SDP-like

def main():
    # --- Load data ---
    try:
        df = pd.read_csv(IN_CSV)
        print(f"Loaded '{IN_CSV}' with shape {df.shape}.")
    except FileNotFoundError:
        print(f"ERROR: '{IN_CSV}' not found.")
        return

    # Clean column names
    df.columns = df.columns.astype(str).str.strip()
    first_col = df.columns[0]
    feat_cols = list(df.columns[1:])
    print(f"First column: '{first_col}'. Feature columns: {len(feat_cols)}")

    # --- Drop NEG rows (case-insensitive, leading/trailing spaces ignored) ---
    mask_neg = (
        df[first_col].astype(str).str.strip().str.upper().str.startswith("NEG")
    )
    removed_neg = int(mask_neg.sum())
    if removed_neg:
        print(f"Removing {removed_neg} rows starting with 'NEG' in '{first_col}'.")
    df = df[~mask_neg].copy()

    # --- Two-pass dictionary renaming on the first column ---
    renamed1 = renamed2 = 0
    renamed_any_mask = pd.Series(False, index=df.index)

    try:
        d2 = pd.read_csv(DICT_CSV, dtype=str)
        d2.columns = d2.columns.astype(str).str.strip()
        if d2.shape[1] < 2:
            raise ValueError("d2-dictionary.csv must have at least two columns.")
        # Replacement (col1) and match source (col2)
        d2_col1 = d2.iloc[:, 0].fillna("").astype(str).str.strip()
        d2_col2 = d2.iloc[:, 1].fillna("").astype(str).str.strip()

        # PASS 1: match by prefix before the first period
        d2_keys1 = d2_col2.apply(before_dot)
        map1 = (
            pd.DataFrame({"key": d2_keys1, "val": d2_col1})
            .drop_duplicates(subset="key", keep="first")
        )
        mapping1 = dict(zip(map1["key"], map1["val"]))

        curr_vals = df[first_col].fillna("").astype(str)
        keys1 = curr_vals.apply(before_dot)
        mask1 = keys1.isin(mapping1)
        if mask1.any():
            df.loc[mask1, first_col] = keys1[mask1].map(mapping1).values
            renamed1 = int(mask1.sum())
            renamed_any_mask |= mask1

        # PASS 2: match by first 30 characters (only rows not renamed in pass 1)
        d2_keys2 = d2_col2.str[:30]
        map2 = (
            pd.DataFrame({"key": d2_keys2, "val": d2_col1})
            .drop_duplicates(subset="key", keep="first")
        )
        mapping2 = dict(zip(map2["key"], map2["val"]))

        curr_vals = df[first_col].fillna("").astype(str)  # after pass 1
        keys2 = curr_vals.str[:30]
        mask2 = (~renamed_any_mask) & keys2.isin(mapping2)
        if mask2.any():
            df.loc[mask2, first_col] = keys2[mask2].map(mapping2).values
            renamed2 = int(mask2.sum())
            renamed_any_mask |= mask2

        print(f"Dictionary renaming: {renamed1} via 'before-dot', {renamed2} via 30-char.")

    except FileNotFoundError:
        print(f"WARNING: '{DICT_CSV}' not found. Skipping dictionary renaming.")
    except Exception as e:
        print(f"WARNING: Dictionary renaming skipped due to error: {e}")

    # --- Normalize SDP tokens ONLY for rows not renamed by dictionary ---
    if (~renamed_any_mask).any():
        df.loc[~renamed_any_mask, first_col] = df.loc[
            ~renamed_any_mask, first_col
        ].apply(normalize_sdp_token)

    # --- Coerce all feature columns to numeric; non-numeric -> 0 ---
    if feat_cols:
        feats_numeric = df[feat_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
        df.loc[:, feat_cols] = feats_numeric

    # --- Write output (no aggregation; keep all 2..n columns) ---
    out_path = os.path.abspath(OUT_CSV)
    try:
        df.to_csv(out_path, index=False)
        print("\n--- Success ---")
        print(f"Wrote {len(df)} rows × {len(df.columns)} cols to:\n{out_path}")
        if removed_neg:
            print(f"Removed {removed_neg} 'NEG' rows.")
        print(f"Renamed total via dictionary: {int(renamed_any_mask.sum())}")
        print("Kept ALL feature columns (no aggregation).")
    except Exception as e:
        print(f"ERROR: Failed to write '{OUT_CSV}': {e}")

if __name__ == "__main__":
    main()
