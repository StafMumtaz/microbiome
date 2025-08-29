import pandas as pd
import re

MATCH_FILE = "a3-genus-matcher.csv"

def load_matcher(path):
    """
    Build:
      - src_names: ordered unique list from column 1 (trimmed, non-empty)
      - rename_map: {col1 -> col2} for rows where col2 has at least one letter
    """
    m = pd.read_csv(path, dtype=str)
    col1 = m.iloc[:, 0].fillna("").astype(str).str.strip()

    if m.shape[1] >= 2:
        col2 = m.iloc[:, 1].fillna("").astype(str).str.strip()
    else:
        col2 = pd.Series([""] * len(col1))

    valid = col1 != ""
    col1 = col1[valid]
    col2 = col2[valid]

    has_alpha = col2.str.contains(r"[A-Za-z]")  # col2 qualifies if it contains any letter

    src_names = pd.unique(col1)  # ordered unique
    rename_map = {c1: c2 for c1, c2, ok in zip(col1, col2, has_alpha) if ok and c2}

    return list(src_names), rename_map

def collapse_duplicate_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    If the header row has duplicates, collapse them into one column by summing
    values down the rows. Non-numeric values are treated as 0 during the sum.
    The first occurrence position is preserved.
    """
    cols = list(df.columns)
    seen = set()
    new_cols = []
    series_list = []

    for i, name in enumerate(cols):
        if name in seen:
            continue
        # find all columns with this name
        idxs = [j for j, n in enumerate(cols) if n == name]
        if len(idxs) == 1:
            series_list.append(df.iloc[:, idxs[0]])
            new_cols.append(name)
        else:
            sub = df.iloc[:, idxs]
            # coerce to numeric; non-numeric -> NaN -> 0
            sub_num = sub.apply(pd.to_numeric, errors="coerce").fillna(0)
            summed = sub_num.sum(axis=1)
            series_list.append(summed)
            new_cols.append(name)
        seen.add(name)

    new_df = pd.concat(series_list, axis=1)
    new_df.columns = new_cols
    return new_df

def add_missing_rename_and_collapse(source_csv: str, out_csv: str,
                                    src_names: list[str], rename_map: dict[str, str]) -> None:
    df = pd.read_csv(source_csv)

    # Normalize header whitespace to avoid accidental mismatches
    df.columns = df.columns.astype(str).str.strip()

    # 1) Ensure every src name exists as a column (zeros if missing)
    for name in src_names:
        if name not in df.columns:
            df[name] = 0

    # 2) Apply renames even if destination already exists;
    #    duplicates will be collapsed in step 3.
    if rename_map:
        df = df.rename(columns=rename_map)

    # 3) Collapse duplicate columns (same header) by summing row-wise
    df = collapse_duplicate_columns(df)

    df.to_csv(out_csv, index=False)

def main():
    src_names, rename_map = load_matcher(MATCH_FILE)

    # a1-new-data.csv -> c1-new-data.csv
    add_missing_rename_and_collapse("a1-new-data.csv", "c1-new-data.csv", src_names, rename_map)

    # a2-old-data.csv -> c2-new-data.csv
    add_missing_rename_and_collapse("a2-old-data.csv", "c2-new-data.csv", src_names, rename_map)

if __name__ == "__main__":
    main()
