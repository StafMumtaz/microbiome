import pandas as pd

MATCH_FILE = "a3-genus-matcher.csv"

def add_missing_columns(source_csv: str, match_csv: str, out_csv: str) -> None:
    # Load the base data
    df = pd.read_csv(source_csv)

    # Load the matcher list: take the FIRST column, skip header automatically,
    # trim whitespace, drop blanks/NaNs, de-duplicate
    matcher = pd.read_csv(match_csv)
    names = (
        matcher.iloc[:, 0]            # first column
        .dropna()
        .astype(str)
        .str.strip()
    )
    names = pd.unique(names)          # de-duplicate, preserves order seen

    # For any matcher name not already a column, add it as a zero-filled column
    for name in names:
        if name and name not in df.columns:
            df[name] = 0

    # Write out the augmented CSV
    df.to_csv(out_csv, index=False)

def main():
    # a1-new-data.csv → b1-new-data.csv
    add_missing_columns("a1-new-data.csv", MATCH_FILE, "b1-new-data.csv")

    # a2-old-data.csv → b2-old-data.csv
    add_missing_columns("a2-old-data.csv", MATCH_FILE, "b2-old-data.csv")

if __name__ == "__main__":
    main()
