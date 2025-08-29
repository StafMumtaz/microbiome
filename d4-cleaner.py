import pandas as pd
import re
import os # Import the 'os' module to work with file paths

# Define input and output file names
IN_CSV = "d3-renamed-aggregated.csv"
OUT_CSV = "f1-cleaned-ready-general.csv"

# Regex to match strings starting with SDP, followed by non-digits, then digits.
# It correctly handles variations like "SDP 1", "SDP-1", "SDP1", and "SDP 1.00".
# Example: "SDP 1.00" -> group(1) will be "1"
# Example: "SDP-11" -> group(1) will be "11"
SDP_PATTERN = re.compile(r"^SDP\D*(\d+)(?:\..*)?$")

def normalize_token(s: str) -> str:
    """
    Cleans and normalizes a string token.
    - Converts to uppercase and strips whitespace.
    - If it matches the SDP_PATTERN, formats it to SDPXX (e.g., SDP01, SDP11).
    - Otherwise, returns the cleaned string.
    """
    if s is None:
        return ""
    # Ensure the input is a string, then strip whitespace and convert to uppercase
    s = str(s).strip().upper()
    
    # Try to match the SDP pattern from the beginning of the string
    m = SDP_PATTERN.match(s)
    
    # If a match is found, reformat the string
    if m:
        # Extract the captured number (group 1) and convert to an integer
        num = int(m.group(1))
        # Format with a leading zero if the number is less than 10
        return f"SDP{num:02d}" if num < 10 else f"SDP{num}"
        
    # If no match, return the cleaned-up original string
    return s

def main():
    """
    Main function to read, clean, and write the CSV data.
    """
    try:
        # Read the input CSV file
        df = pd.read_csv(IN_CSV)
        print(f"Successfully read {IN_CSV}. Shape: {df.shape}")
    except FileNotFoundError:
        print(f"ERROR: The input file '{IN_CSV}' was not found.")
        print("Please make sure the script is in the same directory as the CSV file, or provide the full path.")
        return

    # Clean up column names by stripping any extra whitespace
    df.columns = df.columns.astype(str).str.strip()
    # Get the name of the first column
    first_col = df.columns[0]
    print(f"Processing the first column: '{first_col}'")

    # 1. Create a boolean mask to identify rows to remove.
    rows_to_remove = df[first_col].astype(str).str.strip().str.upper().str.startswith("NEG")
    removed = rows_to_remove.sum()

    # 2. Keep only the rows that DO NOT start with "NEG" by inverting the mask.
    df = df[~rows_to_remove].copy()

    # 3. Apply the normalization function directly to the filtered DataFrame's column.
    df[first_col] = df[first_col].apply(normalize_token)

    # --- FINAL WRITING STEP ---
    try:
        # Get the absolute path for the output file to be sure where it's being saved
        output_filepath = os.path.abspath(OUT_CSV)
        
        # Write the cleaned DataFrame to the new CSV file
        df.to_csv(output_filepath, index=False)
        
        # Print a confirmation message with the full path
        print(f"\n--- Success! ---")
        print(f"Successfully wrote {len(df)} rows and removed {removed} 'NEG' rows.")
        print(f"Cleaned file saved to: {output_filepath}")

    except Exception as e:
        print(f"\nERROR: Failed to write to {OUT_CSV}. Reason: {e}")

if __name__ == "__main__":
    main()
