# Code to extract unique rfam ids from pdb_full_region.txt file
import pandas as pd

def extract_unique_rfam_ids(file_path, output_excel="unique_rfam_ids.xlsx"):
    rfam_ids = set()

    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if line:
                rfam_id = line.split("\t")[0]
                rfam_ids.add(rfam_id)

    df = pd.DataFrame(sorted(rfam_ids), columns=["Unique Rfam ID"])
    df.to_excel(output_excel, index=False)
    print(f"Saved {len(df)} unique Rfam IDs to {output_excel}")

if __name__ == "__main__":
    input_file = r"C:\Users\gundl\Downloads\pdb_full_region.txt"
    extract_unique_rfam_ids(input_file)