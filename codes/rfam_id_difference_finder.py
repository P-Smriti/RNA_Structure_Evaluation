# Compare Rfam Coverage: Identify Families Not Present in PDB List.
import pandas as pd

def find_rfam_id_difference(master_excel, subset_excel, output_excel="rfam_id_difference.xlsx"):
    # Load Rfam IDs from both files
    master_df = pd.read_excel(master_excel, usecols=["Rfam ID"])
    subset_df = pd.read_excel(subset_excel, usecols=["Unique Rfam ID"])

    # Normalize column names
    master_ids = set(master_df["Rfam ID"].dropna().astype(str).str.strip())
    subset_ids = set(subset_df["Unique Rfam ID"].dropna().astype(str).str.strip())

    # Find difference
    missing_ids = sorted(master_ids - subset_ids)

    # Save to Excel
    pd.DataFrame(missing_ids, columns=["Rfam ID not in Unique List"]).to_excel(output_excel, index=False)
    print(f"Saved {len(missing_ids)} Rfam IDs not found in unique list to {output_excel}")

if __name__ == "__main__":
    master_file = r"C:\Users\gundl\Downloads\rfam_id_name_list.xlsx"
    subset_file = r"C:\Users\gundl\Downloads\unique_rfam_ids.xlsx"
    find_rfam_id_difference(master_file, subset_file)
