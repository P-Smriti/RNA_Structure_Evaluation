
# Extract all Rfam IDs and Family Names from Seed File to Excel.

import pandas as pd

def extract_rfam_ids_and_names(seed_file, output_excel="rfam_id_name_list.xlsx"):
    rfam_data = []
    current_id = None
    current_name = None

    with open(seed_file, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()

            if line.startswith("#=GF AC"):
                current_id = line.split()[2]
            elif line.startswith("#=GF ID"):
                current_name = " ".join(line.split()[2:])
            elif line.startswith("//"):
                if current_id and current_name:
                    rfam_data.append({"Rfam ID": current_id, "Family Name": current_name})
                current_id = None
                current_name = None

    df = pd.DataFrame(rfam_data)
    df.to_excel(output_excel, index=False)
    print(f"Saved {len(df)} Rfam families to {output_excel}")

if __name__ == "__main__":
    seed_file = r"C:\Users\gundl\Downloads\Rfam.seed\Rfam.seed"
    extract_rfam_ids_and_names(seed_file)
