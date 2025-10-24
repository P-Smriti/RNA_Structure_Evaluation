

import os
import pandas as pd

input_dir = r"C:\Users\gundl\Downloads\data\split_fa"
output_excel = r"C:\Users\gundl\Desktop\fasta_mapping_with_length.xlsx"

data = []

for filename in sorted(os.listdir(input_dir)):
    if filename.endswith(".fasta"):
        filepath = os.path.join(input_dir, filename)
        with open(filepath, "r") as f:
            lines = f.read().splitlines()
            if not lines:
                continue

            header = lines[0].strip()
            if header.startswith(">"):
                header = header[1:]
            clean_name = header.replace("/", "_").replace("|", "_")

            seq = "".join(lines[1:]).replace(" ", "").replace("\n", "")
            seq_length = len(seq)

            data.append({
                #"fasta_file": filename,
                "file_name": clean_name,
                "sequence_length": seq_length,
                "num_motifs_farfar2": "N/A",     # placeholder
                "num_motifs_rhofold": "N/A",     # placeholder
                "num_motifs_alphafold3": "N/A"   # placeholder
            })

df = pd.DataFrame(data)
df.to_excel(output_excel, index=False)

print(f"Excel saved successfully to: {output_excel}")
