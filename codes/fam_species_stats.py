
import pandas as pd
import re
from collections import defaultdict
from ete3 import NCBITaxa

ncbi = NCBITaxa()

def load_covered_fams(excel_path):
    df = pd.read_excel(excel_path, usecols=["Rfam ID"])
    rfam_ids = set(df["Rfam ID"].dropna().astype(str).str.strip())
    return rfam_ids

def resolve_species_names(tax_ids):
    numeric_ids = [int(tid.strip("_")) for tid in tax_ids]
    translator = ncbi.get_taxid_translator(numeric_ids)
    return {f"_{tid}": name for tid, name in translator.items()}

def extract_species_stats(seed_path, excel_path, output_file="rfam_species_stats.txt"):
    covered_fams = load_covered_fams(excel_path)
    family_stats = {}
    current_family = None
    skip_family = False
    all_tax_ids = set()

    with open(seed_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()

            if line.startswith("#=GF AC"):
                fam_ac = line.split()[2]
                current_family = fam_ac
                skip_family = fam_ac in covered_fams
                continue

            if skip_family or not current_family:
                continue

            if not line.startswith("#") and line:
                parts = line.split()
                if len(parts) == 2:
                    seq_id, seq = parts
                    tax_match = re.search(r"_(\d+)", seq_id)
                    if tax_match:
                        tax_id = f"_{tax_match.group(1)}"
                        all_tax_ids.add(tax_id)
                        fam_data = family_stats.setdefault(current_family, {
                            "total": 0,
                            "human": 0,
                            "species": defaultdict(int)
                        })
                        fam_data["total"] += 1
                        if tax_id == "_9606":
                            fam_data["human"] += 1
                        fam_data["species"][tax_id] += 1

            elif line.startswith("//"):
                current_family = None
                skip_family = False

    species_names = resolve_species_names(all_tax_ids)

    with open(output_file, "w", encoding="utf-8") as out:
        for fam, stats in family_stats.items():
            out.write(f"Family: {fam}\n")
            out.write(f"  Total sequences: {stats['total']}\n")
            human_label = species_names.get("_9606", "Homo sapiens")
            out.write(f"  Human sequences ({human_label}, _9606): {stats['human']}\n")
            out.write(f"  Other species (sorted by count):\n")
            sorted_species = sorted(
                [(tax_id, count) for tax_id, count in stats["species"].items() if tax_id != "_9606"],
                key=lambda x: x[1],
                reverse=True
            )
            for tax_id, count in sorted_species:
                species_name = species_names.get(tax_id, "Unknown species")
                out.write(f"    {species_name} ({tax_id}): {count}\n")
            out.write("\n")

    print(f" Species breakdown written for {len(family_stats)} families. Saved to {output_file}")

if __name__ == "__main__":
    seed_file = r"C:\Users\gundl\Downloads\Rfam.seed\Rfam.seed"
    excel_file = r"C:\Users\gundl\Downloads\Rfam\Rfam_Final_combined_3d_List.xlsx"
    extract_species_stats(seed_file, excel_file, "rfam_species_stats.txt")

