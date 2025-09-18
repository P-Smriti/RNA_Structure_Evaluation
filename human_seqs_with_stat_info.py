
import pandas as pd
import re

def load_covered_fams(excel_path):
    df = pd.read_excel(excel_path, usecols=["Rfam ID"])
    rfam_ids = set(df["Rfam ID"].dropna().astype(str).str.strip())
    return rfam_ids

def gc_content(seq):
    ungapped = seq.replace("-", "")
    gc = ungapped.count("G") + ungapped.count("C")
    return round(100 * gc / len(ungapped), 2) if ungapped else 0

def extract_human_sequences(file_path, excel_path, output_file="human_remaining_final.txt"):
    covered_fams = load_covered_fams(excel_path)
    human_seqs = {}
    current_family = None
    skip_family = False

    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()

            # Detect family accession
            if line.startswith("#=GF AC"):
                fam_ac = line.split()[2]
                current_family = fam_ac
                skip_family = fam_ac in covered_fams
                continue

            if skip_family or not current_family:
                continue

            # Sequence lines
            if not line.startswith("#") and line:
                parts = line.split()
                if len(parts) == 2:
                    seq_id, seq = parts
                    if "_9606" in seq_id:
                        human_seqs.setdefault(current_family, []).append((seq_id, seq))

            # End of family
            elif line.startswith("//"):
                current_family = None
                skip_family = False

    # Save results with stats
    with open(output_file, "w", encoding="utf-8") as out:
        for fam, seqs in human_seqs.items():
            raw_lengths = [len(seq) for _, seq in seqs]
            ungapped_lengths = [len(seq.replace("-", "")) for _, seq in seqs]
            gc_values = [gc_content(seq) for _, seq in seqs]

            out.write(f"Family: {fam}\n")
            out.write(f"Human sequences: {len(seqs)}\n")
            out.write(f"  Avg length (with gaps): {sum(raw_lengths) // len(raw_lengths)} bp\n")
            out.write(f"  Avg length (no gaps): {sum(ungapped_lengths) // len(ungapped_lengths)} bp\n")
            out.write(f"  Min ungapped length: {min(ungapped_lengths)} bp\n")
            out.write(f"  Max ungapped length: {max(ungapped_lengths)} bp\n")
            out.write(f"  Avg GC content: {round(sum(gc_values) / len(gc_values), 2)}%\n\n")

            for seq_id, seq in seqs:
                raw_len = len(seq)
                ungapped_len = len(seq.replace("-", ""))
                gc = gc_content(seq)
                out.write(f"  {seq_id}\n")
                out.write(f"    Length (with gaps): {raw_len} bp\n")
                out.write(f"    Length (no gaps): {ungapped_len} bp\n")
                out.write(f"    GC content: {gc}%\n")
                out.write(f"    Sequence: {seq}\n\n")

            out.write("-" * 50 + "\n")

    print(f"✅ Extracted {sum(len(seqs) for seqs in human_seqs.values())} human sequences "
          f"from {len(human_seqs)} families (excluding those listed in Excel). Saved to {output_file}")

if __name__ == "__main__":
    seed_file = r"C:\Users\gundl\Downloads\Rfam.seed\Rfam.seed"
    excel_file = r"C:\Users\gundl\Downloads\Rfam\Rfam_Final_combined_3d_List.xlsx"
    extract_human_sequences(seed_file, excel_file, "human_remaining_with_stat_info.txt")
