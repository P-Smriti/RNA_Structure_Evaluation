#!/usr/bin/env python3
"""
Fill missing chain IDs in FARFAR2 PDBs using Biopython and save to a flat 'str' directory.

Input root:
  /home/s081p868/scratch/RNA_Structure_Evaluation/farfar2/preds

Output dir:
  /home/s081p868/scratch/RNA_Structure_Evaluation/predictions/farfar2/str

Each input PDB (possibly nested in subfolders) is written to the output dir with the
same basename. Existing non-blank chain IDs are preserved; blank chains become A, B, C, ...
"""

import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Tuple, List

from Bio.PDB import PDBParser, PDBIO

IN_ROOT = Path("/home/s081p868/scratch/RNA_Structure_Evaluation/farfar2/preds")
OUT_DIR = Path("/home/s081p868/scratch/RNA_Structure_Evaluation/predictions/farfar2/str")

# Adjust parallelism for your node
MAX_WORKERS = 16

# Chain ID pool to assign for blanks
CHAIN_POOL: List[str] = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz")


def assign_chain_ids(in_pdb: Path, out_pdb: Path) -> Tuple[str, str]:
    """
    Fix blank chain IDs in a single PDB.
    Returns (status, message). status in {"ok", "copied", "skip", "error"}.
    """
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(in_pdb.stem, str(in_pdb))

        # Collect blank and used IDs
        used = set()
        blanks = []
        has_atoms = False
        for model in structure:
            for chain in model:
                # If chain contains at least 1 residue/atom, consider it
                for _res in chain:
                    has_atoms = True
                    break
                cid = chain.id
                if cid is None or cid == " ":
                    blanks.append(chain)
                else:
                    used.add(cid)

        if not has_atoms:
            # No atoms? Just copy.
            shutil.copy2(in_pdb, out_pdb)
            return ("copied", f"{in_pdb.name}: no atoms, copied")

        if not blanks:
            # Nothing to fix; copy to preserve speed and format
            shutil.copy2(in_pdb, out_pdb)
            return ("copied", f"{in_pdb.name}: no blank chains, copied")

        # Assign new IDs to blank chains avoiding collisions
        pool_idx = 0
        for ch in blanks:
            while pool_idx < len(CHAIN_POOL) and CHAIN_POOL[pool_idx] in used:
                pool_idx += 1
            if pool_idx == len(CHAIN_POOL):
                return ("error", f"{in_pdb.name}: ran out of chain IDs to assign")
            new_id = CHAIN_POOL[pool_idx]
            pool_idx += 1
            ch.id = new_id
            used.add(new_id)

        # Save fixed structure
        out_pdb.parent.mkdir(parents=True, exist_ok=True)
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(out_pdb))
        return ("ok", f"{in_pdb.name}: fixed {len(blanks)} blank chain(s)")

    except Exception as e:
        return ("error", f"{in_pdb.name}: {e}")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Find all .pdb files under the input root (recursive)
    pdb_files = sorted(IN_ROOT.rglob("*.pdb"))
    if not pdb_files:
        print(f"No PDBs found under {IN_ROOT}")
        return

    print(f"Found {len(pdb_files)} PDB(s). Writing to {OUT_DIR}")

    futures = []
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
        for in_pdb in pdb_files:
            out_pdb = OUT_DIR / in_pdb.name  # flat output with same basename
            futures.append(ex.submit(assign_chain_ids, in_pdb, out_pdb))

        ok = copied = skipped = errors = 0
        for fut in as_completed(futures):
            status, msg = fut.result()
            print(msg)
            if status == "ok":
                ok += 1
            elif status == "copied":
                copied += 1
            elif status == "skip":
                skipped += 1
            elif status == "error":
                errors += 1

    print(f"\nSummary → fixed: {ok}, copied: {copied}, skipped: {skipped}, errors: {errors}")
    print(f"Output dir: {OUT_DIR}")


if __name__ == "__main__":
    main()
