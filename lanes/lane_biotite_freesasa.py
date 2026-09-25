"""biotite + freesasa lane: secondary structure and solvent accessibility of
the committed receptor structures (structure QC for the docking control set)."""
import glob, json
import biotite.structure as struc
import biotite.structure.io.pdb as pdb
import freesasa

def main():
    rows = []
    for f in sorted(glob.glob("data/raw/structures/*.pdb")):
        name = f.split("/")[-1][:-4]
        pf = pdb.PDBFile.read(f)
        arr = pf.get_structure(model=1)
        ca = arr[arr.atom_name == "CA"]
        sse = struc.annotate_sse(arr)
        n = len(sse) or 1
        res = freesasa.calc(freesasa.Structure(f))
        rows.append({"structure": name, "n_residues_ca": int(len(ca)),
                     "helix_frac": round(float((sse == "a").mean()), 3),
                     "sheet_frac": round(float((sse == "b").mean()), 3),
                     "total_sasa_A2": round(res.totalArea(), 1)})
        print(name, rows[-1]["helix_frac"], rows[-1]["total_sasa_A2"], flush=True)
    import csv
    with open("results/structure_qc_biotite_freesasa.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    json.dump({"n_structures": len(rows)}, open("results/structure_qc.json", "w"))

if __name__ == "__main__":
    main()
