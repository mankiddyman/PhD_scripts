#!/usr/bin/env python3
import csv
from pathlib import Path

INPUT_FASTA = Path("/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/paradoxa/Files/genomes/raw/out_JBAT_V2.FINAL.fa")
MAPPING_TSV = Path("/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/paradoxa/Files/genomes/metadata/scaffold_to_chr_assignment.tsv")

OUT_RENAMED = Path("/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/paradoxa/Files/genomes/derived/renamed/D_paradoxa_renamed.fa")
OUT_HAP1 = Path("/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/paradoxa/Files/genomes/derived/haplotypes/D_paradoxa_hap1_chr.fasta")
OUT_HAP2 = Path("/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/paradoxa/Files/genomes/derived/haplotypes/D_paradoxa_hap2_chr.fasta")


def read_fasta(path):
    header = None
    seq_chunks = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            yield header, "".join(seq_chunks)


def wrap_seq(seq, width=60):
    for i in range(0, len(seq), width):
        yield seq[i:i+width]


def load_mapping(path):
    mapping = {}
    seen_new = set()

    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"old_id", "new_id", "haplotype"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise ValueError(
                "Mapping file must have tab-separated columns: old_id, new_id, haplotype"
            )

        for row in reader:
            old_id = row["old_id"].strip()
            new_id = row["new_id"].strip()
            hap = row["haplotype"].strip()

            if hap not in {"hap1", "hap2"}:
                raise ValueError(f"Invalid haplotype '{hap}' for {old_id}. Use hap1 or hap2.")

            if old_id in mapping:
                raise ValueError(f"Duplicate old_id in mapping: {old_id}")

            if new_id in seen_new:
                raise ValueError(f"Duplicate new_id in mapping: {new_id}")

            mapping[old_id] = {"new_id": new_id, "haplotype": hap}
            seen_new.add(new_id)

    return mapping


def main():
    mapping = load_mapping(MAPPING_TSV)

    OUT_RENAMED.parent.mkdir(parents=True, exist_ok=True)
    OUT_HAP1.parent.mkdir(parents=True, exist_ok=True)
    OUT_HAP2.parent.mkdir(parents=True, exist_ok=True)

    found = set()

    with open(OUT_RENAMED, "w") as renamed, open(OUT_HAP1, "w") as hap1, open(OUT_HAP2, "w") as hap2:
        for old_id, seq in read_fasta(INPUT_FASTA):
            if old_id not in mapping:
                continue

            new_id = mapping[old_id]["new_id"]
            hap = mapping[old_id]["haplotype"]
            found.add(old_id)

            renamed.write(f">{new_id}\n")
            for chunk in wrap_seq(seq):
                renamed.write(chunk + "\n")

            out_handle = hap1 if hap == "hap1" else hap2
            out_handle.write(f">{new_id}\n")
            for chunk in wrap_seq(seq):
                out_handle.write(chunk + "\n")

    missing = set(mapping.keys()) - found
    if missing:
        raise ValueError(
            "These scaffolds were listed in the mapping file but not found in the FASTA:\n"
            + "\n".join(sorted(missing))
        )

    print("Done.")
    print(f"Renamed FASTA: {OUT_RENAMED}")
    print(f"Hap1 FASTA:    {OUT_HAP1}")
    print(f"Hap2 FASTA:    {OUT_HAP2}")


if __name__ == "__main__":
    main()
