#!/usr/bin/env python3

import argparse
import pandas as pd
from Bio import SeqIO
# import ipdb

def main():
    parser = argparse.ArgumentParser(
        description="Deduplicate sequences with identical strain names AND identical sequences"
    )
    parser.add_argument("--in-sequences", required=True, help="Input FASTA")
    parser.add_argument("--metadata", required=True, help="Input metadata TSV")
    parser.add_argument("--out-sequences", required=True, help="Output FASTA")
    parser.add_argument("--id-field",default="accession",help="Metadata column matching FASTA record IDs (default: accession)")
    parser.add_argument("--strain-field",default="strain",help="Metadata column containing strain names for deduplication (default: strain)")
    parser.add_argument("--threshold",type=float,default=1,help="Percent identity threshold to consider sequences identical (default: 1)")

    args = parser.parse_args()
    # ------------------
    # Read metadata
    # ------------------
    meta = pd.read_csv(args.metadata, sep="\t", dtype=str)

    if args.id_field not in meta.columns:
        raise ValueError(
            f"ID field '{args.id_field}' not found in metadata columns"
        )

    # ------------------
    # Read FASTA
    # ------------------
    records = list(SeqIO.parse(args.in_sequences, "fasta"))

    fasta_df = pd.DataFrame(
        {
            args.id_field: [rec.id for rec in records],
            "sequence": [str(rec.seq) for rec in records],
        }
    )

    # Metadata accession, strain, length
    meta_df = meta[[args.id_field, args.strain_field, "NCBI_length_genome"]].copy()
    meta_df = meta_df.rename(
        columns={args.id_field: "accession",args.strain_field: "strain", "NCBI_length_genome": "length"}
    )
    meta_df = meta_df.loc[meta_df["accession"].isin(fasta_df[args.id_field])]
    meta_df["length"] = pd.to_numeric(meta_df["length"], errors="coerce")

    # Duplicated strains sorted by length
    dup_strains = meta_df[meta_df.duplicated(subset=["strain"])].sort_values("length", ascending=False)

    # Array with accessions for each duplicated strains
    acc_arrays = (
        dup_strains.groupby("strain")["accession"]
        .apply(list)
        .reset_index(name="accession_array")
    )


    pid=0
    sequences_to_remove = set()
    # Percent identical of sequences with same strain name
    for strain,acc in acc_arrays.itertuples(index=False):
        seqs = fasta_df.loc[fasta_df[args.id_field].isin(acc),"sequence"].tolist()
        ref_seq = seqs[0]
        for seq in seqs[1:]:
            matches = sum(1 for a, b in zip(ref_seq, seq) if a == b)
            pid = matches / max(len(ref_seq), len(seq))
            # print(f"Strain: {strain}, PID between {acc[0]} and {acc[1]}: {pid*100:.2f}%")
        if pid > args.threshold:
           sequences_to_remove.update(acc[1:])  # keep first accession (longest)


    # ------------------
    # Deduplicate
    # ------------------
    dedup_fasta_df = fasta_df[
    ~fasta_df[args.id_field].isin(sequences_to_remove)
    ]

    keep_ids = set(dedup_fasta_df[args.id_field])

    with open(args.out_sequences, "w") as out_fasta:
        SeqIO.write(
            (rec for rec in records if rec.id in keep_ids),
            out_fasta,
            "fasta",
        )


    # ------------------
    # Summary
    # ------------------
    print(f"Input sequences: {len(records)}")
    print(f"Deduplicated sequences: {len(dedup_fasta_df)}")
    print(f"Removed duplicates: {len(records) - len(dedup_fasta_df)}")


if __name__ == "__main__":
    main()