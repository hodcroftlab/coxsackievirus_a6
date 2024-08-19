import argparse
import pandas as pd
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--blast', required=True)
    parser.add_argument('--seqs', required=True)
    parser.add_argument('--out_seqs', required=True)
    parser.add_argument('--range', required=True, choices=["vp1","whole_genome"])
    parser.add_argument('--match_length', required=False, type=int)
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Load BLAST results
    blast_results = pd.read_csv(args.blast, 
                                names=["qseqid","sseqid","pident","length","mismatch","gapopen",
                                        "qstart","qend","sstart","send","evalue","bitscore","qcovs"])
    # Load sequences
    sequences = list(SeqIO.parse(args.seqs, "fasta"))
    
    # Filter blast results and remove duplicates
    blast_results=blast_results.sort_values(['length','evalue'],ascending=[False,True]).drop_duplicates(subset='qseqid', keep='first')
    blast_results.loc[:,["qseqid","length"]].to_csv("vp1/results/blast_vp1_length.csv",index=False,header=False) # needed later for adding to metadata
    selected_seqs = []

    # Process sequences based on the specified length range
    if args.range == "vp1":
        r=(600, 900)
        for seq_record in sequences:
            if (seq_record.id in blast_results.qseqid.unique()) and (r[0] <= len(seq_record)):
                # Get the relevant BLAST result
                blast_hit = blast_results[blast_results.qseqid == seq_record.id].iloc[0]
                
                # Extract region from sequence
                reg_seq = seq_record.seq[blast_hit["qstart"] - 1:blast_hit["qend"]]
                
                # Create new sequence record
                new_seq_record = seq_record[:0]  # Copy metadata
                new_seq_record.seq = reg_seq
                selected_seqs.append(new_seq_record)
    
    elif args.range == "whole_genome":
        r=(6400, 8000)
        for seq_record in sequences:
            seq_length = len(seq_record.seq)
            if r[0] <= seq_length <= r[1]:
                selected_seqs.append(seq_record)

    # Write selected sequences
    SeqIO.write(selected_seqs, args.out_seqs,"fasta")
    print("Filtering retained {} out of {} sequences in this file".format(len(selected_seqs), len(sequences)))

if __name__ == "__main__":
    main()
