###############
# Snakemake execution templates:

# To run a default VP1 run(<600bp):
# snakemake  auspice/coxsackievirus_A6_vp1.json --cores 9

# To run a default whole genome run ( <6400bp):
# snakemake auspice/coxsackievirus_A6_whole-genome.json --cores 9

import os
from datetime import date

# Load config file
if not config:
    configfile: "config/config.yaml"

# Load environment variables
# Try to load .env, but don't fail if it doesn't exist (for Actions)
try:
    from dotenv import load_dotenv
    load_dotenv(".env")
except:
    pass

TAXID = "86107"
REMOTE_GROUP = os.getenv("REMOTE_GROUP")
UPLOAD_DATE = date.today().isoformat()

DOWNLOAD_INGEST=True

###############
wildcard_constraints:
    seg="vp1|whole_genome|P1",
    gene="|-5utr|-vp4|-vp2|-vp3|-vp1|-2A|-2B|-2C|-3A|-3B|-3C|-3D|-3utr"
   
#     #from: https://bitbucket.org/snakemake/snakemake/issues/910/empty-wildcard-assignment-works-only-if

# Define segments to analyze
segments = ['vp1', 'whole-genome', 'P1'] # add more segments if you want to analyze them separately
GENES=["-5utr","-vp4", "-vp2", "-vp3", "-vp1", "-2A", "-2B", "-2C", "-3A", "-3B", "-3C", "-3D","-3utr"]
CODING_GENES = ["VP4", "VP2", "VP3", "VP1", "2A", "2B", "2C", "3A", "3B", "3C", "3D"]


# Rule to handle configuration files
rule files:
    input:
        sequence_length =   "{seg}",
        dropped_strains =   "config/exclude.txt",
        incl_strains =      "config/include.txt",
        reference =         "config/reference_sequence.gb",
        gff_reference =     "{seg}/config/annotation.gff3",
        lat_longs =         "config/lat_longs.tsv",
        auspice_config =    "{seg}/config/auspice_config.json",
        colors =            "config/colors.tsv",
        clades =            "{seg}/config/clades_genome.tsv",
        regions=            "config/geo_regions.tsv",
        meta_public=        "data/meta_public.tsv",
        meta_collab =       "data/meta_collab.tsv",
        last_updated_file = "data/date_last_updated.txt",
        local_accn_file =   "data/local_accn.txt",
        SEQUENCES =         "data/fetch/sequences.fasta",
        METADATA =          "data/fetch/metadata.tsv",

files = rules.files.input
##############################

# Expand augur JSON paths
rule all:
    input:
        augur_jsons = expand("auspice/coxsackievirus_A6_{segs}.json", segs=segments),
        meta = files.METADATA,
        seq = files.SEQUENCES

rule all_genes:
    input:
        augur_jsons = expand("auspice/coxsackievirus_A6_gene_{genes}.json", genes=GENES),
        meta = files.METADATA,
        seq = files.SEQUENCES

rule next_update:
    """Final rule to generate all required JSON outputs for a monthly run"""
    input:
        expand("auspice/coxsackievirus_A6_{segs}.json", segs=segments),
        # expand("auspice/coxsackievirus_A6_gene_{genes}.json", genes=["-vp1", "-3D"]),

    threads: workflow.cores


##############################
# Download from NBCI Virus with ingest snakefile
###############################
if DOWNLOAD_INGEST==True:
    rule fetch:
        input:
            dir = "ingest"
        output:
            sequences=files.SEQUENCES,
            metadata=files.METADATA
        threads: workflow.cores
        shell:
            """
            cd {input.dir}
            snakemake --cores {threads} all
            cd ../
            """

##############################

# # This rule is very slow. Only give accessions as input where you are certain that they have GenBank metadata.
rule fetch_metadata:
    message:
        """
        Retrieving GenBank metadata for the specified accessions.
        """
    input:
        accessions="data/metadata/FRA.txt",
        config="config/config.yaml" # include symptom list and isolation source mapping
    output:
        metadata="data/metadata/FRA.tsv",
    params:
        virus="Coxsackievirus A6",
        columns = "accession	country	location	region	subgenogroup	lineage	date	collection_yr	sex	age_yrs	age_mo	diagnosis	isolation	origin	strain	doi",
        genbank_metadata="data/genbank_metadata.tsv"
    log:
        "logs/fetch_metadata.log"
    shell:
        """
        python scripts/fetch_genbank_metadata.py \
            --virus "{params.virus}" \
            --accession_file {input.accessions} \
            --output {output.metadata} \
            --genbank {params.genbank_metadata} \
            --config {input.config} \
            --columns "{params.columns}" \
            2> {log}
        """

##############################
# Change the format of the dates in the metadata
# Attention: `augur curate` only accepts iso 8 formats; please make sure that you save e.g. Excel files in the correct format
###############################

rule curate:
    message:
        """
        Cleaning up metadata with augur curate
        """
    input:
        metadata = files.meta_public,  # Path to input metadata file
        meta_collab = files.meta_collab,  # Data shared with us by collaborators
        meta_genbank = "data/genbank_metadata.tsv"
    params:
        strain_id_field=config["id_field"],
        date_fields=config["curate"]["date_fields"],
        expected_date_formats=config["curate"]["expected_date_formats"],
        tmp = temp("temp/merged_meta.tsv"),  # Final output file for publications metadata
        meta = temp("temp/curated_meta.tsv")  # Final output file for publications metadata
    output:
        meta="data/curated/all_meta.tsv"  # Final merged output file
    shell:
        """
        mkdir -p temp

        augur curate format-dates \
            --metadata {input.metadata} \
            --date-fields {params.date_fields} \
            --no-mask-failure \
            --expected-date-formats {params.expected_date_formats} \
            --id-column {params.strain_id_field} \
            --output-metadata {params.meta}

        # Merge curated metadata
        augur merge --metadata metadata={params.meta} meta_collab={input.meta_collab} genbank={input.meta_genbank}\
            --metadata-id-columns {params.strain_id_field} \
            --output-metadata {params.tmp}

        # Normalize strings for publication metadata
        augur curate normalize-strings \
            --id-column {params.strain_id_field} \
            --metadata {params.tmp} \
        | augur curate format-dates \
            --date-fields {params.date_fields} \
            --no-mask-failure \
            --expected-date-formats {params.expected_date_formats} \
            --id-column {params.strain_id_field} \
            --output-metadata {output.meta}
        echo "Curated metadata saved to {output.meta}"
        """

##############################
# Add additional sequences
# if you have sequences that are not on NCBI Virus
###############################

rule update_sequences:
    input:
        sequences = files.SEQUENCES,
        metadata = files.METADATA,
        extra_metadata = rules.curate.output.meta
    output:
        sequences = "data/all_sequences.fasta",
    params:
        strain_id_field = config["id_field"],
        file_ending = "data/*.fas*",
        temp = temp("temp/sequences.fasta"),
        date_last_updated = files.last_updated_file,
        local_accn = files.local_accn_file,
    shell:
        """
        set -euo pipefail
        shopt -s nullglob

        mkdir -p temp
        tmpdir=$(mktemp -d temp/update_seq.XXXXXX)
        tmp="$tmpdir/merged_sequences.fasta"
        dedup="$tmpdir/dedup_sequences.fasta"
        trap 'rm -rf "$tmpdir"' EXIT

        # Concatenate ingest sequences + any additional fasta files (if present)
        cat {input.sequences} {params.file_ending} > "$tmp"

        python scripts/update_sequences.py --in_seq "$tmp" --dates {params.date_last_updated} \
            --local_accession {params.local_accn} --meta {input.metadata} --add {input.extra_metadata} \
            --ingest_seqs {input.sequences} --out_seq {params.temp}

        # Deduplicate FASTA headers (keep first occurrence) and atomically replace
        seqkit rmdup {params.temp} > {output.sequences}
        """


##############################
# BLAST
# blast fasta files for your specific proteins
# cut out your protein from fasta sequences
###############################

rule extract:
    input: 
        genbank_file = files.reference
    output: 
        extracted_fasta = "{seg}/config/reference.fasta",    
        extracted_genbank = "{seg}/config/reference.gbk",
        annotation = "{seg}/config/annotation.gff3"
    params:
        product_name = "{seg}",
        taxid = TAXID,
        # annotation = lambda wildcards: f'{wildcards.seg}/config/annotation.gff3' if wildcards.seg != "whole_genome" else ""

    shell:
        """
        python scripts/extract_gene_from_whole_genome.py \
        --genbank_file {input.genbank_file} \
        --output_fasta {output.extracted_fasta} \
        --product_name {params.product_name} \
        --output_genbank {output.extracted_genbank} \
        --taxid {params.taxid} \
        --output_gff {output.annotation}
        """

rule blast:
    input: 
        blast_db_file = rules.extract.output.extracted_fasta,  
        seqs_to_blast = rules.update_sequences.output.sequences,
    output:
        blast_out = "temp/{seg}/blast_out.csv"
    params:
        blast_db =  "temp/{seg}/blast_database"
    shell:
        """
        sed -i 's/-//g' {input.seqs_to_blast}
        makeblastdb -in {input.blast_db_file} -out {params.blast_db} -dbtype nucl
        blastn -task blastn -query {input.seqs_to_blast} -db {params.blast_db} \
        -outfmt '10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -out {output.blast_out} -evalue 0.0005
        """

rule blast_sort:
    input:
        blast_result = rules.blast.output.blast_out, # output blast (for your protein)
        input_seqs = rules.update_sequences.output.sequences,
    output:
        sequences = "{seg}/results/sequences.fasta"
        
    params:
        range = "{seg}",  # Determines which protein (or whole genome) is processed
        min_length = lambda wildcards: {"vp1": 600, "whole_genome": 6400, "P1": 2000}[wildcards.seg],  # Min length
        max_length = lambda wildcards: {"vp1": 950, "whole_genome": 8000, "P1": 2650}[wildcards.seg]  # Max length
    shell:
        """
        python scripts/blast_sort.py --blast {input.blast_result} \
            --seqs {input.input_seqs} \
            --out_seqs {output.sequences} \
            --range {params.range} \
            --min_length {params.min_length} \
            --max_length {params.max_length}
        """

##############################
# Merge all metadata files (NCBI download and own files) and clean them up
# potentially use augur merge: but not the same output can be achieved with augur
###############################

rule add_metadata:
    message:
        """
        Cleaning data in metadata
        """
    input:
        metadata = files.METADATA,
        new_data = rules.curate.output.meta,
        regions = ancient(files.regions),
        last_updated = files.last_updated_file,
        local_accn = files.local_accn_file,
    params:
        strain_id_field = config["id_field"],
    output:
        metadata = "data/all_metadata.tsv"
    shell:
        """
        python scripts/add_metadata.py \
            --input {input.metadata} \
            --add {input.new_data} \
            --local {input.local_accn} \
            --update {input.last_updated}\
            --regions {input.regions} \
            --id {params.strain_id_field} \
            --output {output.metadata}
            """

## Deduplicate sequences that have identical strain names and sequences
rule deduplicate:
    message:
        """
        Deduplicating sequences with identical strain names and sequences
        """
    input:
        sequences = rules.blast_sort.output.sequences,
        metadata = rules.add_metadata.output.metadata
    params:
        id_field = config["id_field"],
        threshold = 0.98 # percent identity threshold to consider sequences as duplicates
    output:
        sequences = "{seg}/results/deduplicated_sequences.fasta",
    shell:
        """
        python scripts/deduplicate.py \
            --in-sequences {input.sequences} \
            --metadata {input.metadata} \
            --id-field {params.id_field} \
            --threshold {params.threshold} \
            --out-sequences {output.sequences} 
        """


##############################
# Create an index of sequence composition for filtering & filter
###############################
rule index_sequences:
    message:
        """
        Creating an index of sequence composition for filtering
        """
    input:
        sequences = rules.deduplicate.output.sequences
    output:
        sequence_index = "{seg}/results/sequence_index.tsv"
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.sequence_index}
        """

rule filter:
    message:
        """
        Filtering to
          - {params.sequences_per_group} sequence(s) per {params.group_by!s}
          - from {params.min_date} onwards
          - excluding strains in {input.exclude}
        """
    input:
        sequences = rules.deduplicate.output.sequences,
        sequence_index = rules.index_sequences.output.sequence_index,
        metadata = rules.add_metadata.output.metadata,
        exclude = files.dropped_strains,
        include = files.incl_strains,
    output:
        sequences = "{seg}/results/filtered.fasta",
        reason ="{seg}/results/reasons.tsv",
    log:
        "logs/filter.{seg}.log"
    params:
        group_by = "country year",
        sequences_per_group = 500, # set lower if you want to have a max sequences per group
        strain_id_field= config["id_field"],
        min_date = 1949,  # G-10 was collected in 1952
        min_length = lambda wildcards: {"vp1": 600, "whole_genome": 6400, "P1": 2000}[wildcards.seg], # to be safe
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --exclude {input.exclude} \
            --include {input.include} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --min-date {params.min_date} \
            --min-length {params.min_length} \
            --output-sequences {output.sequences}\
            --output-log {output.reason} \
            >> {log} 2>&1

        echo "Filtered sequences saved to {output.sequences}"
        """

##############################
# Reference sequence &
# Alignment
###############################
rule align: 
    message:
        """
        Segment: {wildcards.seg}
        Aligning sequences to {input.reference} using Nextclade run.
        """
    input:
        gff_reference = files.gff_reference,
        sequences = rules.filter.output.sequences,
        reference = rules.extract.output.extracted_fasta,
    output:
        alignment = "{seg}/results/aligned.fasta",
        tsv = "{seg}/results/nextclade.tsv",    
    params:
        penalty_gap_extend = config["align"]["penalty_gap_extend"],
        penalty_gap_open = config["align"]["penalty_gap_open"],
        penalty_gap_open_in_frame = config["align"]["penalty_gap_open_in_frame"],
        penalty_gap_open_out_of_frame = config["align"]["penalty_gap_open_out_of_frame"],
        kmer_length = config["align"]["kmer_length"],
        kmer_distance = config["align"]["kmer_distance"],
        min_match_length = config["align"]["min_match_length"],
        allowed_mismatches = config["align"]["allowed_mismatches"],
        min_length = config["align"]["min_length"]
    threads: workflow.cores
    shell:
        """
        nextclade3 run \
        -j {threads} \
        {input.sequences} \
        --input-ref {input.reference} \
        --input-annotation {input.gff_reference} \
        --penalty-gap-open {params.penalty_gap_open} \
        --penalty-gap-extend {params.penalty_gap_extend} \
        --penalty-gap-open-in-frame {params.penalty_gap_open_in_frame} \
        --penalty-gap-open-out-of-frame {params.penalty_gap_open_out_of_frame} \
        --kmer-length {params.kmer_length} \
        --kmer-distance {params.kmer_distance} \
        --min-match-length {params.min_match_length} \
        --allowed-mismatches {params.allowed_mismatches} \
        --min-length {params.min_length} \
        --include-reference false \
        --output-tsv {output.tsv} \
        --output-translations "{wildcards.seg}/results/translations/cds_{{cds}}.translation.fasta" \
        --output-fasta {output.alignment}
        """

# potentially add one-by-one genes
# use wildcards
# rule sub_alignments:
#     input:
#         alignment=rules.align.output.alignment,
#         reference=files.reference
#     output:
#         alignment = "{seg}/results/aligned{gene}.fasta"
#     benchmark:
#         "benchmark/sub_alignments.{seg}{gene}.log"
#     run:
#         from Bio import SeqIO
#         from Bio.Seq import Seq

#         real_gene = wildcards.gene.replace("-", "", 1)

#         # Extract boundaries from the reference GenBank file
#         gene_boundaries = {}
#         with open(input.reference) as handle:
#             for record in SeqIO.parse(handle, "genbank"):
#                 for feature in record.features:
#                     if feature.type == "CDS" and 'Name' in feature.qualifiers:
#                         product = feature.qualifiers['Name'][0].upper()
#                         if product == real_gene.upper():
#                             # Corrected: Use .start and .end directly
#                             gene_boundaries[product] = (feature.location.start, feature.location.end)

#         if real_gene.upper() not in gene_boundaries:
#             raise ValueError(f"Gene {real_gene} not found in reference file.")

#         b = gene_boundaries[real_gene.upper()]

#         alignment = SeqIO.parse(input.alignment, "fasta")
#         with open(output.alignment, "w") as oh:
#             for record in alignment:
#                 sequence = Seq(record.seq)
#                 gene_keep = sequence[b[0]:b[1]]
#                 if set(gene_keep) in [{"N"}, {"-"}, set()]:
#                     continue  # Skip sequences that are entirely masked
#                 sequence = len(sequence) * "-"
#                 sequence = sequence[:b[0]] + gene_keep + sequence[b[1]:]
#                 record.seq = Seq(sequence)
#                 SeqIO.write(record, oh, "fasta")

##############################
# Tree building
###############################
rule tree:
    message:
        """
        Segment: {wildcards.seg} {wildcards.gene}
        Creating a maximum likelihood tree
        """
    input:
        alignment = rules.align.output.alignment,
        # alignment = rules.sub_alignments.output.alignment,
    output:
        # tree = "{seg}/results/tree_raw.nwk"
        tree = "{seg}/results/tree_raw{gene}.nwk"
    threads: workflow.cores    
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --nthreads {threads} \
            --output {output.tree}
        """

##############################
# Refine tree &
# Ancestral sequence reconstruction
# & Translation
###############################

rule refine:
    message:
        """
        Segment: {wildcards.seg} {wildcards.gene}
        Refining tree by rerooting and resolving polytomies
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
          - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output.alignment,
        # alignment = rules.sub_alignments.output.alignment,
        metadata =  rules.add_metadata.output.metadata,
    output:
        # tree = "{seg}/results/tree.nwk",
        # node_data = "{seg}/results/branch_lengths.json"
        tree = "{seg}/results/tree{gene}.nwk",
        node_data = "{seg}/results/branch_lengths{gene}.json"
    params:
        coalescent = "opt",
        date_inference = "marginal",
        clock_filter_iqd = lambda w: 11 if getattr(w, "seg", "") == "vp1" else 3,
        strain_id_field = config["id_field"],
        clock_rate = 0.0039, # estimated with clockor: VP1 = 3.882 x 10^-3, WHOLE-GENOME = 4.033 x 10^-3
        clock_std_dev = 0.0015
        # clock_rate_string = lambda wildcards: f"--clock-rate 0.004 --clock-std-dev 0.0015" if wildcards.gene or wildcards.quart else ""
    log:
        "logs/refine.{seg}{gene}.log"
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent {params.coalescent} \
            --date-confidence \
            --stochastic-resolve \
            --clock-rate {params.clock_rate}\
            --clock-std-dev {params.clock_std_dev} \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd} \
            2>&1 | (grep -i "pruning leaf" || cat > /dev/null) > {log}

        echo "Refined tree saved to {output.tree}"
        """
        #            


rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        # alignment = rules.sub_alignments.output.alignment,
        alignment = rules.align.output.alignment,
        annotation = rules.extract.output.extracted_genbank,
    output:
        node_data = "{seg}/results/muts{gene}.json",
    params:
        inference = "joint",
        genes = lambda wildcards: (
            CODING_GENES if wildcards.seg == "whole_genome" 
            else [wildcards.seg.upper()] if not wildcards.gene
            else []
        ),
        translation_template= r"{seg}/results/translations/cds_%GENE.translation.fasta",
        output_translation_template=r"{seg}/results/translations/cds_%GENE.ancestral.fasta",
        root = "{seg}/results/ancestral_sequences.fasta",

    run:
        # Check if this is for a specific gene (wildcards.gene is not empty)
        if wildcards.gene != "":
            # Running for a specific gene
            shell("""
                augur ancestral \
                --tree {input.tree} \
                --alignment {input.alignment} \
                --output-node-data {output.node_data} \
                --keep-ambiguous \
                --inference {params.inference}
            """)
        else:
            # Running for whole genome with translation
            shell("""
                augur ancestral \
                --tree {input.tree} \
                --alignment {input.alignment} \
                --annotation {input.annotation} \
                --genes {params.genes} \
                --translations {params.translation_template} \
                --output-node-data {output.node_data} \
                --output-translations {params.output_translation_template} \
                --output-sequences {params.root} \
                --skip-validation
            """)

            # --root-sequence {input.annotation} \  -> assigns mutations to the root relative to the reference, not wanted here


##############################
# Clade assignment
###############################

rule clades: 
    message: "Assigning clades according to nucleotide mutations"
    input:
        tree=rules.refine.output.tree,
        muts = rules.ancestral.output.node_data,
        clades = files.clades #"vp1/config/vp1_clades.tsv" 
    output:
        # clade_data = "{seg}/results/clades.json"
        clade_data = "{seg}/results/clades{gene}.json"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.muts} \
            --clades {input.clades} \
            --output-node-data {output.clade_data}
        """

rule traits:
    message: "Inferring ancestral traits for {params.traits!s}"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.add_metadata.output.metadata
    output:
        # node_data = "{seg}/results/traits.json"
        node_data = "{seg}/results/traits{gene}.json",
    params:
        traits = "country",
        strain_id_field= config["id_field"]
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --output-node-data {output.node_data} \
            --columns {params.traits} \
            --confidence
        """

rule clade_published:
    message: "Assigning clades from publications"
    input:
        metadata = rules.add_metadata.output.metadata,
        subgenotypes = "data/clades_vp1.tsv",
        alignment="vp1/results/aligned.fasta"
    params:
        strain_id_field= config["id_field"]
    output:
        final_metadata = "data/final_metadata.tsv"
    run:
        import pandas as pd
        from Bio import SeqIO
        import numpy as np

        # Load the input data files
        metadata_df = pd.read_csv(input.metadata, sep="\t")
        subgenotypes_df = pd.read_csv(input.subgenotypes, sep="\t")

        # Merge the dataframes on the specified column
        merged_df = pd.merge(metadata_df, subgenotypes_df, on=params.strain_id_field, how="left")

        # Read alignment
        seqs = list(SeqIO.parse(input.alignment, "fasta"))

        # Calculate VP1 lengths for each sequence without gaps ("-") and Ns
        ids = [record.id for record in seqs]
        vp1_lengths = [len(record.seq.replace("N", "").replace("-", "")) for record in seqs]

        # Create a DataFrame with sequence IDs and VP1 lengths
        len_df = pd.DataFrame({"accession": ids, "length_vp1": vp1_lengths})

        # Add the length to the metadata
        merged_df = pd.merge(merged_df, len_df, left_on=params.strain_id_field, right_on="accession", how="left")

        # add url with genbank accession
        merged_df['url'] = "https://www.ncbi.nlm.nih.gov/nuccore/" + merged_df['accession']

        merged_df.rename(columns={"has_age":"Age available"}, inplace=True)
        merged_df.rename(columns={"has_diagnosis":"Diagnosis available"}, inplace=True)
        
        # Save the merged dataframe to the output file
        merged_df.to_csv(output.final_metadata, sep="\t", index=False)

#########################
#  EXPORT
#########################

rule export:
    message: "Creating auspice JSONs"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.clade_published.output.final_metadata,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.node_data,
        colors = files.colors,
        lat_longs = files.lat_longs,
        auspice_config = files.auspice_config,
        clades = rules.clades.output.clade_data,
        muts = rules.ancestral.output.node_data,
    params:
        strain_id_field= config["id_field"],
        muts_flag = lambda wildcards: "" if wildcards.gene else f"{wildcards.seg}/results/muts.json",
    output:
        auspice_json = "auspice/coxsackievirus_A6_{seg}{gene}.json"        
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --node-data {input.branch_lengths} {input.traits} {params.muts_flag} {input.clades} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json}
        """


# ##############################

rule rename_whole_genome:
    message: 
        "Rename whole-genome built"
    input: 
        json="auspice/coxsackievirus_A6_whole_genome.json"
    output:
        json="auspice/coxsackievirus_A6_whole-genome.json" # easier view in auspice
    shell:
        """
        mv {input.json} {output.json}
        """

rule rename_genes:
    message: 
        "Rename the single genome builts"
    input: 
        json="auspice/coxsackievirus_A6_whole_genome{gene}.json"
    output:
        json="auspice/coxsackievirus_A6_gene_{gene}.json" # easier view in auspice
    shell:
        """
        mv {input.json} {output.json}
        """

rule clean:
    message: 
        """
        Removing previous results and temporary files in the following directories:
        {params.targets}
        """
    params:
        targets = [
            "*/results/*",
            "auspice/*.json",
            "ingest/data/*.*",
            "temp/*",
            "logs/*",
            "benchmark/*",
            "data/fetch/*",
            "data/all_*.*",
            "data/curated/*"
        ]
    shell: 
        """
        for dir in {params.targets}; do
            rm -rf "$dir" 2>/dev/null || true
        done
        """

rule upload: ## make sure you're logged in to Nextstrain
    message: "Uploading auspice JSONs to Nextstrain"
    input:
        jsons = expand("auspice/coxsackievirus_A6_{segs}.json", segs=segments)
    params:
        remote_group=REMOTE_GROUP,
        date=UPLOAD_DATE,
        USERNAME=os.getenv("NEXTSTRAIN_USERNAME"),

    shell:
        """
        nextstrain login --username {params.USERNAME}
        nextstrain remote upload \
            nextstrain.org/groups/{params.remote_group}/ \
            {input.jsons}
        nextstrain logout
        mkdir -p auspice/{params.date}
        cp {input.jsons} auspice/{params.date}/
        """