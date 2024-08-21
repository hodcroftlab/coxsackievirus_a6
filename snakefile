###############
# Snakemake execution templates:

# To run a default VP1 run(<600bp):
# snakemake  vp1/auspice/cva6_vp1.json --cores 1

# To run a default whole genome run ( <6400bp):
# snakemake whole_genome/auspice/cva6_whole_genome.json --cores 1

###############
# Define segments to analyze
segments = ['vp1', 'whole_genome']

# Expand augur JSON paths
augur_jsons = expand("auspice/cva6_{seg}.json", seg=segments)

rule all:
    input:
        augur_jsons

# Rule to handle configuration files
rule files:
    input:
        sequence_length =   "{seg}",
        dropped_strains =   "config/dropped_strains.txt",
        reference =         "{seg}/config/reference_sequence.gb",
        lat_longs =         "config/lat_longs.tsv",
        auspice_config =    "{seg}/config/auspice_config.json",
        colors =            "config/colors.tsv",
        clades =            "{seg}/config/clades_genome.tsv",
        regions=            "config/geo_regions.tsv",
        extended_metafile=  "data/assign_publications_corrected.tsv"

files = rules.files.input

##############################
# Download from NBCI Virus with ingest snakefile
###############################

rule fetch:
    input:
        dir = "ingest"
    output:
        sequences="data/sequences.fasta",
        metadata="data/metadata.tsv"
    params:
        seq="ingest/data/sequences.fasta",
        meta="ingest/data/metadata.tsv"
    shell:
        """
        cd {input.dir} 
        snakemake --cores 9 all
        cd ../
        cp -u {params.seq} {output.sequences}
        cp -u {params.meta} {output.metadata}
        """

##############################
# Update strain names
###############################

rule update_strain_names:
    message:
        """
        Updating strain information in metadata.
        """
    input:
        file_in = "data/metadata.tsv"
    output:
        file_out = "data/updated_strain_names.tsv"
    shell:
        """
        time bash scripts/update_strain.sh {input.file_in} {output.file_out}
        """

##############################
# Add additional sequences
# if you have sequences that are not on NCBI Virus
###############################

rule update_sequences:
    input:
        sequences = "data/sequences.fasta",
        date_last_updated = "data/date_last_updated.txt",
        local_accn = "data/local_accn.txt",
        add_metadata = files.extended_metafile
    output:
        sequences = "data/sequences_added.fasta"
    params:
        file_ending = "data/*.fas*",
        temp = "data/temp_sequences_added.fasta"
    shell:
        """
        touch {params.temp} && rm {params.temp}
        cat {params.file_ending} > {params.temp}
        python scripts/update_sequences.py --in_seq {params.temp} --out_seq {output.sequences} --dates {input.date_last_updated} --local_accession {input.local_accn} --add {input.add_metadata}
        rm {params.temp}
        awk '/^>/{{if (seen[$1]++ == 0) print; next}} !/^>/{{print}}' {output.sequences} > {params.temp} && mv {params.temp} {output.sequences}
        """


##############################
# BLAST
# blast fasta files for vp1 
###############################

rule blast:
    input: 
        blast_db_file = "data/references/reference_vp1_blast.fasta",
        seqs_to_blast = rules.update_sequences.output.sequences
    output:
        blast_out = "vp1/temp/blast_out.csv"
    params:
        blast_db = "vp1/temp/entero_db_vp1"
    shell:
        """
        sed -i 's/-//g' {input.seqs_to_blast}
        makeblastdb -in {input.blast_db_file} -out {params.blast_db} -dbtype nucl
        blastn -task blastn -query {input.seqs_to_blast} -db {params.blast_db} -outfmt '10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -out {output.blast_out} -evalue 0.0005
        """

rule blast_sort:
    input:
        blast_result = rules.blast.output.blast_out, # output blast (vp1)
        input_seqs = rules.update_sequences.output.sequences
    output:
        sequences = "{seg}/results/sequences.fasta"
        
    params:
        matchLen = 300,
        range="{seg}"
    shell:
        """
        python scripts/blast_sort.py --blast {input.blast_result} \
            --seqs {input.input_seqs} \
            --out_seqs {output.sequences} \
            --range {params.range}
        """

##############################
# Change the format of the dates in the metadata
# Attention: ```augur curate``` only accepts iso 8 formats; please make sure that you save e.g. Excel files in the correct format
###############################

rule curate_meta_dates:
    message:
        """
        Cleaning up metadata with augur curate
        """
    input:
        metadata =  files.extended_metafile, #files.extended_metafile
        
    params:
        strain_id_field= "strain",
        date_column = 'date',
        format=['%Y', '%m.%Y', '%d.%m.%Y'],
        temp_metadata="data/temp_curated.tsv"  # Temporary file
    output:
        metadata="data/assign_publications_curated.tsv"  # Final output file
    shell:
        """
        augur curate normalize-strings --metadata {input.metadata} \
            --id-column {params.strain_id_field} \
            --output-metadata {params.temp_metadata} 
        augur curate format-dates \
            --metadata {params.temp_metadata} \
            --date-fields {params.date_column} \
            --expected-date-formats {params.format} \
            --id-column {params.strain_id_field} \
            --output-metadata {output.metadata}
        rm {params.temp_metadata}
        """

##############################
# Merge all metadata files (NCBI download and own files)
###############################

rule add_metadata:
    message:
        """
        Cleaning dates in metadata
        """
    input:
        metadata="data/metadata.tsv",
        new_data=rules.curate_meta_dates.output.metadata,
        regions=ancient(files.regions),
        local_accn="data/local_accn.txt",
        renamed_strains=rules.update_strain_names.output.file_out
    params:
        strain_id_field="accession"
    output:
        metadata="data/final_metadata.tsv"
    shell:
        """
        python scripts/add_metadata.py \
            --input {input.metadata} \
            --add {input.new_data} \
            --rename {input.renamed_strains} \
            --local {input.local_accn} \
            --regions {input.regions} \
            --id {params.strain_id_field} \
            --output {output.metadata}
        """

##############################
# Rest of the augur pipeline
###############################

rule index_sequences:
    message:
        """
        Creating an index of sequence composition for filtering
        """
    input:
        sequences = rules.blast_sort.output.sequences
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
        sequences = rules.blast_sort.output.sequences,
        sequence_index = rules.index_sequences.output.sequence_index,
        metadata = rules.add_metadata.output.metadata,
        exclude = files.dropped_strains
    output:
        sequences = "{seg}/results/filtered.fasta"
    params:
        group_by = "country",
        sequences_per_group = 15000, # set lower if you want to have a max sequences per group
        strain_id_field= "accession",
        min_date = 1950  # Gdula was collected in 1949
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --exclude {input.exclude} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --min-date {params.min_date} \
            --output {output.sequences}
        """

rule reference_gb_to_fasta:
    message:
        """
        Converting reference sequence from genbank to fasta format
        """
    input:
        reference = files.reference

    output:
        reference = "{seg}/results/reference_sequence.fasta"
    shell:
        """
        python scripts/reference_genbank_to_fasta.py \
            --input {input.reference} \
            --output {output.reference}
        """

rule align: 
    message:
        """
        Aligning sequences to {input.reference} using Nextalign
        """
    input:
        sequences = rules.filter.output.sequences,
        reference = rules.reference_gb_to_fasta.output.reference
    output:
        alignment = "{seg}/results/aligned.fasta"

    params:
        nuc_mismatch_all = 10,
        nuc_seed_length = 30
    shell:
        """
        nextclade run \
        {input.sequences}  \
        --input-ref {input.reference}\
        --allowed-mismatches {params.nuc_mismatch_all} \
        --min-length {params.nuc_seed_length} \
        --include-reference false \
        --output-fasta {output.alignment} 
        """

rule fix_align_codon:
    input:
        sequences = rules.align.output.alignment
    output:
        alignment = "{seg}/results/aligned_fixed.fasta"
    shell:
        """
        Rscript scripts/fixAlignmentGaps.R {input.sequences} {output.alignment}
        """

# potentially add one-by-one genes
# use wildcards

rule tree:
    message:
        """
        Creating a maximum likelihood tree
        """
    input:
        alignment = rules.fix_align_codon.output.alignment
    output:
        tree = "{seg}/results/tree_raw.nwk"
    threads: 9
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --nthreads {threads} \
            --output {output.tree}
        """

rule refine:
    message:
        """
        Refining tree by rerooting and resolving polytomies
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
          - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.fix_align_codon.output.alignment,
        metadata =  rules.add_metadata.output.metadata,
        reference = rules.reference_gb_to_fasta.output.reference
    output:
        tree = "{seg}/results/tree.nwk",
        node_data = "{seg}/results/branch_lengths.json"
    params:
        coalescent = "opt",
        date_inference = "marginal",
        clock_filter_iqd = 6, # was 3
        strain_id_field ="accession",
        clock_rate = 0.004, # leave it empty for estimation?
        clock_std_dev = 0.0015
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
            --clock-rate {params.clock_rate} \
            --clock-std-dev {params.clock_std_dev} \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd}
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output.alignment,
        reference = rules.reference_gb_to_fasta.output.reference
    output:
        node_data = "{seg}/results/nt_muts.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = files.reference
    output:
        node_data = "{seg}/results/aa_muts.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data}
        """

rule clades: 
    message: "Assigning clades according to nucleotide mutations"
    input:
        tree=rules.refine.output.tree,
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = files.clades
    output:
        clade_data = "{seg}/results/clades.json"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {output.clade_data}
        """

rule traits:
    message: "Inferring ancestral traits for {params.traits!s}"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.add_metadata.output.metadata
    output:
        node_data = "{seg}/results/traits.json",
    params:
        traits = "country",
        strain_id_field= "accession"
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

rule export:
    message: "Creating auspice JSONs"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.add_metadata.output.metadata,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        clades = rules.clades.output.clade_data,
        colors = files.colors,
        lat_longs = files.lat_longs,
        auspice_config = files.auspice_config
    params:
        strain_id_field= "accession"
    output:
        auspice_json = "auspice/cva6_{seg}-accession.json"
        
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} {input.clades} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json}
        """

##############################
# Change from accession to strain name view in tree
###############################

rule rename_json:
    input:
        auspice_json= rules.export.output.auspice_json,
        metadata = rules.add_metadata.output.metadata,
    output:
        auspice_json="auspice/cva6_{seg}.json"
    params:
        strain_id_field="accession",
        display_strain_field= "strain"
    shell:
        """
        python3 scripts/set_final_strain_name.py --metadata {input.metadata} \
                --metadata-id-columns {params.strain_id_field} \
                --input-auspice-json {input.auspice_json} \
                --display-strain-name {params.display_strain_field} \
                --output {output.auspice_json}
        """

rule clean:
    message: "Removing directories: {params}"
    params:
        # "results ",
        "auspice"
    shell:
        "rm -rfv {params}"


rule rename_whole_genome:
    message: 
        "Rename whole-genome built"
    input: 
        json="auspice/cva6_whole_genome.json"
    output:
        json="auspice/cva6_whole-genome.json" # easier view in auspice
    shell:
        """
        mv {input.json} {output.json}
        """
