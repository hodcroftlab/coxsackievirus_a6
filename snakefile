###############
# Snakemake execution templates:

# To run a default VP1 run(<600bp):
# snakemake  vp1/auspice/cva6_vp1.json --cores 1

# To run a default whole genome run ( <6400bp):
# snakemake whole_genome/auspice/cva6_whole_genome.json --cores 1

###############
wildcard_constraints:
    seg="vp1|whole_genome",
    gene="|-vp4|-vp2|-vp3|-vp1|-2A|-2B|-2C|-3A|-3B|-3C|-3D"
   
#     #from: https://bitbucket.org/snakemake/snakemake/issues/910/empty-wildcard-assignment-works-only-if

# Define segments to analyze
segments = ['vp1', 'whole_genome']
GENES=["-vp4", "-vp2", "-vp3", "-vp1", "-2A", "-2B", "-2C", "-3A", "-3B", "-3C", "-3D"]

# Expand augur JSON paths
rule all:
    input:
        augur_jsons = expand("auspice/cva6_{segs}.json", segs=segments)

rule all_genes:
    input:
        augur_jsons = expand("auspice/cva6_whole_genome{genes}.json", genes=GENES)


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
        meta=               "data/metadata.tsv",
        extended_metafile=  "data/assign_publications_corrected.tsv",
        last_updated_file = "data/date_last_updated.txt",
        local_accn_file =   "data/local_accn.txt"

files = rules.files.input

##############################
# Download from NBCI Virus with ingest snakefile
###############################

rule fetch:
    input:
        dir = "ingest"
    output:
        sequences="data/sequences.fasta",
        metadata=files.meta
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

# rule update_strain_names:
#     message:
#         """
#         Updating strain information in metadata.
#         """
#     input:
#         file_in = files.meta
#     output:
#         file_out = "data/updated_strain_names.tsv"
#     shell:
#         """
#         time bash scripts/update_strain.sh {input.file_in} {output.file_out}
#         """

##############################
# Add additional sequences
# if you have sequences that are not on NCBI Virus
###############################

rule update_sequences:
    input:
        sequences = "data/sequences.fasta",
        metadata=files.meta,
        add_metadata = files.extended_metafile
    output:
        sequences = "data/sequences_added.fasta"
    params:
        file_ending = "data/*.fas*",
        temp = "data/temp_sequences_added.fasta",
        date_last_updated = files.last_updated_file,
        local_accn = files.local_accn_file,
    shell:
        """
        touch {params.temp} && rm {params.temp}
        cat {params.file_ending} > {params.temp}
        python scripts/update_sequences.py --in_seq {params.temp} --out_seq {output.sequences} --dates {params.date_last_updated} --local_accession {params.local_accn} --meta {input.metadata} --add {input.add_metadata}
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
        metadata=files.extended_metafile,  # Path to input metadata file
        genbank_meta="data/metadata/genbank_metadata_additional.tsv"  # Generated with bin/extract_genbank_metadata.py
    params:
        strain_id_field="accession",
        date_column="date",
        format=['%Y', '%m.%Y', '%d.%m.%Y', "%b-%Y", "%d-%b-%Y","%Y-%m-%d"],
        temp_metadata="data/temp_curated.tsv"  # Temporary file
    output:
        metadata="data/assign_publications_curated.tsv",  # Final output file for metadata
        genbank_meta="data/metadata/genbank_metadata_curated.tsv",  # Curated genbank metadata
        final_metadata="data/assign_publications_corr_fetched.tsv"  # Final merged output file
    shell:
        """
        # Normalize strings for metadata
        augur curate normalize-strings --metadata {input.metadata} \
            --id-column {params.strain_id_field} \
            --output-metadata {params.temp_metadata}

        # Format dates for metadata
        augur curate format-dates \
            --metadata {params.temp_metadata} \
            --date-fields {params.date_column} \
            --no-mask-failure \
            --expected-date-formats {params.format} \
            --id-column {params.strain_id_field} \
            --output-metadata {output.metadata}
        
        # Remove temporary file
        rm {params.temp_metadata}

        # Normalize strings for genbank metadata
        augur curate normalize-strings --metadata {input.genbank_meta} \
            --id-column {params.strain_id_field} \
            --output-metadata {params.temp_metadata}

        # Format dates for genbank metadata
        augur curate format-dates \
            --metadata {params.temp_metadata} \
            --date-fields {params.date_column} \
            --no-mask-failure \
            --expected-date-formats {params.format} \
            --id-column {params.strain_id_field} \
            --output-metadata {output.genbank_meta}
        
        # Remove temporary file
        rm {params.temp_metadata}

        # Merge curated metadata
        augur merge --metadata meta={output.metadata} genbank_meta={output.genbank_meta}\
            --metadata-id-columns {params.strain_id_field} \
            --output-metadata {output.final_metadata}
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
        metadata=files.meta,
        new_data=rules.curate_meta_dates.output.final_metadata,
        regions=ancient(files.regions),
        local_accn=files.local_accn_file,
        last_updated=files.last_updated_file,
        renamed_strains="data/updated_strain_names.tsv"
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
            --update {input.last_updated}\
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
    run:
        from Bio import SeqIO 
        SeqIO.convert(input.reference, "genbank", output.reference, "fasta")

        # """
        # python scripts/reference_genbank_to_fasta.py \
        #     --input {input.reference} \
        #     --output {output.reference}
        # """

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
rule sub_alignments:
    input:
        alignment=rules.fix_align_codon.output.alignment,
        reference=files.reference
    output:
        # alignment = "{seg}/results/aligned.fasta"
        alignment = "{seg}/results/aligned_fixed{gene}.fasta"
    run:
        from Bio import SeqIO
        from Bio.Seq import Seq

        real_gene = wildcards.gene.replace("-", "", 1)

        # Extract boundaries from the reference GenBank file
        gene_boundaries = {}
        with open(input.reference) as handle:
            for record in SeqIO.parse(handle, "genbank"):
                for feature in record.features:
                    if feature.type == "CDS" and 'Name' in feature.qualifiers:
                        product = feature.qualifiers['Name'][0].upper()
                        if product == real_gene.upper():
                            # Corrected: Use .start and .end directly
                            gene_boundaries[product] = (feature.location.start, feature.location.end)

        if real_gene.upper() not in gene_boundaries:
            raise ValueError(f"Gene {real_gene} not found in reference file.")

        b = gene_boundaries[real_gene.upper()]

        alignment = SeqIO.parse(input.alignment, "fasta")
        with open(output.alignment, "w") as oh:
            for record in alignment:
                sequence = Seq(record.seq)
                gene_keep = sequence[b[0]:b[1]]
                if set(gene_keep) == {"N"} or len(gene_keep) == 0:
                    continue  # Skip sequences that are entirely masked
                sequence = len(sequence) * "N"
                sequence = sequence[:b[0]] + gene_keep + sequence[b[1]:]
                record.seq = Seq(sequence)
                SeqIO.write(record, oh, "fasta")


rule tree:
    message:
        """
        Creating a maximum likelihood tree
        """
    input:
        # alignment = rules.fix_align_codon.output.alignment,
        alignment = rules.sub_alignments.output.alignment
    output:
        # tree = "{seg}/results/tree_raw.nwk"
        tree = "{seg}/results/tree_raw{gene}.nwk"
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
        # alignment = rules.fix_align_codon.output.alignment,
        alignment = rules.sub_alignments.output.alignment,
        metadata =  rules.add_metadata.output.metadata,
        reference = rules.reference_gb_to_fasta.output.reference
    output:
        # tree = "{seg}/results/tree.nwk",
        tree = "{seg}/results/tree{gene}.nwk",
        # node_data = "{seg}/results/branch_lengths.json"
        node_data = "{seg}/results/branch_lengths{gene}.json"
    params:
        coalescent = "opt",
        date_inference = "marginal",
        clock_filter_iqd = 3, # was 3
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
        # alignment = rules.align.output.alignment,
        alignment = rules.sub_alignments.output.alignment,
        reference = rules.reference_gb_to_fasta.output.reference
    output:
        # node_data = "{seg}/results/nt_muts.json"
        node_data = "{seg}/results/nt_muts{gene}.json"
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
        node_data = "{seg}/results/aa_muts{gene}.json"
        # node_data = "{seg}/results/aa_muts.json"
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
        # clade_data = "{seg}/results/clades.json"
        clade_data = "{seg}/results/clades{gene}.json"
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
        # node_data = "{seg}/results/traits.json"
        node_data = "{seg}/results/traits{gene}.json",
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
        auspice_json = "auspice/cva6_{seg}{gene}-accession.json"
        # auspice_json = rules.all.input.augur_jsons
        
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

# rule export_genes:
#     message: "Creating auspice JSONs"
#     input:
#         tree = rules.refine.output.tree,
#         metadata = rules.add_metadata.output.metadata,
#         branch_lengths = rules.refine.output.node_data,
#         traits = rules.traits.output.node_data,
#         nt_muts = rules.ancestral.output.node_data,
#         aa_muts = rules.translate.output.node_data,
#         clades = rules.clades.output.clade_data,
#         colors = files.colors,
#         lat_longs = files.lat_longs,
#         auspice_config = files.auspice_config
#     params:
#         strain_id_field= "accession"
#     output:
#         auspice_json = "auspice/cva6_whole_genome{gene}-accession.json"
#         # auspice_json = rules.all_genes.input.augur_jsons
        
#     shell:
#         """
#         augur export v2 \
#             --tree {input.tree} \
#             --metadata {input.metadata} \
#             --metadata-id-columns {params.strain_id_field} \
#             --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} {input.clades} \
#             --colors {input.colors} \
#             --lat-longs {input.lat_longs} \
#             --auspice-config {input.auspice_config} \
#             --output {output.auspice_json}
#         """

##############################
# Change from accession to strain name view in tree
###############################

rule rename_json:
    input:
        auspice_json= rules.export.output.auspice_json,
        metadata = rules.add_metadata.output.metadata,
    output:
        # auspice_json = rules.all.input.augur_jsons
        auspice_json="auspice/cva6_{seg}{gene}.json"
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

# rule rename_json_all:
#     input:
#         auspice_json= rules.export_genes.output.auspice_json,
#         metadata = rules.add_metadata.output.metadata,
#     output:
#         # auspice_json = rules.all_genes.input.augur_jsons
#         auspice_json="auspice/cva6_whole_genome{gene}.json"
#     params:
#         strain_id_field="accession",
#         display_strain_field= "strain"
#     shell:
#         """
#         python3 scripts/set_final_strain_name.py --metadata {input.metadata} \
#                 --metadata-id-columns {params.strain_id_field} \
#                 --input-auspice-json {input.auspice_json} \
#                 --display-strain-name {params.display_strain_field} \
#                 --output {output.auspice_json}
#         """

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
