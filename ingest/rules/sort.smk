"""
This part of the workflow handles sorting downloaded sequences and metadata
by aligning them to reference sequences.

It produces output files as

    metadata = "results/metadata.tsv"
    sequences = "results/sequences.fasta"

"""

rule sort:
    input:
        sequences = rules.curate.output.sequences
    output:
        "results/sequences.fasta",
    shell:
        '''
        seqkit rmdup {input.sequences} > {output}
        '''

rule metadata:
    input:
        metadata = rules.subset_metadata.output.subset_metadata,
        sequences = "results/sequences.fasta"
    output:
        metadata = "results/metadata_raw.tsv"
    run:
        import pandas as pd
        from Bio import SeqIO

        strains = [s.id for s in SeqIO.parse(input.sequences, 'fasta')]
        d = pd.read_csv(input.metadata, sep='\t', index_col='accession').loc[strains].drop_duplicates()
        d.to_csv(output.metadata, sep='\t')

rule nextclade:
    input:
        sequences = "results/sequences.fasta",
        ref = "data/references/reference.fasta"
    output:
        nextclade = "results/nextclade.tsv"
    params:
        dataset = "data/references/",
        min_seed = "0.25"
    threads: workflow.cores
    log:
        "logs/run_nextclade.txt"
    benchmark:
        "benchmarks/run_nextclade.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        nextclade3 run \
            -D {params.dataset:q} -j {threads} \
            --min-seed-cover {params.min_seed:q} \
            --output-tsv {output.nextclade:q} \
            {input.sequences:q} 
        """

rule nextclade_metadata:
    input:
        nextclade="results/nextclade.tsv",
    output:
        nextclade_metadata=temp("results/nextclade_metadata.tsv"),
    params:
        nextclade_id_field=config["nextclade"]["id_field"],
        nextclade_field_map=[f"{old}={new}" for old, new in config["nextclade"]["field_map"].items()],
        nextclade_fields=",".join(config["nextclade"]["field_map"].values()),
    benchmark:
        "benchmarks/nextclade_metadata.txt"
    log:
        "logs/nextclade_metadata.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        augur curate rename \
            --metadata {input.nextclade:q} \
            --id-column {params.nextclade_id_field:q} \
            --field-map {params.nextclade_field_map:q} \
            --output-metadata - \
        | csvtk cut -t --fields {params.nextclade_fields:q} \
        > {output.nextclade_metadata:q}
        """


rule join_metadata_and_nextclade:
    input:
        metadata="data/subset_metadata.tsv",
        nextclade_metadata="results/nextclade_metadata.tsv",
    output:
        metadata = "results/metadata.tsv",
    params:
        metadata_id_field=config["curate"]["output_id_field"],
        nextclade_id_field=config["nextclade"]["id_field"],
    benchmark:
        "benchmarks/join_metadata_and_nextclade.txt"
    log:
        "logs/join_metadata_and_nextclade.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        augur merge \
            --metadata \
                metadata={input.metadata:q} \
                nextclade={input.nextclade_metadata:q} \
            --metadata-id-columns \
                metadata={params.metadata_id_field:q} \
                nextclade={params.nextclade_id_field:q} \
            --output-metadata {output.metadata:q} \
            --no-source-columns
        """

