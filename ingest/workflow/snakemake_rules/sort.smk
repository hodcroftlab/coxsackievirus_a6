"""
This part of the workflow handles sorting downloaded sequences and metadata
by aligning them to reference sequences.

It produces output files as

    metadata = "data/metadata.tsv"
    sequences = "data/sequences.fasta"

"""

rule sort:
    input:
        sequences = rules.curate.output.sequences
    output:
        "data/sequences.fasta",
    shell:
        '''
        seqkit rmdup {input.sequences} > {output}
        '''

rule metadata:
    input:
        metadata = rules.curate.output.metadata,
        sequences = "data/sequences.fasta"
    output:
        metadata = "data/metadata_raw.tsv"
    run:
        import pandas as pd
        from Bio import SeqIO

        strains = [s.id for s in SeqIO.parse(input.sequences, 'fasta')]
        d = pd.read_csv(input.metadata, sep='\t', index_col='accession').loc[strains].drop_duplicates()
        d.to_csv(output.metadata, sep='\t')

rule nextclade:
    input:
        sequences = "data/sequences.fasta",
        ref = "data/references/reference.fasta"
    output:
        nextclade = "data/references/nextclade.tsv"
    params:
        dataset = "data/references/",
        output_columns = "seqName clade qc.overallScore qc.overallStatus alignmentScore  alignmentStart  alignmentEnd  coverage dynamic"
    threads: 8
    shell:
        """
        nextclade3 run -D {params.dataset}  -j {threads} \
                          --output-columns-selection {params.output_columns} \
                          --output-tsv {output.nextclade} \
                          {input.sequences}
        """

rule extend_metadata: 
    input:
        nextclade = "data/references/nextclade.tsv",
        metadata = "data/metadata_raw.tsv"
    output:
        metadata = "data/metadata.tsv"
    shell:
        """
        python3 bin/extend-metadata.py --metadata {input.metadata} \
                                       --id-field accession \
                                       --nextclade {input.nextclade} \
                                       --output {output.metadata}
        """

