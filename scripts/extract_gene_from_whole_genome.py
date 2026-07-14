#To run this script use: 
# python extract_gene_from_whole_genome.py --genbank_file ../ingest/data/references/nl63_full_reference.gb --output_directory ../data/references --product_names "spike protein" "membrane protein" "nucleocapsid protein"
#use  nargs='+', if want multiple protein names
import os
import argparse
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import re
# import ipdb

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genbank_file', required=True, help='Input reference file in genbank format (.gb)')
    parser.add_argument('--output_fasta', required=True, help='Output directory for FASTA files')
    parser.add_argument('--product_name',  required=True, help='Product names of proteins to extract - ensure its the same as CDS names in reference.gb file')
    parser.add_argument('--output_genbank', required=True, help='Output genbank file')
    parser.add_argument('--output_gff3', required=False, help='Output GFF3 annotation file')
    parser.add_argument('--taxid', required=True, help='NCBI Taxonomy ID for the organism')
    
    return parser.parse_args()

def clean_product_name(name):
    """Strip generic descriptor words, keeping only the specific identifier (e.g. VP1)."""
    cleaned = re.sub(r'\b(protein|polyprotein|precursor|protease)\b', '', name, flags=re.IGNORECASE)
    cleaned = re.sub(r'\s+', ' ', cleaned).strip(' _-')
    return cleaned if cleaned else name  # fallback if regex empties it

def extract_protein(genbank_file, output_fasta, product_name, output_genbank, output_gff3=None, taxid=None):
    os.makedirs(os.path.dirname(output_fasta), exist_ok=True)  
    
    if product_name.lower() == "whole_genome":
        for record in SeqIO.parse(genbank_file, "genbank"):
            with open(output_fasta, "w") as output_handle:
                output_handle.write(f">{record.id}\n{record.seq}\n")
            print(f"Whole genome sequence saved to {output_fasta}")
            
            with open(output_genbank, "w") as output_handle:
                SeqIO.write(record, output_handle, "genbank")
            print(f"Whole genome GenBank file saved to {output_genbank}")

            if output_gff3:
                create_gff3(record, output_gff3, product_name, taxid)
        return

    found = False
    for record in SeqIO.parse(genbank_file, "genbank"):
        
        locus = record.annotations.get('locus', 'Unknown_Locus')
        date = record.annotations.get('date', 'Unknown_Date')
        accession = record.id
        version = record.annotations.get('sequence_version', 'Unknown_Version')
        source = record.annotations.get('source', 'Unknown_Source')
        organism = record.annotations.get('organism', 'Unknown_Organism')
                    
        for feature in record.features:
            if feature.type == "CDS":
                product_list = feature.qualifiers.get("product", [])
                if not product_list:
                    product_list = feature.qualifiers.get("gene", [])
                if any(re.search(r'\b' + re.escape(product_name.lower()) + r'\b', p.lower()) for p in product_list):
                    nucleotide_sequence = feature.location.extract(record.seq)

                    with open(output_fasta, "w") as output_handle:
                        output_handle.write(f">{accession} {product_name}\n{nucleotide_sequence}\n")
                    print(f"Protein {product_name} saved to {output_fasta}")
                    
                    extracted_record = record[feature.location.start:feature.location.end]
                    extracted_record.id = f"{product_name}_extracted_{accession}"
                    extracted_record.description = f"{product_name} gene sequence ({organism}, {source}, {date}, Accession: {accession}, Version: {version})"

                    extracted_record.annotations['locus'] = locus
                    extracted_record.annotations['date'] = date
                    extracted_record.annotations['accession'] = accession
                    extracted_record.annotations['version'] = version
                    extracted_record.annotations['source'] = source
                    extracted_record.annotations['organism'] = organism
                    
                    # Find the original source feature
                    original_source = next((f for f in record.features if f.type == "source"), None)

                    # Use its qualifiers if available, or set default
                    if original_source:
                        new_qualifiers = original_source.qualifiers.copy()
                    else:
                        new_qualifiers = {
                            "organism": [organism],
                            "mol_type": ["genomic RNA"],
                            "country": ["Unknown"],
                            "isolate": ["Unknown"],
                            "strain": ["Unknown"],
                            "source": [source]
                        }

                    # Create a new source feature for the extracted region
                    source_feature = SeqFeature(
                        FeatureLocation(0, len(extracted_record.seq)), 
                        type="source",
                        qualifiers=new_qualifiers
                    )

                    extracted_record.features.append(source_feature)

                    with open(output_genbank, "w") as output_handle:
                        SeqIO.write(extracted_record, output_handle, "genbank")
                    print(f"Extracted gene GenBank file saved to {output_genbank}")

                    # Create GFF3 file if requested
                    if output_gff3:
                        create_gff3(extracted_record, output_gff3, product_name, taxid)
                    
                    found = True
                    return  
    
    if not found:
        print(f"Protein {product_name} not found in the GenBank file.")


def create_gff3(record, output_gff3, product_name, taxid):
    """Create a GFF3 annotation file from a SeqRecord"""
    with open(output_gff3, "w") as f:
        # Write GFF3 header
        f.write("##gff-version 3\n")
        f.write("#!gff-spec-version 1.21\n")
        f.write("#!processor NCBI annotwriter\n")
        f.write(f"##sequence-region {record.id} 1 {len(record.seq)}\n")
        if taxid:
            f.write(f"##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={taxid}\n")
        
        # Write source feature first
        source_feature = next((f for f in record.features if f.type == "source"), None)
        if source_feature:
            start = 1
            end = len(record.seq)
            strand = "+"
            
            attributes = f"ID={record.id}:1..{end}"
            if "db_xref" in source_feature.qualifiers:
                for xref in source_feature.qualifiers["db_xref"]:
                    attributes += f";Dbxref={xref}"
            if "organism" in source_feature.qualifiers:
                attributes += f";organism={source_feature.qualifiers['organism'][0]}"
            if "mol_type" in source_feature.qualifiers:
                attributes += f";mol_type={source_feature.qualifiers['mol_type'][0]}"
            if "strain" in source_feature.qualifiers:
                attributes += f";strain={source_feature.qualifiers['strain'][0]}"
            if "nat_host" in source_feature.qualifiers:
                attributes += f";nat-host={source_feature.qualifiers['nat_host'][0]}"
            
            attributes += ";gbkey=Src"
            
            f.write(f"{record.id}\tGenbank\tregion\t{start}\t{end}\t.\t{strand}\t.\t{attributes}\n")
        
        # Write CDS features
        for i, feature in enumerate(record.features):
            if feature.type == "CDS":
                # Adjust coordinates: features in extracted record are relative to the extracted region
                start = int(feature.location.start) + 1  # GFF3 is 1-based
                end = int(feature.location.end)
                strand = "+" if feature.location.strand == 1 else ("-" if feature.location.strand == -1 else ".")
                
                # Build attributes - use product name as the gene name to avoid duplicates
                raw_product = feature.qualifiers.get("product", [product_name])[0]
                product = clean_product_name(raw_product)
                
                attributes = f"Name={product};gbkey=Prot;product={product}"
                
                # Add optional qualifiers
                if "note" in feature.qualifiers:
                    attributes += f";Note={feature.qualifiers['note'][0]}"
                if "protein_id" in feature.qualifiers:
                    attributes += f";ID=id-{feature.qualifiers['protein_id'][0]}"
                else:
                    attributes += f";ID=id-{product}_{start}..{end}"
                
                f.write(f"{record.id}\tGenbank\tCDS\t{start}\t{end}\t.\t{strand}\t.\t{attributes}\n")
    
    print(f"GFF3 annotation file saved to {output_gff3}")


if __name__ == "__main__":
    args = parse_args()
    extract_protein(args.genbank_file, args.output_fasta, args.product_name, args.output_genbank, args.output_gff3, args.taxid)