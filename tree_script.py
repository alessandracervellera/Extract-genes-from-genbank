import os
import shutil
from Bio import SeqIO
import re
import argparse

parser = argparse.ArgumentParser(description="Process GenBank files to extract gene sequences, for each species.")
parser.add_argument('-i', '--input', required=True, metavar='\b', help='input folder path')
parser.add_argument('-o', '--output', required=True, metavar='\b', help='output folder path')

args = parser.parse_args()
input_dir = args.input
output_dir = args.output

# Equivalent gene names
same_genes = [
    ("nad1", "nd1"),
    ("nad2", "nd2"),
    ("nad3", "nd3"),
    ("nad4", "nd4"),
    ("nad4l", "nd4l"),
    ("nad5", "nd5"),
    ("nad6", "nd6"),
    ("cytb", "cob"),
    ("rrn12", "rrnL", "12s ribosomal RNA", "rrnL ribosomal RNA", "small subunit ribosomal RNA", "s-RNA"),
    ("rrn16", "rrnS", "16s ribosomal RNA", "rrnS ribosomal RNA", "large subunit ribosomal RNA", "l-RNA")
]

# Create a map for equivalent genes
gene_map = {}
for group in same_genes:
	for gene in group:
		gene_map[gene.lower()] = group

def get_canonical_file_path(gene_name):
	gene_group = gene_map.get(gene_name.lower())
	if gene_group:
		for alt_gene in gene_group:
			alt_gene_file = os.path.join(output_dir, f"{alt_gene}.fasta")
			if os.path.exists(alt_gene_file):
				return alt_gene_file
	return os.path.join(output_dir, f"{gene_name}.fasta")


# Create the output directory if it doesn't exist
if not os.path.exists(output_dir):
	os.makedirs(output_dir, exist_ok=True)

# Process each GenBank file
for file_name in os.listdir(input_dir):
	if file_name.endswith(".gb"):
		file_path = os.path.join(input_dir, file_name)
		base_name = os.path.splitext(os.path.basename(file_path))[0]
		inter_dir = os.path.join(input_dir, f"{base_name}_genes")

		# Remove intermediate directory if it exists and create a new one
		if os.path.exists(inter_dir):
			shutil.rmtree(inter_dir)
			os.makedirs(inter_dir, exist_ok=True)

		# Extract genes from the GenBank file
		for gene in SeqIO.parse(file_path, "genbank"):
			for feature in gene.features:
				if feature.type == "gene":
					gene_name = feature.qualifiers.get('gene', ['unknown'])[0]
					if re.match(r'^(trna|trn).*', gene_name.lower()):
						continue
					gene_sequence = feature.extract(gene.seq)
					gene_file_path = os.path.join(inter_dir, f"{gene_name}.fasta")
					if not os.path.exists(gene_file_path):
						with open(gene_file_path, "w") as gene_file:
							gene_file.write(f">{gene_name}\n{gene_sequence}\n")
				elif feature.type == "rRNA":
					rrna_name = (feature.qualifiers.get('product') or
					feature.qualifiers.get('rrna') or
					feature.qualifiers.get('gene', ['unknown']))[0]
					rrna_sequence = feature.extract(gene.seq)
					rrna_file_path = os.path.join(inter_dir, f"{rrna_name}.fasta")
					if not os.path.exists(rrna_file_path):
						with open(rrna_file_path, "w") as rrna_file:
							rrna_file.write(f">{rrna_name}\n{rrna_sequence}\n")

# Remove output directory if it exists and create a new one
if os.path.exists(output_dir):
	shutil.rmtree(output_dir)
os.makedirs(output_dir, exist_ok=True)

# Process each subdirectory containing gene sequences
for dir_name in os.listdir(input_dir):
	dir_path = os.path.join(input_dir, dir_name)
	if os.path.isdir(dir_path) and dir_name.endswith("_genes"):
		species = dir_name[:-6]  # Extract species name by removing "_genes"

		# Iterate through the FASTA files in the current directory
		for fasta_file in os.listdir(dir_path):
			input_file_path = os.path.join(dir_path, fasta_file)
			gene_name = fasta_file.lower().split(".")[0]
			output_file_path = get_canonical_file_path(gene_name)

			# Check if the file with the gene name does not already exist in the output directory
			if not os.path.exists(output_file_path):
				with open(output_file_path, "a") as output_file:
					for record in SeqIO.parse(input_file_path, "fasta"):
						sequence = record.seq
						output_file.write(f">{species}\n{sequence}\n")
			else:
				# Otherwise, check if the species name is not already present, then append the sequence
				with open(output_file_path, "a") as output_file:
					for record in SeqIO.parse(input_file_path, "fasta"):
						sequence = record.seq
						output_file.write(f"\n>{species}\n{sequence}\n")
