import os
import gzip
import csv

# Define columns to keep
columns_needed = [
    "chr", "pos", "ref", "alt", "rs_dbSNP151",
    "ANNOVAR_refseq_Gene_ID",
    "ANNOVAR_refseq_Closest_gene(intergenic_only)",
    "SnpEff_refseq_Gene_ID",
    "VEP_refseq_Gene_Name",
    "ANNOVAR_ensembl_Gene_ID",
    "ANNOVAR_ensembl_Closest_gene(intergenic_only)",
    "SnpEff_ensembl_Gene_ID",
    "VEP_ensembl_Gene_ID"
]

# Paths
input_folder = r"C:\Users\bryan\OneDrive - University of Southern California\Research\Huaiyu Mi\AnnoQ\AnnoQ_data_all"
output_folder = r"C:\Users\bryan\Desktop\AnnoQ_short"

# Make output folder if needed
os.makedirs(output_folder, exist_ok=True)

# Process each .gz file
for file in os.listdir(input_folder):
    if file.endswith(".gz"):
        input_path = os.path.join(input_folder, file)
        output_path = os.path.join(output_folder, file.replace(".gz", "_short.gz"))

        print(f"🔄 Working on: {file}")

        with gzip.open(input_path, 'rt', encoding='utf-8') as fin, \
             gzip.open(output_path, 'wt', newline='', encoding='utf-8') as fout:

            reader = csv.DictReader(fin, delimiter='\t')
            writer = csv.DictWriter(fout, fieldnames=columns_needed, delimiter='\t')
            writer.writeheader()

            for row in reader:
                writer.writerow({col: row.get(col, "") for col in columns_needed})

        print(f"✅ Saved: {os.path.basename(output_path)}\n")
