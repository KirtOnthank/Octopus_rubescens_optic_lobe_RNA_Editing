from Bio import SeqIO

# Define file paths
file1_path = 'AtoG_aa_edited_zf_targets_N.fastq'
file2_path = 'AtoG_aa_unedited_zf_targets_N.fastq'
output_path = 'non_identical_headers.txt'

# Read the sequences from both files
sequences1 = {record.id: record.seq for record in SeqIO.parse(file1_path, "fastq")}
sequences2 = {record.id: record.seq for record in SeqIO.parse(file2_path, "fastq")}

# Find headers of non-identical sequences
non_identical_headers = []
for header in sequences1:
    if header in sequences2 and sequences1[header] != sequences2[header]:
        non_identical_headers.append(header)

# Write the headers to the output file
with open(output_path, 'w') as output_file:
    for header in non_identical_headers:
        output_file.write(header + '\n')

print(f"Non-identical headers have been written to {output_path}")
