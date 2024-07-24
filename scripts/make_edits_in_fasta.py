import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import sys
import os

def apply_edits(edits_df, sequences_dict, edit_type):
    for index, row in edits_df.iterrows():
        orf = row['orf']
        pos = int(row['pos'])
        mrna_con = row['mrna_con'].upper()
        upstream_base = row['upstream_base'].upper()
        downstream_base = row['downstream_base'].upper()
        if edit_type == 'unedited':
            new_base = row['gdna_con'].upper()
        elif edit_type == 'edited':
            new_base = row['edited'].upper()
        else:
            raise ValueError("Invalid edit type. Choose 'unedited' or 'edited'.")

        # Find the matching sequence in the FASTA file
        matching_headers = [header for header in sequences_dict if header.startswith(orf)]
        if not matching_headers:
            continue

        header = matching_headers[0]
        sequence_record = sequences_dict[header]
        sequence = list(sequence_record.seq.upper())

        # Check for the correct position (+1, 0, -1)
        for offset in [-1, 0, 1]:
            actual_pos = pos + offset
            if (sequence[actual_pos] == mrna_con and
                sequence[actual_pos - 1] == upstream_base and
                sequence[actual_pos + 1] == downstream_base):
                sequence[actual_pos] = new_base
                break

        # Update the sequence in the dictionary
        sequence_record.seq = Seq(''.join(sequence))
        sequences_dict[header] = sequence_record

    return list(sequences_dict.values())

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 make_edits_in_fasta.py <path_to_csv> <path_to_fasta>")
        sys.exit(1)

    csv_file_path = sys.argv[1]
    fasta_file_path = sys.argv[2]

    # Load the CSV file
    edits_df = pd.read_csv(csv_file_path)

    # Load the FASTA file
    sequences = list(SeqIO.parse(fasta_file_path, "fasta"))

    # Create a dictionary for fast access to sequences by header key
    sequences_dict = {record.id: record for record in sequences}

    # Copy the original sequences to ensure we start with the original unmodified sequences for both edits
    unedited_sequences_dict = {record.id: record[:] for record in sequences}
    edited_sequences_dict = {record.id: record[:] for record in sequences}

    # Apply edits for "unedited" file
    unedited_sequences = apply_edits(edits_df, unedited_sequences_dict, 'unedited')

    # Apply edits for "edited" file
    edited_sequences = apply_edits(edits_df, edited_sequences_dict, 'edited')

    # Generate output file paths
    base_name = os.path.splitext(fasta_file_path)[0]
    unedited_output_path = f"{base_name}_unedited.fasta"
    edited_output_path = f"{base_name}_edited.fasta"

    # Save the results to new FASTA files
    SeqIO.write(unedited_sequences, unedited_output_path, "fasta")
    SeqIO.write(edited_sequences, edited_output_path, "fasta")

    print(f"Unedited file saved to: {unedited_output_path}")
    print(f"Edited file saved to: {edited_output_path}")

if __name__ == "__main__":
    main()
