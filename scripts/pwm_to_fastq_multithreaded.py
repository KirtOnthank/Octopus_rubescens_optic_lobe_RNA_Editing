import numpy as np
import concurrent.futures
import os

def read_pwm_file(file_path):
    with open(file_path, 'r') as file:
        content = file.readlines()
    
    pwm_data = {}
    header = None
    for line in content:
        line = line.strip()
        if line.startswith('>'):
            header = line[1:]
            pwm_data[header] = []
        elif header:
            pwm_data[header].append(list(map(float, line.split())))
    
    return pwm_data

def pwm_to_sequence(pwm):
    bases = 'ACGT'
    sequence = []
    qualities = []
    for position in zip(*pwm):
        max_prob = max(position)
        max_base = bases[position.index(max_prob)]
        if max_prob == 1.0:
            phred_score = 42  # Assign the highest possible PHRED score
        else:
            phred_score = int(-10 * np.log10(1 - max_prob + 1e-10))
        sequence.append(max_base)
        qualities.append(phred_score)
    return ''.join(sequence), qualities

def process_chunk(chunk):
    chunk_data = {}
    for header, pwm in chunk:
        sequence, qualities = pwm_to_sequence(pwm)
        chunk_data[header] = (sequence, qualities)
    return chunk_data

def write_fastq(file_path, data):
    with open(file_path, 'w') as file:
        for header, (sequence, qualities) in data.items():
            file.write(f"@{header}\n")
            file.write(f"{sequence}\n")
            file.write(f"+\n")
            file.write(''.join(chr(qual + 33) for qual in qualities) + '\n')

def main(pwm_file_path, fastq_file_path, num_threads):
    # Read the input PWM file
    pwm_data = read_pwm_file(pwm_file_path)

    # Split data into chunks for multithreading
    pwm_items = list(pwm_data.items())
    chunk_size = len(pwm_items) // num_threads
    chunks = [pwm_items[i:i + chunk_size] for i in range(0, len(pwm_items), chunk_size)]

    # Process chunks in parallel
    fastq_data = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = [executor.submit(process_chunk, chunk) for chunk in chunks]
        for future in concurrent.futures.as_completed(futures):
            fastq_data.update(future.result())

    # Write to a FASTQ file
    write_fastq(fastq_file_path, fastq_data)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Convert PWM file to FASTQ format with multithreading")
    parser.add_argument("pwm_file", help="Path to the input PWM file")
    parser.add_argument("fastq_file", help="Path to the output FASTQ file")
    parser.add_argument("--threads", type=int, default=os.cpu_count(), help="Number of threads to use")

    args = parser.parse_args()

    main(args.pwm_file, args.fastq_file, args.threads)