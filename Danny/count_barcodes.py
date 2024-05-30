# Description: Analyze sequencing data for sgRNA library distribution
# Author: Juan M. Martinez-Villalobos / Danny Gallant
# Adapted from: Joung, J., et al. (2016)
# Last modified: 2024-05-30

import csv
from collections import OrderedDict
import argparse
import numpy as np
from Bio import SeqIO

KEY = "TTGGCTTTATATATCTTGTGGAAAGGACGA"  # identifies barcode sequence to determine barcode position
KEY_LENGTH = 30
BARCODE_LENGTH = 24

def count_barcodes(input_file, fastq_file, output_file, guide_g, stats_file):
    """
    Creates a dictionary with barcode counts from fastq_file, writes to output_file.
    fastq_file: forward read fastq file
    output_file: csv file to write barcode dictionary to
    dictionary: barcode sequence as key, barcode count as entry
    """
    num_reads = 0  # total number of reads processed
    perfect_matches = 0  # barcodes with perfect match to library
    non_perfect_matches = 0  # number of reads without a perfect match to any barcode

    # Open library sequences and initiate dictionary of read counts for each barcode
    try:
        with open(input_file, mode='r', newline='') as infile:
            reader = csv.reader(infile)
            next(reader)  # Skip the header row if there is one
            dictionary = {rows[0]: {'count': 0, 'name': rows[1]} for rows in reader}
    except Exception as e:
        print(f'Could not open {input_file}: {e}')
        return

    # Open fastq file
    try:
        with open(fastq_file, "r", newline='') as handle:
            readiter = SeqIO.parse(handle, "fastq")
            for record in readiter:  # Contains the seq and Qscore etc.
                num_reads += 1
                read_sequence = str(record.seq).upper()
                key_index = read_sequence.find(KEY)
                if key_index >= 0:
                    barcode = read_sequence[key_index + KEY_LENGTH : key_index + KEY_LENGTH + BARCODE_LENGTH]
                    if barcode in dictionary:
                        dictionary[barcode]['count'] += 1
                        perfect_matches += 1
                    else:
                        non_perfect_matches += 1
    except Exception as e:
        print(f"Could not find fastq file {fastq_file}: {e}")
        return

    # Create ordered dictionary with barcodes and respective counts and output as a CSV file
    dict_sorted = OrderedDict(sorted(dictionary.items(), key=lambda t: t[0]))
    try:
        with open(output_file, 'w', newline='') as csvfile:
            mywriter = csv.writer(csvfile, delimiter=',')
            mywriter.writerow(['Barcode', 'Count', 'Gene'])
            for barcode, data in dict_sorted.items():
                mywriter.writerow([barcode, data['count'], data['name']])
    except Exception as e:
        print(f'Could not write to {output_file}: {e}')
        return

    # Percentage of barcodes that matched perfectly
    if perfect_matches + non_perfect_matches > 0:
        percent_matched = round(perfect_matches / float(perfect_matches + non_perfect_matches) * 100, 1)
    else:
        percent_matched = 0.0

    # Percentage of undetected barcodes with no read counts
    barcodes_with_reads = sum(1 for data in dictionary.values() if data['count'] > 0)
    barcodes_no_reads = len(dictionary) - barcodes_with_reads
    percent_no_reads = round(barcodes_no_reads / float(len(dictionary)) * 100, 1)

    # Skew ratio of top 10% to bottom 10% of barcode counts
    counts = [data['count'] for data in dictionary.values()]
    top_10 = np.percentile(counts, 90)
    bottom_10 = np.percentile(counts, 10)
    skew_ratio = top_10 / bottom_10 if top_10 != 0 and bottom_10 != 0 else 'Not enough perfect matches to determine skew ratio'

    # Write analysis statistics to statistics file
    try:
        with open(stats_file, 'w', newline='') as infile:
            infile.write(f'Number of perfect barcode matches: {perfect_matches}\n')
            infile.write(f'Number of nonperfect barcode matches: {non_perfect_matches}\n')
            infile.write(f'Number of reads processed: {num_reads}\n')
            infile.write(f'Percentage of barcodes that matched perfectly: {percent_matched}\n')
            infile.write(f'Percentage of undetected barcodes: {percent_no_reads}\n')
            infile.write(f'Skew ratio of top 10% to bottom 10%: {skew_ratio}\n')
    except Exception as e:
        print(f'Could not write to {stats_file}: {e}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze sequencing data for sgRNA library distribution')
    parser.add_argument('-f', '--fastq', type=str, dest='fastq_file', help='fastq file name', default='NGS.fastq')
    parser.add_argument('-o', '--output', type=str, dest='output_file', help='output file name', default='library_count.csv')
    parser.add_argument('-i', '--input', type=str, dest='input_file', help='input file name', default='library_sequences.csv')
    parser.add_argument('-s', '--stats', type=str, dest='stats_file', help='statistics file name', default='statistics.txt')
    parser.add_argument('-no-g', dest='guide_g', help='presence of guanine before spacer', action='store_false')
    parser.set_defaults(guide_g=True)
    args = parser.parse_args()

    count_barcodes(args.input_file, args.fastq_file, args.output_file, args.guide_g, args.stats_file)
