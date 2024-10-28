import pysam
import argparse
import statistics
import gzip


def make_pair_d(pairs):
    # This function makes a dictionary to enable quick finding of the read position which corresponds to a given
    # reference position
    d = {}
    for pair in pairs:
        seq_pos = pair[0]
        ref_pos = pair[1]
        if ref_pos is None:
            continue
        d[ref_pos] = seq_pos

    return d


def phred_to_numeric(phred_string):
    return [ord(char) - 33 for char in phred_string]


def read_fasta_to_dict(fasta_file):
    fasta_dict = {}
    with open(fasta_file, 'r') as file:
        identifier = None
        sequence = []
        for line in file:
            line = line.strip()
            if line.startswith(">"):  # Header line
                if identifier:  # Save the previous entry
                    fasta_dict[identifier] = ''.join(sequence)
                identifier = line[1:]  # Remove the ">" and use the rest as identifier
                sequence = []  # Reset the sequence list for the new entry
            else:
                sequence.append(line)  # Collect sequence lines
        if identifier:  # Save the last entry
            fasta_dict[identifier] = ''.join(sequence)
    return fasta_dict


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bam", required=True, help="Bam file")
    parser.add_argument("-o", "--output", required=True, help="Output filename")
    parser.add_argument('-f', '--fasta', required=True, help="Reference fasta used for alignment")
    parser.add_argument('-bc', '--bc_seq', required=True, help="The randomers used to represent the barcode "
                                                               "in the fasta file. Eg NNNNNNN or HHHHH")
    parser.add_argument("--min_mapq", required=False, type=float, help="Minimum mapping quality",
                        default=0)
    parser.add_argument('--min_phred', required=False, default=0, help="Minimum phred score of barcode")
    parser.add_argument('--mean_phred', required=False, default=0, help="Minimum average phred score of barcode")
    args = parser.parse_args()

    # First, read in the fasta file and search for barcode positions for each
    fasta_dict = read_fasta_to_dict(args.fasta)
    barcode_position_d = {}
    for name, sequence in fasta_dict.items():
        assert args.bc_seq in sequence, "Can't find barcode in fasta reference sequence: " + name
        bc_start = sequence.find(args.bc_seq)
        bc_end = bc_start + len(args.bc_seq)
        barcode_position_d[name] = {'start': bc_start, 'end': bc_end}

    # Now go through each record in the bam file
    to_write = ["qname,barcode,min_phred,mean_phred"]
    skipped_due_to_phred = 0
    with pysam.AlignmentFile(args.bam, "rb") as samfile:
        for record in samfile:
            if record.is_secondary:
                continue
            if record.is_supplementary:
                continue
            if record.mapping_quality < args.min_mapq:
                continue
            if record.is_unmapped:
                continue

            assert record.rname in fasta_dict.keys(), 'Reference fasta does not contain sequence: ' + record.rname

            pair_d = make_pair_d(record.get_aligned_pairs())

            try:
                barcode_start = pair_d[fasta_dict[record.rname['start']]]
                barcode_end = pair_d[fasta_dict[record.rname['end']]]
            except:
                continue

            if barcode_start is None or barcode_end is None:
                continue

            # Check barcode is the right length
            if barcode_end - barcode_start != fasta_dict[record.rname['end']] - fasta_dict[record.rname['start']]:
                continue

            barcode = record.query_sequence[barcode_start:barcode_end]
            phred_string = record.query_qualities[barcode_start:barcode_end]
            phred = phred_to_numeric(phred_string)

            min_phred = min(phred)
            if min_phred < args.min_phred:
                skipped_due_to_phred += 1
                continue

            mean_phred = statistics.mean(phred)
            if mean_phred < args.mean_phred:
                skipped_due_to_phred += 1
                continue

            to_write.append(f"{record.qname},{barcode},{min_phred},{mean_phred}")

    with gzip.open(args.output, 'wb') as out:
        out.write(('\n'.join(to_write) + '\n').encode())

    print(f"{skipped_due_to_phred} reads skipped due to failing phred-score threshold")


if __name__ == "__main__":
    main()