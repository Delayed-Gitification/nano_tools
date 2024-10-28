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


def reverse_complement(dna_sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(dna_sequence))


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
    parser.add_argument('--search_flanks', action='store_true', default=False, help='Also use a simpler approach of'
                                                                                    'searching for the barcode flanking sequences. May '
                                                                                    'moderately increase barcode numbers for records '
                                                                                    'where the alignment failed.')
    parser.add_argument('--flank_size', default=8, type=int, help="When using --search_flanks, this sets how long "
                                                                  "the flanks are.")
    args = parser.parse_args()

    # First, read in the fasta file and search for barcode positions for each
    fasta_dict = read_fasta_to_dict(args.fasta)

    barcode_position_d = {}
    upstream_flanks = set()
    downstream_flanks = set()
    bc_with_flanks_found = False  # default

    for name, sequence in fasta_dict.items():
        assert args.bc_seq in sequence, "Can't find barcode in fasta reference sequence: " + name
        bc_start = sequence.find(args.bc_seq)
        bc_end = bc_start + len(args.bc_seq)
        barcode_position_d[name] = {'start': bc_start, 'end': bc_end}

        if args.search_flanks:
            upstream_flanks.add(sequence[bc_start - args.flank_size:bc_start])
            downstream_flanks.add(sequence[bc_end:bc_end + args.flank_size])

    if args.search_flanks:
        assert len(upstream_flanks) == 1 and len(downstream_flanks) == 1, (
            'Cannot use --search_flanks because flanks are '
            'inconsistent in different reference sequences')
        upstream_flank = upstream_flanks.pop()
        downstream_flank = downstream_flanks.pop()

    # Now go through each record in the bam file
    to_write = ["query_name,barcode,min_phred,mean_phred"]
    skipped_due_to_phred = 0
    skipped_due_to_wrong_size = 0
    skipped_due_to_barcode_not_found = 0

    with pysam.AlignmentFile(args.bam, "rb") as samfile:
        for record in samfile:
            if record.is_secondary:
                continue
            if record.is_supplementary:
                continue

            if args.search_flanks:
                bc_with_flanks_found = False  # assume not
                barcode_start = record.query_sequence.find(upstream_flank) + len(upstream_flank)
                barcode_end = record.query_sequence.find(downstream_flank)
                if barcode_end - barcode_start == len(args.bc_seq):
                    barcode = record.query_sequence[barcode_start:barcode_end]
                    phred = record.query_qualities[barcode_start:barcode_end]
                    bc_with_flanks_found = True

                if not bc_with_flanks_found:

                    sequence_rev_c = reverse_complement(record.query_sequence)

                    barcode_start = sequence_rev_c.find(upstream_flank) + len(upstream_flank)
                    barcode_end = sequence_rev_c.find(downstream_flank)
                    if barcode_end - barcode_start == len(args.bc_seq):
                        barcode = reverse_complement(sequence_rev_c[barcode_start:barcode_end])
                        phred = record.query_qualities[::-1][barcode_start:barcode_end]
                        bc_with_flanks_found = True

            if not bc_with_flanks_found or not args.search_flanks:
                if record.mapping_quality < args.min_mapq:
                    continue
                if record.is_unmapped:
                    continue

                assert record.reference_name in fasta_dict.keys(), 'Reference fasta does not contain sequence: ' + record.reference_name

                pair_d = make_pair_d(record.get_aligned_pairs())

                try:
                    barcode_start = pair_d[barcode_position_d[record.reference_name]['start']]
                    barcode_end = pair_d[barcode_position_d[record.reference_name]['end']]
                except:
                    continue

                if barcode_start is None or barcode_end is None:
                    skipped_due_to_barcode_not_found += 1
                    continue

                # Check barcode is the right length
                if barcode_end - barcode_start != len(args.bc_seq):
                    skipped_due_to_wrong_size += 1
                    continue

                barcode = record.query_sequence[barcode_start:barcode_end]
                phred = record.query_qualities[barcode_start:barcode_end]

            min_phred = min(phred)
            if min_phred < args.min_phred:
                skipped_due_to_phred += 1
                continue

            mean_phred = statistics.mean(phred)
            if mean_phred < args.mean_phred:
                skipped_due_to_phred += 1
                continue

            to_write.append(f"{record.query_name},{barcode},{min_phred},{mean_phred}")

    with gzip.open(args.output, 'wb') as out:
        out.write(('\n'.join(to_write) + '\n').encode())

    print(f"{skipped_due_to_barcode_not_found} reads skipped due to barcode not found")
    print(f"{skipped_due_to_wrong_size} reads skipped due to wrong barcode length")
    print(f"{skipped_due_to_phred} reads skipped due to failing phred-score threshold")
    print(f"{len(to_write) - 1} barcodes extracted successfully")


if __name__ == "__main__":
    main()