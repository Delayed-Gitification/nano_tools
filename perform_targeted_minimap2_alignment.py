import mappy as mp
import pysam
import dnaio
import argparse
import pandas as pd

def fasta_to_dict(fasta_path):
    """
    Reads a FASTA file into a dictionary without external libraries.

    Parameters:
        fasta_path (str): Path to the FASTA file.

    Returns:
        dict: {sequence_id: sequence_string}
    """
    seq_dict = {}
    with open(fasta_path, 'r') as f:
        seq_id = None
        seq_chunks = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                # Save previous record
                if seq_id is not None:
                    seq_dict[seq_id] = ''.join(seq_chunks)
                seq_id = line[1:].split()[0]  # take first word as ID
                seq_chunks = []
            else:
                seq_chunks.append(line)
        # Save last record
        if seq_id is not None:
            seq_dict[seq_id] = ''.join(seq_chunks)
    return seq_dict


def csv_to_dict(csv_path, key_col, value_col):
    """
    Reads a CSV file and returns a dictionary with keys and values
    from the specified columns.

    Parameters:
        csv_path (str): Path to the CSV file.
        key_col (str): Name of the column to use as dictionary keys.
        value_col (str): Name of the column to use as dictionary values.

    Returns:
        dict: Dictionary mapping key_col -> value_col
    """
    df = pd.read_csv(csv_path)

    if key_col == 'qname':
        if 'query_name' in df.columns:
            key_col = 'query_name'

    if value_col == 'qname':
        if 'query_name' in df.columns:
            value_col = 'query_name'

    if key_col not in df.columns or value_col not in df.columns:
        raise ValueError(f"CSV must contain columns: {key_col}, {value_col}")

    return dict(zip(df[key_col], df[value_col]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Align FASTQ to FASTA and output BAM, using barcode-target mapping.")

    parser.add_argument('--fasta', required=True, help='Path to reference FASTA file')
    parser.add_argument('--fastq', required=True, help='Path to input FASTQ file (i.e. of the RNA-seq file)')
    parser.add_argument('--output', required=True, help='Path to output BAM file')
    parser.add_argument('--plasmid_barcode_alignment_csv', required=True,
                        help='CSV containing two columns: barcode and reference. These should match the reference names in the supplied fasta file, obviously...')
    parser.add_argument('--RNA_barcode_csv', required=True,
                        help='CSV containing two columns: barcode and qname. This contains the barcode for each read in the RNA seq fastq')
    parser.add_argument("--minimap2_preset", default="splice")

    args = parser.parse_args()

    plasmid_csv_d = csv_to_dict(args.plasmid_barcode_alignment_csv, 'barcode', 'reference')

    RNA_csv_d = csv_to_dict(args.RNA_barcode_csv, 'qname', 'barcode')

    fasta_d = fasta_to_dict(args.fasta)

    # make a dictionary of aligners
    aligner_d = {}
    for key, value in fasta_d.items():
        aligner_d[key] = mp.Aligner(seq=value.upper(), preset=args.minimap2_preset)

    # Prepare a bam file
    header = {
        'HD': {'VN': '1.6'},
        'SQ': [{'LN': len(value), 'SN': key} for key, value in fasta_d.items()]
    }

    bamfile = pysam.AlignmentFile(args.output, "wb", header=header)

    total_records = 0
    total_written = 0

    with dnaio.open(args.fastq) as fastq:
        for record in fastq:
            total_records += 1

            record_name = record.name.split(" ")[0]

            if record_name not in RNA_csv_d.keys():
                continue

            barcode = RNA_csv_d[record_name]

            if barcode not in plasmid_csv_d.keys():
                continue

            target_name = plasmid_csv_d[barcode]

            if target_name not in fasta_d.keys():
                continue

            aligner = aligner_d[target_name]

            hits = aligner.map(str(record.sequence))

            hit = next(hits, None)  # Get the first (best) hit

            if hit is not None:
                aln = pysam.AlignedSegment()
                aln.query_name = record_name

                if hit.strand == -1:
                    aln.flag = 16
                    aln.query_sequence = mp.revcomp(str(record.sequence)[hit.q_st:hit.q_en])
                    aln.query_qualities = pysam.qualitystring_to_array(str(record.qualities)[hit.q_st:hit.q_en])[::-1]
                else:
                    aln.flag = 0
                    aln.query_sequence = str(record.sequence)[hit.q_st:hit.q_en]
                    aln.query_qualities = pysam.qualitystring_to_array(str(record.qualities)[hit.q_st:hit.q_en])

                aln.reference_id = bamfile.get_tid(target_name)
                aln.reference_start = hit.r_st
                aln.mapping_quality = hit.mapq
                aln.query_name = target_name
                aln.cigarstring = hit.cigar_str
                aln.next_reference_id = -1
                aln.next_reference_start = -1
                aln.template_length = 0

                bamfile.write(aln)
                total_written += 1

    bamfile.close()
    print(f"{total_written} records written to bam from a total of {total_records}")

