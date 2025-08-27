import argparse
import pysam
import gzip

def read_qnames(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        return set(line.strip() for line in f if line.strip())

def main():
    ap = argparse.ArgumentParser(description="Filter BAM file by QNAME list")
    ap.add_argument("--input_bam", required=True, help="Input BAM file")
    ap.add_argument("--qname_file", required=True, help="Text file with qnames (one per line, may be gzipped)")
    ap.add_argument("--output_bam", required=True, help="Output BAM file")
    ap.add_argument("--against", action="store_true",
                    help="Filter AGAINST the qnames instead of FOR")
    args = ap.parse_args()

    qnames = read_qnames(args.qname_file)

    with pysam.AlignmentFile(args.input_bam, "rb") as infile, \
         pysam.AlignmentFile(args.output_bam, "wb", template=infile) as outfile:

        for read in infile:
            in_list = read.query_name in qnames
            if (in_list and not args.against) or (not in_list and args.against):
                outfile.write(read)

if __name__ == "__main__":
    main()
