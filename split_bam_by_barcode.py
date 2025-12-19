import pysam
import argparse
import sys
import os
from rapidfuzz import fuzz
from collections import defaultdict
import multiprocessing
from functools import partial

# --- Helper Functions ---

def rev_c(seq):
    """
    Simple function that reverse complements a given sequence
    """
    tab = str.maketrans("ACTGN", "TGACN")
    seq = seq[::-1]
    seq = seq.translate(tab)
    return seq

def fast_fuzz(s1, s2):
    # s1 is the barcode, s2 is the read sequence chunk
    # If the barcode is perfectly inside the chunk, return 100
    if (s1 in s2) or (s2 in s1):
        return 100
    else:
        return fuzz.partial_ratio(s1, s2)

def find_barcode_aggregated(candidates, barcodes):
    """
    Scans ALL provided sequence chunks (candidates) for every barcode.
    Returns dict: { "barcode": max_score }
    """
    barcode_max_scores = {}
    for bc_seq in barcodes:
        current_max = 0
        for chunk in candidates:
            score = fast_fuzz(bc_seq, chunk)
            if score > current_max:
                current_max = score

        barcode_max_scores[bc_seq] = current_max
    return barcode_max_scores

def merge_max_value(dict1, dict2):
    # Start with a copy of the first dictionary
    result = dict1.copy()
    for key, value in dict2.items():
        if key in result:
            result[key] = max(result[key], value)
        else:
            result[key] = value
    return result

def read_barcodes_list(file_path):
    barcodes = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if "barcode" in line.lower() and "A" not in line and "T" not in line:
                continue
            barcodes.append(line.upper())
    return barcodes

# --- Multiprocessing Logic ---

# Global variable to hold barcodes in worker processes 
# (avoids pickling/copying the list for every single task)
global_barcodes = []

def init_worker(barcodes):
    """Initializer to set global barcodes in worker processes"""
    global global_barcodes
    global_barcodes = barcodes

def process_single_sequence(args_tuple):
    """
    Worker function to process a single sequence string.
    Returns the matched barcode or "no_match".
    """
    seq, n5, n3, min_score, max_ambiguity = args_tuple
    
    # 1. Helper to extract regions (inline for worker isolation)
    def get_search_zones(sequence):
        zones = []
        if n5 == 0 and n3 == 0:
            zones.append(sequence)
        else:
            if n5 > 0:
                zones.append(sequence[:n5])
            if n3 > 0:
                zones.append(sequence[-n3:])
        return zones

    # 2. Search Forward
    forward_zones = get_search_zones(seq)
    match_seq_f = find_barcode_aggregated(forward_zones, global_barcodes)

    # 3. Search Reverse
    seq_rc = rev_c(seq)
    rc_zones = get_search_zones(seq_rc)
    match_seq_r = find_barcode_aggregated(rc_zones, global_barcodes)

    # 4. Merge
    combined_match_scores = merge_max_value(match_seq_f, match_seq_r)

    # 5. Determine Winner
    match_seq = "no_match"
    if combined_match_scores:
        scores_above_ambiguity = [v for _, v in combined_match_scores.items() if v >= max_ambiguity]
        
        # Only proceed if not ambiguous
        if len(scores_above_ambiguity) <= 1:
            max_score = max(combined_match_scores.values())
            if max_score >= min_score:
                keys = [k for k, v in combined_match_scores.items() if v == max_score]
                if len(keys) == 1:
                    match_seq = keys[0]
                    
    return match_seq

# --- Main ---

def main():
    parser = argparse.ArgumentParser(description="Split BAM file by barcodes using multicore processing.")
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file")
    parser.add_argument("-c", "--barcodes", required=True, help="Text file with one barcode per line")
    parser.add_argument("-o", "--output_prefix", required=True, help="Prefix for output BAM files")
    parser.add_argument("--min_score", type=float, default=100, help="Minimum matching score")
    parser.add_argument("--max_ambiguity", type=float, default=95, help="Ambiguity threshold")
    parser.add_argument("--n5", type=int, default=0, help="5' search window size")
    parser.add_argument("--n3", type=int, default=0, help="3' search window size")
    parser.add_argument("--threads", type=int, default=os.cpu_count(), help="Number of threads to use")
    
    args = parser.parse_args()

    # Load Barcodes
    barcodes = read_barcodes_list(args.barcodes)
    print(f"Loaded {len(barcodes)} unique barcodes.")
    print(f"Using {args.threads} CPU cores.")

    # Open Input BAM
    try:
        in_bam = pysam.AlignmentFile(args.bam, "rb")
    except ValueError:
        print("Error: Could not open BAM file.")
        sys.exit(1)

    read_buckets = defaultdict(list)
    processed_count = 0
    matched_count = 0

    print("Starting processing (Buffering reads to memory)...")

    # Define batch size (tuning this helps performance)
    BATCH_SIZE = 10000
    
    # Initialize Pool
    pool = multiprocessing.Pool(processes=args.threads, initializer=init_worker, initargs=(barcodes,))

    batch_reads = []
    batch_seqs_args = []

    for read in in_bam.fetch(until_eof=True):
        if read.is_secondary or read.is_supplementary:
            continue

        seq = read.query_sequence
        if not seq:
            continue
        
        # Accumulate batch
        batch_reads.append(read)
        # Prepare arguments for the worker: (seq, n5, n3, min, max)
        batch_seqs_args.append((seq, args.n5, args.n3, args.min_score, args.max_ambiguity))

        # Process batch when full
        if len(batch_reads) >= BATCH_SIZE:
            # Map workers to sequences
            results = pool.map(process_single_sequence, batch_seqs_args)
            
            # Store results
            for r, match_seq in zip(batch_reads, results):
                processed_count += 1
                if match_seq != "no_match":
                    matched_count += 1
                    read_buckets[match_seq].append(r)
            
            print(f"Scanned {processed_count} reads...", end='\r')
            
            # Reset batch
            batch_reads = []
            batch_seqs_args = []

    # Process remaining reads in the final batch
    if batch_reads:
        results = pool.map(process_single_sequence, batch_seqs_args)
        for r, match_seq in zip(batch_reads, results):
            processed_count += 1
            if match_seq != "no_match":
                matched_count += 1
                read_buckets[match_seq].append(r)

    print(f"\nScanning complete. Writing {matched_count} reads to {len(read_buckets)} output files...")

    pool.close()
    pool.join()

    # Write buckets to disk
    file_counter = 0
    for barcode, reads in read_buckets.items():
        file_counter += 1
        out_filename = f"{args.output_prefix}_{barcode}.bam"
        try:
            with pysam.AlignmentFile(out_filename, "wb", template=in_bam) as out_f:
                for r in reads:
                    out_f.write(r)
        except Exception as e:
            print(f"Error writing file for {barcode}: {e}")

    in_bam.close()
    print(f"Done. Processed {processed_count} reads. Assigned {matched_count} reads.")

if __name__ == "__main__":
    main()