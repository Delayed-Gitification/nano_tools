import argparse
import dnaio
import random
from rapidfuzz import fuzz
from os.path import exists
import gzip
import numpy as np
import time


def rev_c(seq):
	"""
	simple function that reverse complements a given sequence
	"""
	tab = str.maketrans("ACTGN", "TGACN")
	# first reverse the sequence
	seq = seq[::-1]
	# and then complement
	seq = seq.translate(tab)
	return seq


def fast_fuzz(s1, s2):
	if (s1 in s2) or (s2 in s1):
		return 100
	else:
		return fuzz.partial_ratio(s1, s2)


def find_barcode(seq, barcodes, min_score, max_ambiguity):
	scores = []
	names = []

	for name, bc in barcodes.items():
		score = fast_fuzz(seq, bc)
		scores.append(score)
		names.append(name)

	# Reject if similar to >1 barcode (ambiguous)
	if len([a for a in scores if a > max_ambiguity]) > 1:
		return "no_match"

	# We know none are ambiguous, so which check if the highest scoring passes the threshold
	elif max(scores) >= min_score:
		return str(names[scores.index(max(scores))])

	else:
		return "no_match"


def initialise_d(forward_primers, reverse_primers):
	d = {}

	if len(forward_primers) > 0 and len(reverse_primers) > 0:
		for i in list(forward_primers.keys()) + ["no_match"]:
			for j in list(reverse_primers.keys()) + ["no_match"]:
				d[str(i) + "_" + str(j)] = ""
	else:
		if len(forward_primers) > 0:
			for i in list(forward_primers.keys()) + ["no_match"]:
				d[str(i)] = ""

		if len(reverse_primers) > 0:
			for i in list(reverse_primers.keys()) + ["no_match"]:
				d[str(i)] = ""

	return d


def write_out_d(d, output):
	for key, value in d.items():

		filename = output + "_" + key + ".fastq.gz"

		if exists(filename):
			mode = "ab"
		else:
			mode = "wb"

		with gzip.open(filename, mode) as file:
			file.write(value.encode())


def read_primers(filename):
	# Read in primers
	forward_primers = {}
	reverse_primers = {}
	with open(filename) as file:
		for i, line in enumerate(file):
			line2 = line.rstrip().split(",")
			this_name = line2[0]
			this_bc = line2[1].upper()
			f_or_r = line2[2].upper()

			assert f_or_r == "F" or f_or_r == "R"

			if f_or_r == "R":
				this_bc = rev_c(this_bc)

			if f_or_r == "F":
				forward_primers[this_name] = this_bc
			else:
				reverse_primers[this_name] = this_bc

	assert len(forward_primers) > 0 or len(reverse_primers) > 0, "No primers found in csv"

	return forward_primers, reverse_primers


def gen_random_sequences(length, number):
	# Generate random integers representing bases (0: A, 1: T, 2: G, 3: C)
	random_bases = np.random.randint(0, 4, size=(number, length))

	# Create a lookup array for mapping integer codes to DNA bases
	bases = np.array(['A', 'T', 'G', 'C'])

	# Map the integer codes to DNA bases using advanced indexing
	random_sequences = bases[random_bases]

	# Convert to strings
	random_sequences = [''.join(a) for a in random_sequences]

	return random_sequences


def get_fdr(barcodes, length, number, min_score, ambiguity):
	"""
	This function calculates the discovery rate of barcodes in random sequences
	"""
	seqs = gen_random_sequences(length, number)

	matches = 0

	for seq in seqs:
		result = find_barcode(seq, barcodes, min_score, ambiguity)
		if result != "no_match":
			matches += 1

	return 100*matches/number


def get_on_target_success_rate(barcodes, length, number, min_score, ambiguity):
	"""
	This function calculates the fraction of true barcode containing sequences that are successfully detected
	"""
	barcode_list = list(barcodes.values())

	bc_l = len(barcode_list[0])
	upstream_l = (length - bc_l)//2
	downstream_l = length - bc_l - upstream_l

	upstream_seqs = gen_random_sequences(upstream_l, number)
	downstream_seqs = gen_random_sequences(downstream_l, number)

	matches = 0

	for i in range(number):
		result = find_barcode(upstream_seqs[i]+random.choice(barcode_list)+downstream_seqs[i],
		                      barcodes, min_score, ambiguity)
		if result != "no_match":
			matches += 1

	return 100*matches/number


def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--fastq", required=True)
	parser.add_argument("-p", "--primers", required=True, help="This should be a headerless, three column csv with "
	                                                           "columns 1=name,2=barcode,3=f or r")
	parser.add_argument("-s", "--min_score", type=float, default=90, help="When matching barcodes")
	parser.add_argument("--max_ambiguity", type=float, default=84, help="If another barcode has this score or higher, "
	                                                                    "then ignore this read because it's ambiguous")
	parser.add_argument("-l", "--length", type=int, default=100, help="length to look at start and end of reads "
	                                                                 "for barcode")
	parser.add_argument("--ignore_rc", action="store_true", default=False, help="by default it looks at the read and "
	                                                                            "its reverse complement. Use this "
	                                                                            "command to only look at the forward "
	                                                                            "read")
	parser.add_argument("--just_check_fdr", action="store_true", default=False, help="Check FDR of primers")
	parser.add_argument("-o", "--output", required=True)
	parser.add_argument("--chunk", default=500, type=int)
	parser.add_argument("--min_length", default=1, type=int)
	parser.add_argument("--fdr_n", default=10_000, type=int, help="Number of random seqs to calculate FDR")
	args = parser.parse_args()

	assert args.min_score >= args.max_ambiguity

	return args


def rationalise_paired_indexing(forward_results, reverse_results, ignore_rc):
	"""
	This function takes forward and, optionally, reverse results, and tries to rationalise them
	"""

	# If we're ignoring the the reverse results then it's simple
	if ignore_rc:
		bc1 = forward_results[0]
		bc2 = forward_results[1]
		return bc1, bc2

	# If we're looking in both directions then we need to think more
	else:
		if forward_results == ["no_match", "no_match"] and reverse_results != ["no_match", "no_match"]:
			# forward results are bad reverse results are better
			bc1 = reverse_results[0]
			bc2 = reverse_results[1]
		elif forward_results != ["no_match", "no_match"] and reverse_results == ["no_match", "no_match"]:
			# Reverse results are bad and forward results are better
			bc1 = forward_results[0]
			bc2 = forward_results[1]
		else:
			# Both are terrible
			bc1 = "no_match"
			bc2 = "no_match"
		return bc1, bc2


def rationalise_single_index(forward_result, reverse_result, ignore_rc):
	"""
	This function rationalises a single index
	"""

	if ignore_rc:  # then it's very simple
		return forward_result

	else:
		if forward_result == "no_match":
			return reverse_result  # may or may not be no_match
		else:
			return forward_result  # definitely a match


def main():
	# Parse args and get primer sequences
	args = get_args()

	forward_primers, reverse_primers = read_primers(args.primers)
	print("Forward primers:")
	print(forward_primers)
	print("\nReverse primers:")
	print(reverse_primers)
	print("")

	if len(forward_primers) == 0 or len(reverse_primers) == 0:
		single_index = True
		print("Using single indexing")
		if len(forward_primers) > 0:
			use_only_forward = True
		else:
			use_only_forward = False
	else:
		print("Using paired indexing")
		use_only_forward = False
		single_index = False

	# Check FDR with current settings
	if single_index:
		if use_only_forward:
			fdr = get_fdr(forward_primers, args.length, args.fdr_n, args.min_score, args.max_ambiguity)
			true_p = get_on_target_success_rate(forward_primers, args.length, args.fdr_n, args.min_score, args.max_ambiguity)
		else:
			fdr = get_fdr(reverse_primers, args.length, args.fdr_n, args.min_score, args.max_ambiguity)
			true_p = get_on_target_success_rate(reverse_primers, args.length, args.fdr_n, args.min_score, args.max_ambiguity)
		print("FDR is " + str(fdr) + "%")
		print("True positive rate is " + str(true_p) + "%")

		if fdr > 0.1:
			print("WARNING - FDR UNACCEPTABLY HIGH!")
	else:
		fdr1 = get_fdr(forward_primers, args.length, args.fdr_n, args.min_score, args.max_ambiguity)
		fdr2 = get_fdr(reverse_primers, args.length, args.fdr_n, args.min_score, args.max_ambiguity)
		true_p1 = get_on_target_success_rate(forward_primers, args.length, args.fdr_n, args.min_score,
		                                    args.max_ambiguity)
		true_p2 = get_on_target_success_rate(reverse_primers, args.length, args.fdr_n, args.min_score,
		                                    args.max_ambiguity)

		print("FDR is " + str(fdr1) + "% and " + str(fdr2) + "%")
		print("True positive rate is " + str(true_p1) + "% and " + str(true_p2) + "%")

		if fdr1 > 0.1 or fdr2 > 0.1:
			print("WARNING - FDR UNACCEPTABLY HIGH!")

	if args.just_check_fdr:
		return 0

	to_write_d = initialise_d(forward_primers, reverse_primers)

	t1 = time.time()

	with dnaio.open(args.fastq) as file:
		since_written = 0

		# Set defaults for bc1 and bc2
		bc1 = "no_match"
		bc2 = "no_match"
		forward_results = [bc1, bc2]
		reverse_results = [bc1, bc2]

		good = 0

		for i, record in enumerate(file):

			if len(record.sequence) < args.min_length:
				continue

			if i % 1000 == 0 and i > 0:
				print(f"Number of reads: {i}", end="\r")

			# Identify the barcodes
			for rc in [True, False]:
				if rc:
					if args.ignore_rc:
						continue

					seq = rev_c(str(record.sequence))
				
				else:
					seq = str(record.sequence)

				s1 = seq[0:args.length]
				s2 = seq[-args.length:]

				if use_only_forward or not single_index:
					bc1 = find_barcode(s1, forward_primers, args.min_score, args.max_ambiguity)
				if not use_only_forward or not single_index:
					bc2 = find_barcode(s2, reverse_primers, args.min_score, args.max_ambiguity)

				if rc:
					reverse_results = [bc1, bc2]
				else:
					forward_results = [bc1, bc2]

			# Rationalise the barcodes
			if single_index:
				if use_only_forward:
					bc = rationalise_single_index(forward_results[0], reverse_results[0], args.ignore_rc)
				else:
					bc = rationalise_single_index(forward_results[1], reverse_results[1], args.ignore_rc)

				if "no_match" not in bc:
					good += 1

			else:
				bc1, bc2 = rationalise_paired_indexing(forward_results, reverse_results, args.ignore_rc)

				if "no_match" not in [bc1, bc2]:
					good += 1

			# Write out the barcodes
			if single_index:
				to_write_d[bc] += "\n".join(
					["@" + record.name, record.sequence, "+", record.qualities]) + "\n"
			else:
				to_write_d[bc1 + "_" + bc2] += "\n".join(["@"+record.name, record.sequence, "+", record.qualities]) + "\n"
			since_written += 1

			if since_written == args.chunk:
				write_out_d(to_write_d, args.output)
				to_write_d = initialise_d(forward_primers, reverse_primers)
				since_written = 0

	write_out_d(to_write_d, args.output)

	t2 = time.time()

	time_diff = round(t2-t1, 1)

	print(f"{i} Reads demultiplexed in {time_diff} seconds.")
	print(f"{round(100*good/i, 2)}% of reads succesfully assigned")


if __name__ == "__main__":
	print("### Demultiplex_nanopore_barcodes.py v0.1 ###\n")
	main()




