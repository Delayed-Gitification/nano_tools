import pysam
import argparse
from shared_functions import *
import gzip


def add_to_dict(d, key):
	if key in d.keys():
		d[key] += 1
	else:
		d[key] = 1
	return d


def make_key(rname, position, value, insertion=0):
	return ','.join([str(rname), str(position), str(insertion), str(value)])


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-b", "--bam", help="input bam file")
	parser.add_argument("-o", "--output", help="output file name")
	parser.add_argument("-l", "--min_aligned_length", default=20, type=int,
	                    help="Alignments shorter than this are ignored (Default 20)")
	parser.add_argument("--primary", action="store_true", default=False,
	                    help="Add this option to only consider "
	                         "primary alignments. I.e. ignore reads that are 'not primary alignment' or 'secondary "
	                         "alignment'")
	parser.add_argument("--min_intron_length", type=int, default=50, help="If 'deletions' in the read are this "
	                                                                      "length or longer, it is considered to be "
	                                                                      "an intron and are therefore ignored. "
	                                                                      "Set to very larger number to disable. "
	                                                                      "(Default 50)")

	args = parser.parse_args()

	assert args.min_intron_length > 0, "Can't have negative intron length!"
	print("Ignoring deletions of " + str(args.min_intron_length) + " or more")

	flag_d = make_flag_d()
	d = {}  # a dictionary to store pileup results
	n_reads = 0
	skipped = 0
	with pysam.AlignmentFile(args.bam, 'rb') as file:
		for record in file:
			flag_string = flag_d[record.flag]

			if args.primary:
				if "not primary alignment" in flag_string or "secondary alignment" in flag_string:
					continue

			if record.is_unmapped:
				continue

			positions = record.get_reference_positions(full_length=True)
			seq = record.query_sequence
			rname = record.reference_id

			if seq is None:
				skipped += 1
				continue

			first_pos = min([i for i in positions if i is not None])
			last_pos = max([i for i in positions if i is not None])

			if last_pos - first_pos < args.min_aligned_length:
				continue

			n_reads += 1

			# Iterate though all the aligned positions in the read
			started = False
			insertion_counter = 0
			for seq_pos, p in enumerate(positions):
				# seq_pos = position in read sequence, p = position of reference aligned to

				# Find the first alignment position that isn't "None"
				if p == first_pos:
					started = True
					prev_pos = first_pos

				if not started:
					continue  # true alignment hasn't started yet so continue

				if p is not None:  # therefore genuinely aligned to a reference position
					# It is aligned therefore not a position of a read insertion, therefore reset counter
					insertion_counter = 0

					# Add this position to the pileup as it is aligned
					d = add_to_dict(d, make_key(rname, p, seq[seq_pos]))

					# Check if there is a deletion in the read between this position and the previous position
					if args.min_intron_length > p - prev_pos > 1:  # short deletion in read
						# Add deletion to pileup for every single position in deleted region
						for i in range(prev_pos + 1, p):
							d = add_to_dict(d, make_key(rname, i, "del"))

					prev_pos = p

				else:  # i.e. p == None, therefore it's an insertion in the read that isn't in the reference
					insertion_counter += 1
					this_pos = prev_pos
					d = add_to_dict(d, make_key(rname, this_pos, seq[seq_pos], insertion=insertion_counter))

				if p == last_pos:
					break

	print("Skipped " + str(skipped) + " records")
	if args.output[-3:] != ".gz":
		args.output += ".gz"
	with gzip.open(args.output, 'wb') as file:
		file.write(b"reference_name,position,insertion_number,nt,n,total_reads\n")
		to_write = []
		for key, value in d.items():
			to_write.append(key + "," + str(value) + "," + str(n_reads))

		file.write('\n'.join(to_write).encode())


if __name__ == "__main__":
	print("### perform_enhanced_pileup.py v0.1 ###\n")
	main()
