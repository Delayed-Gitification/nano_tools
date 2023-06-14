import argparse
import pysam
import gzip
from shared_functions import *


def analyse_record(record, flag_d, min_intron_length):
	"""
	This function analyses a single bam record
	"""

	flag_string = flag_d[record.flag]

	positions = record.get_reference_positions()

	first_pos = min(positions)
	last_pos = max(positions)

	mapping_quality = record.mapping_quality

	ref = record.reference_name
	junctions = []

	if record.is_reverse:
		strand = "-"
	else:
		strand = "+"

	for i in range(len(positions) - 1):
		distance = positions[i + 1] - positions[i]
		if distance >= min_intron_length:
			junctions.append(str(positions[i]) + "-" + str(positions[i + 1]))

	analysed_string = ','.join(
		[ref, str(mapping_quality), str(record.flag), flag_string, strand, str(first_pos), str(last_pos)])

	analysed_string += "," + ";".join(junctions)

	return analysed_string


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-b", "--bam", required=True, help="Bam file")
	parser.add_argument("-o", "--output", required=True, help="Output filename")
	parser.add_argument("--min_intron_length", default=50, type=int, help="Minimum intron length (default 50)")
	parser.add_argument("--early_stop", default=-1, type=int, help="Stop after N records")
	args = parser.parse_args()

	output_d = {}
	skipped = 0
	record_number = 0

	flag_d = make_flag_d()

	with pysam.AlignmentFile(args.bam) as bam:
		for record in bam:
			record_number += 1
			try:
				if record.is_unmapped:
					to_write = ','.join(["NA", "NA", "NA", "NA"])
					to_write += "," + ";".join([])

				else:
					if record_number % 10_000 == 0:
						print(record_number)

					analysed_string = analyse_record(record, flag_d, args.min_intron_length)

				if analysed_string in output_d.keys():
					output_d[analysed_string] += 1
				else:
					output_d[analysed_string] = 1

				if record_number > args.early_stop > 0:
					break
			except:
				skipped += 1

	if args.output[-3:] != ".gz":
		args.output += ".gz"

	with gzip.open(args.output, 'wb') as out:
		out.write(b"reference,mapping_quality,flag,flag_string,strand,first_pos,last_pos,junctions,number_of_reads\n")
		to_write = []

		for key, value in output_d.items():
			to_write.append(key + "," + str(value))

		out.write('\n'.join(to_write).encode())

	print(str(skipped) + " skipped of " + str(record_number))


if __name__ == '__main__':
	print("### extract_splice_junctions_from_bam.py v0.1 ###\n")
	main()
