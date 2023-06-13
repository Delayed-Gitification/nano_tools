import argparse
import pysam


def write_out(file, to_write):
	file.write(to_write)


def get_flags(flag):
	"""
	Creates a human-interpretable string for a given flag tag
	"""
	# Convert to binary representation
	bin_string = bin(flag)[2:]

	components = []
	for i, one_or_zero in bin_string[::-1]:
		if one_or_zero == "1":
			components.append(2**int(i))

	flag_list = []
	if 1 in components:
		flag_list.append("read paired")

	if 2 in components:
		flag_list.append("read mapped in proper pair")

	if 4 in components:
		flag_list.append("read unmapped")

	if 8 in components:
		flag_list.append("mate unmapped")

	if 16 in components:
		flag_list.append("read reverse strand")

	if 32 in components:
		flag_list.append("mate reverse strand")

	if 64 in components:
		flag_list.append("first in pair")

	if 128 in components:
		flag_list.append("second in pair")

	if 256 in components:
		flag_list.append("not primary alignment")

	if 512 in components:
		flag_list.append("read fails platform or vendor quality checks")

	if 1024 in components:
		flag_list.append("read is PCR or optical duplicate")

	if 2048 in components:
		flag_list.append("supplementary alignment")

	return ";".join(flag_list)


def make_flag_d():
	"""
	This function makes a dictionary so that we can rapidly get a string with info on the flags
	"""
	return {i: get_flags(i) for i in range(2*2048)}


def analyse_record(record, flag_d, min_intron_length):
	"""
	This function analyses a single bam record
	"""

	flag_string = flag_d(record.flag)

	if "not primary alignment" in flag_string:
		secondary = True
	else:
		secondary = False

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

	to_write = ','.join(
		[ref, str(mapping_quality), flag_string, strand, str(first_pos), str(last_pos)])

	to_write += "," + ";".join(junctions)

	return to_write


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-b", "--bam", required=True)
	parser.add_argument("--min_intron_length", default=50, type=int)
	parser.add_argument("-o", "--output", required=True)
	parser.add_argument("--early_stop", default=-1, type=int)
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

		with open(args.output, 'w') as out:
			out.write("reference,mapping_quality,flag_string,strand,first_pos,last_pos,junctions,number_of_reads\n")
			to_write = []

			for key, value in output_d.items():
				to_write.append(key + "," + str(value))

			out.write('\n'.join(to_write))

	print(str(skipped) + " skipped of " + str(record_number))


if __name__ == '__main__':
	print("### extract_splice_junctions_from_bam.py v0.1 ###\n")
	main()
