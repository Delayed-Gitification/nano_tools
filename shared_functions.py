def get_flags(flag):
	"""
	Creates a human-interpretable string for a given sam flag
	"""
	if flag == 0:
		return ""
	# Convert to binary representation
	bin_string = bin(flag)[2:]

	components = []
	for i, one_or_zero in enumerate(bin_string[::-1]):
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

	return "; ".join(flag_list)


def make_flag_d():
	"""
	This function makes a dictionary so that we can rapidly get a string with info on the flags
	"""
	return {i: get_flags(i) for i in range(2*2048)}