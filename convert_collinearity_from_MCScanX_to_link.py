#!/usr/bin/python
import sys


def get_col(in_col, in_gff, out_txt):
	id_db = {}
	with open(in_gff, 'r') as f_gff:
		for line in f_gff:
			data = line.strip().split()
			id_db[data[1]] = [data[0], data[2]]

	with open(in_col, 'r') as f_col:
		with open(out_txt, 'w') as f_out:
			for line in f_col:
				if line[0] == "#":
					continue
				data = line.strip().split('\t')
				id_1 = data[1]
				id_2 = data[2]
				f_out.write("%s\t%s\n"%('\t'.join(id_db[id_1]), '\t'.join(id_db[id_2])))
			

if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Usage: python " + sys.argv[0] + " <collinearity_file> <gff_file> <out_link>")
	else:
		in_col = sys.argv[1]
		in_gff = sys.argv[2]
		out_txt = sys.argv[3]
		get_col(in_col, in_gff, out_txt)

