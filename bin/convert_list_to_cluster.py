#!/usr/bin/env python
import sys


def convert_clu(in_list, out_clu):
	with open(in_list, 'r') as fin:
		with open(out_clu, 'w') as fout:
			fout.write("#Group\tnContigs\tContigs\n")
			for line in fin:
				data = line.strip().split()
				sub_data = data[1].split(',')
				fout.write("%s\t%d\t%s\n"%(data[0], len(sub_data), ' '.join(sub_data)))


if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Usage: python %s <in_list> <out_cluster>"%sys.argv[0])
	else:
		in_list, out_clu = sys.argv[1:]
		convert_clu(in_list, out_clu)

