#!/usr/bin/env python
import sys
import os


def convert_ord_to_agp(in_ord_dir, in_ctg, out_agp):
	chr_len_db = {}
	print("Getting contig length")
	with open(in_ctg, 'r') as fin:
		for line in fin:
			if line[0] == '>':
				id = line.strip()[1:]
				chr_len_db[id] = 0
			else:
				chr_len_db[id] += len(line.strip())
	
	print("Reading ordering file and writing agp file")
	i = 1
	with open(out_agp, 'w') as fout:
		for ord_file in os.listdir(in_ord_dir):
			if 'orderings.txt' not in ord_file:
				continue
			grp = "group%d"%i
			sbase = 1
			with open(os.path.join(in_ord_dir, ord_file), 'r') as fin:
				cnt = 1
				for line in fin:
					data = line.strip().split()
					ctg = data[0]
					direct = data[1]
					ebase = sbase+chr_len_db[ctg]-1
					fout.write("%s\t%d\t%d\t%d\tW\t%s\t1\t%d\t%s\n"%(grp, sbase, ebase, cnt, ctg, chr_len_db[ctg], direct))
					sbase += chr_len_db[ctg]
					ebase = sbase+99
					cnt += 1
					fout.write("%s\t%d\t%d\t%d\tU\t100\tcontig\tyes\tmap\n"%(grp, sbase, ebase, cnt))
					sbase += 100
					cnt += 1
	
	print("Finished")					


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Usage: python "+sys.argv[0]+" <in_ordering_dir> <in_ctg_fasta> <out_agp>")
	else:
		in_ord_dir, in_ctg, out_agp = sys.argv[1:]
		convert_ord_to_agp(in_ord_dir, in_ctg, out_agp)

