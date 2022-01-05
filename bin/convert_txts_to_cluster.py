#!/usr/bin/env python
import sys
import os


def convert_txts_to_cluster(in_txt_dir, out_clu):
	clu_db = {}
	for fn in os.listdir(in_txt_dir):
		if fn.split('.')[-1] != 'txt':
			continue
		chrn = fn.split('.')[0]
		clu_db[chrn] = []
		with open(os.path.join(in_txt_dir, fn), 'r') as fin:
			for line in fin:
				if line.strip() == "" or line[0] == '#':
					continue
				clu_db[chrn].append(line.strip().split()[0])
	
	with open(out_clu, 'w') as fout:
		fout.write("#Group\tnContigs\tContigs\n")
		for chrn in sorted(clu_db):
			fout.write("%s\t%d\t%s\n"%(chrn, len(clu_db[chrn]), ' '.join(clu_db[chrn])))


if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Usage: python %s <in_txt_dir> <out_cluster>"%sys.argv[0])
	else:
		in_txt_dir, out_clu = sys.argv[1:]
		convert_txts_to_cluster(in_txt_dir, out_clu)

