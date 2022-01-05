#!/usr/bin/env python
import sys
import os


def convert_tours_to_cluster(in_tour_dir, out_clu):
	clu_db = {}
	for fn in os.listdir(in_tour_dir):
		if fn.split('.')[-1] != 'tour':
			continue
		chrn = fn.split('.')[0]
		clu_db[chrn] = []
		with open(os.path.join(in_tour_dir, fn), 'r') as fin:
			for line in fin:
				continue
			for ctg in line.strip().split():
				clu_db[chrn].append(ctg[:-1])
	with open(out_clu, 'w') as fout:
		fout.write("#Group\tnContigs\tContigs\n")
		for chrn in sorted(clu_db):
			fout.write("%s\t%d\t%s\n"%(chrn, len(clu_db[chrn]), ' '.join(clu_db[chrn])))


if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Usage: python %s <in_tour_dir> <out_cluster>"%sys.argv[0])
	else:
		in_tour_dir, out_clu = sys.argv[1:]
		convert_tours_to_cluster(in_tour_dir, out_clu)

