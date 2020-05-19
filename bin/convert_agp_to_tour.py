#!/usr/bin/env python
import sys
import os


def convert_agp_to_tour(in_agp, out_dir):
	agp_db = {}
	print("Reading agp")
	with open(in_agp, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			if data[4] != 'W':
				continue
			chrn = data[0]
			tig = data[5]
			dir = data[-1]
			if chrn.startswith("tig"):
				continue
			if chrn not in agp_db:
				agp_db[chrn] = []
			agp_db[chrn].append([tig, dir])
	
	if not os.path.exists(out_dir):
		os.mkdir(out_dir)
	
	print("Writing result")
	for chrn in agp_db:
		with open(os.path.join(out_dir, chrn+".tour"), 'w') as ftour:
			tigs = []
			for tig, dir in agp_db[chrn]:
				tigs.append(tig+dir)
				
			ftour.write(" ".join(tigs))

	print("Finished")


if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Usage: python "+sys.argv[0]+" <in_agp> <out_dir>")
	else:
		in_agp, out_dir = sys.argv[1:]
		convert_agp_to_tour(in_agp, out_dir)

