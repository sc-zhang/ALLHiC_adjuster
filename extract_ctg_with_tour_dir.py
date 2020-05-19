#!/usr/bin/env python
import sys
import os


def extract_ctg_with_tour(in_tour_dir, in_ctg_fasta, out_dir):
	if not os.path.exists(out_dir):
		os.mkdir(out_dir)
	print("Loading contig fasta")
	ctg_db = {}
	with open(in_ctg_fasta, 'r') as fin:
		for line in fin:
			if line[0] == '>':
				id = line.strip().split()[0][1:]
				ctg_db[id] = []
			else:
				ctg_db[id].append(line.strip())
	for id in ctg_db:
		ctg_db[id] = ''.join(ctg_db[id])
	used_ctg_id = {}
	print("Exracting")
	for fn in os.listdir(in_tour_dir):
		fn_part = fn.split('.')
		if fn_part[-1] != 'tour':
			continue
		print("\t%s"%fn)
		in_tour = '.'.join(fn_part)
		fn_part[-1] = 'fasta'
		out_fasta = '.'.join(fn_part)
		with open(os.path.join(in_tour_dir, in_tour), 'r') as fin:
			with open(os.path.join(out_dir, out_fasta), 'w') as fout:
				for line in fin:
					continue
				data = line.strip().split()

				for tig in data:
					used_ctg_id[tig] = 1
					fout.write(">%s\n%s\n"%(tig[:-1], ctg_db[tig[:-1]]))
	with open(os.path.join(out_dir, 'unanchored.fasta'), 'w') as fout:
		for tig in sorted(ctg_db):
			if tig not in used_ctg_id:
				fout.write(">%s\n%s\n"%(tig, ctg_db[tig]))

	print("Finished")


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Usage: python "+sys.argv[0]+" <in_tour_dir> <in_ctg_fasta> <out_dir>")
	else:
		in_tour_dir, in_ctg_fasta, out_dir = sys.argv[1:]
		extract_ctg_with_tour(in_tour_dir, in_ctg_fasta, out_dir)
