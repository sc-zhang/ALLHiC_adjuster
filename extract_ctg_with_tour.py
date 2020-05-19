#!/usr/bin/env python
import sys


def extract_ctg_with_tour(in_tour, in_ctg_fasta, out_fasta):
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

	print("Exracting")
	with open(in_tour, 'r') as fin:
		with open(out_fasta, 'w') as fout:
			for line in fin:
				continue
			data = line.strip().split()

			for tig in data:
				fout.write(">%s\n%s\n"%(tig[:-1], ctg_db[tig[:-1]]))

	print("Finished")


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Usage: python "+sys.argv[0]+" <in_tour> <in_ctg_fasta> <out_fasta>")
	else:
		in_tour, in_ctg_fasta, out_fasta = sys.argv[1:]
		extract_ctg_with_tour(in_tour, in_ctg_fasta, out_fasta)
