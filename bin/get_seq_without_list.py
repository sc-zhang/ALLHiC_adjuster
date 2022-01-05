#!/usr/bin/env python
import sys


def get_seq_with_list(in_fa, in_list, out_fa):
	fa_db = {}
	with open(in_fa, 'r') as fin:
		for line in fin:
			if line[0] == '>':
				id = line.strip()[1:]
				fa_db[id] = []
			else:
				fa_db[id].append(line.strip())
	
	for id in fa_db:
		fa_db[id] = ''.join(fa_db[id])
	

	with open(in_list, 'r') as fin:
		list_db = {}
		for line in fin:
			if line[0] == '#':
				continue
			list_db[line.strip().split()[0]] = 1

	with open(out_fa, 'w') as fout:
		for tig in fa_db:
			if tig not in list_db:	
				fout.write(">%s\n%s\n"%(tig, fa_db[tig]))


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Usage: python %s <in_fa> <in_list> <out_fa>"%sys.argv[0])
	else:
		in_fa, in_list, out_fa = sys.argv[1:]
		get_seq_with_list(in_fa, in_list, out_fa)

