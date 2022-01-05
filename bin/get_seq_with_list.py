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
		with open(out_fa, 'w') as fout:
			for line in fin:
				if line[0] == '#':
					continue
				else:
					tig = line.strip().split()[0]
					fout.write(">%s\n%s\n"%(tig, fa_db[tig]))


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Usage: python %s <in_fa> <in_list> <out_fa>"%sys.argv[0])
	else:
		in_fa, in_list, out_fa = sys.argv[1:]
		get_seq_with_list(in_fa, in_list, out_fa)

