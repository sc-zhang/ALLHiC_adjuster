#!/usr/bin/env python
import sys


def split_group(grp_name, brk_ctgs):
	with open(grp_name, 'r') as fin:
		for line in fin:
			if line.strip() != '':
				last_line = line.strip()
	data = last_line.split()
	brk_ctgs = brk_ctgs.split(',')
	ctg_grp = [[]]
	for ctg in data:
		if ctg[:-1] in brk_ctgs:
			ctg_grp.append([])
		ctg_grp[-1].append(ctg)
	
	for i in range(0, len(ctg_grp)):
		with open(grp_name[:-5]+"_"+str(i+1)+".tour", 'w') as fout:
			fout.write(' '.join(ctg_grp[i]))


if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Usage: python "+sys.argv[0]+" <group_tour> <break_contigs>")
	else:
		grp_name, brk_ctgs = sys.argv[1:]
		split_group(grp_name, brk_ctgs)

