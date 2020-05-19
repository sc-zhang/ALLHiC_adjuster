#!/usr/bin/env python
import sys


def split_group(grp_name, brk_ctgs):
	with open(grp_name+'.tour', 'r') as fin:
		for line in fin:
			if line.strip() != '':
				last_line = line.strip()
	data = last_line.split()
	brk_ctgs = brk_ctgs.split(',')
	ctg_grp = [[]]
	for ctg in data:
		if ctg[:-1] in brk_ctgs:
			ctg_grp.append([])
		ctg_grp[-1].append(ctg[:-1])
	ctg_db = {}
	with open(grp_name+'.txt', 'r') as fin:
		for line in fin:
			if line[0] == "#":
				header = line
			else:
				data = line.strip().split()
				ctg_db[data[0]] = line
	for i in range(0, len(ctg_grp)):
		with open(grp_name+"_"+str(i+1)+".txt", 'w') as fout:
			fout.write(header)
			for ctg in ctg_grp[i]:
				fout.write(ctg_db[ctg])


if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Usage: python "+sys.argv[0]+" <group_name> <break_contigs>")
	else:
		grp_name, brk_ctgs = sys.argv[1:]
		split_group(grp_name, brk_ctgs)

