#!/usr/bin/env python
import sys


def remove_redundancy_tigs(in_txt, in_allele, group_name, out_txt):
	print("Loading allele")
	allele_list = set()
	with open(in_allele, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			chrn = data[0]
			if len(data[2:]) < 2:
				continue
			if chrn.lower() != group_name.lower():
				continue
			tigs = '\t'.join(data[2:])
			allele_list.add(tigs)
	
	new_allele_list = []
	for tigs in allele_list:
		data = tigs.split('\t')
		new_allele_list.append(data)
	
	tig_cnt = {}
	for tigs in new_allele_list:
		for tig in tigs:
			if tig not in tig_cnt:
				tig_cnt[tig] = 0
			tig_cnt[tig] += 1
	
	print("Loading contigs")
	tig_info = {}
	with open(in_txt, 'r') as fin:
		for line in fin:
			if line[0] == '#':
				continue
			data = line.strip().split()
			tig_info[data[0]] = data[1:]

	remove_list = {}
	save_list = {}
	tig_cnt = {}
	for tigs in new_allele_list:
		tmp_list = []
		for tig in tigs:
			if tig not in tig_cnt:
				tig_cnt[tig] = 0
			tig_cnt[tig] += 1
			if tig in tig_info:
				tig_len = int(tig_info[tig][-1])
			else:
				tig_len = 0
			tmp_list.append([tig_len, tig])
		sorted_tmp_list = sorted(tmp_list, reverse=True)
		save_list[sorted_tmp_list[0][1]] = 1
		for cnt, tig in sorted_tmp_list[1:]:
			if tig not in save_list:
				remove_list[tig] = 1
				
	print("Writing contigs")
	with open(out_txt, 'w') as fout:
		fout.write("#Contig\tRECounts\tLength\n")
		for tig in tig_info:
			if tig not in remove_list:
				fout.write("%s\t%s\n"%(tig, '\t'.join(tig_info[tig])))

	print("Finished")


if __name__ == "__main__":
	if len(sys.argv) < 5:
		print("Usage: python %s <in_txt> <in_allele> <group_name> <out_txt>"%sys.argv[0])
	else:
		in_txt, in_allele, group_name, out_txt = sys.argv[1:]
		remove_redundancy_tigs(in_txt, in_allele, group_name, out_txt)

