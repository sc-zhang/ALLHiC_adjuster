#!/usr/bin/env python
import sys


def read_bed(bed_file):
	bed_db = {}
	with open(bed_file, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			chrn = data[0]
			sp = int(data[1])
			ep = int(data[2])
			if sp > ep:
				sp, ep = ep, sp
			gene = data[3]
			bed_db[gene] = [chrn, sp, ep]
	return bed_db


def convert_anchors(query_bed, sub_bed, anc_file, out_link):
	qdb = read_bed(query_bed)
	sdb = read_bed(sub_bed)
	link_db = {}
	with open(anc_file, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			if data[0][0] == '#':
				continue
			qg = data[0]
			sg = data[1]
			if qg not in qdb or sg not in sdb:
				continue
			qchr, qsp, qep = qdb[qg]
			schr, ssp, sep = sdb[sg]
			if qchr not in link_db:
				link_db[qchr] = {}
			if schr not in link_db[qchr]:
				link_db[qchr][schr] = []
			link_db[qchr][schr].append([qsp, qep, ssp, sep])
	with open(out_link, 'w') as fout:
		i = 0
		for chrx in sorted(link_db):
			for chry in sorted(link_db[chrx]):
				for xsp, xep, ysp, yep in sorted(link_db[chrx][chry]):
					fout.write("link%d\t%s\t%d\t%d\n"%(i, chrx, xsp, xep))
					fout.write("link%d\t%s\t%d\t%d\n"%(i, chry, ysp, yep))
					i += 1


if __name__ == "__main__":
	if len(sys.argv) < 5:
		print("Usage: python "+sys.argv[0]+" <query_bed> <sub_bed> <anchors_file> <out_link>")
	else:
		query_bed, sub_bed, anc_file, out_link = sys.argv[1:]
		convert_anchors(query_bed, sub_bed, anc_file, out_link)