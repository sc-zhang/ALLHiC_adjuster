#!/usr/bin/env python
import sys
import math


def Euc_dist(a, b):
	return math.sqrt((b[0]-a[0])**2+(b[1]-a[1])**2)


def get_break_blocks(in_bed, dtr, out_blocks):
	bed_db = {}
	chr_len_db = {}
	with open(in_bed, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			qchr = data[0]
			schr = data[2]
			qp = int(float(data[1]))
			sp = int(float(data[-1]))
			if qchr not in bed_db:
				bed_db[qchr] = {}
			if schr not in bed_db[qchr]:
				bed_db[qchr][schr] = []
			bed_db[qchr][schr].append([qp, sp])
			if qchr not in chr_len_db:
				chr_len_db[qchr] = {}
			if schr not in chr_len_db[qchr]:
				chr_len_db[qchr][schr] = [0, 0]
			if qp > chr_len_db[qchr][schr][0]:
				chr_len_db[qchr][schr][0] = qp
			if sp > chr_len_db[qchr][schr][1]:
				chr_len_db[qchr][schr][1] = sp

	block_db = {}
	for qchr in bed_db:
		block_db[qchr] = {}
		for schr in bed_db[qchr]:
			groups = []
			chr_len_merge = Euc_dist([0, 0], chr_len_db[qchr][schr])
			tmp_list = []
			for x, y in bed_db[qchr][schr]:
				if len(groups) == 0:
					groups.append([[x, y]])
				else:
					is_add = False
					for i in range(0, len(groups)):
						tail_x, tail_y = groups[i][-1]
						head_x, head_y = groups[i][0]
						if Euc_dist([x, y], [tail_x, tail_y])*dtr/chr_len_merge < 1:
							groups[i].append([x, y])
							is_add = True
							break
					if is_add == False:
						groups.append([[x, y]])

			block_db[qchr][schr] = []								
			for group in groups:
				x = []
				y = []
				for i in range(0, len(group)):
					x.append(group[i][0])
					y.append(group[i][1])
				min_y = min(y)
				max_y = max(y)
				min_index = y.index(min_y)
				max_index = y.index(max_y)
				min_x = x[min_index]
				max_x = x[max_index]
				sx, sy = group[0]
				ex, ey = group[-1]
				if min_x > max_x:
					tmp = min_x
					min_x = max_x
					max_x = tmp
					tmp = min_y
					min_y = max_y
					max_y = tmp
				tmp_list = []
				tmp_list.extend([sx, sy])
				if min_x > sx and min_x < ex:
					tmp_list.extend([min_x, min_y])
				if max_x > sx and max_x < ex:
					tmp_list.extend([max_x, max_y])
				tmp_list.extend([ex, ey])
				for i in range(0, len(tmp_list)-2, 2):
					block_db[qchr][schr].append([tmp_list[i], tmp_list[i+1], tmp_list[i+2], tmp_list[i+3]])
				
	with open(out_blocks, 'w') as fout:
		for qchr in sorted(block_db):
			for schr in sorted(block_db[qchr]):
				for x1, y1, x2, y2 in block_db[qchr][schr]:
					fout.write("%s\t%s\t%d\t%d\t%d\t%d\n"%(qchr, schr, x1, y1, x2, y2))


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Usage: python "+sys.argv[0]+" <in_link> <dist_threshold> <out_blocks>")
		print("\t<dist_threshold> means 1/dist_threshold of chromosome\'s length")
	else:
		in_bed, dtr, out_blocks = sys.argv[1:]
		dtr = float(dtr)*1.0
		get_break_blocks(in_bed, dtr, out_blocks)