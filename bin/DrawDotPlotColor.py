#!/usr/bin/env python
import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import re
import math


def Euc_dist(a, b):
	return math.sqrt((b[0]-a[0])**2+(b[1]-a[1])**2)


def read_agp(in_agp):
	agp_db = {}
	with open(in_agp, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			if data[4] == 'U':
				continue
			#chr_pre_x, index_x = re.findall(r'(\S.*?)(\d.*)', data[0])[0]
			#chr_x = "%s%02d"%(chr_pre_x, int(index_x))
			chr_x = data[0]
			sp = int(float(data[1]))
			ep = int(float(data[2]))
			ctg = data[5]
			direct = data[-1]
			if chr_x not in agp_db:
				agp_db[chr_x] = []
			agp_db[chr_x].append([sp, ep, ctg, direct])
	return agp_db


def read_table(in_link):
	data_db = {}
	chr_list_x =[]
	chr_list_y =[]
	with open(in_link, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			chr_x = data[0]
			chr_y = data[2]
			pos_x = int(float(data[1]))
			pos_y = int(float(data[3]))
			if chr_x not in chr_list_x:
				chr_list_x.append(chr_x)
			if chr_y not in chr_list_y:
				chr_list_y.append(chr_y)
			if chr_x not in data_db:
				data_db[chr_x] = {}
			if chr_y not in data_db[chr_x]:
				data_db[chr_x][chr_y] = []
			data_db[chr_x][chr_y].append([pos_x, pos_y])
	
	chr_list_x = sorted(chr_list_x)
	chr_list_y = sorted(chr_list_y)
	return chr_list_x, chr_list_y, data_db

def read_block(in_block):
	block_db = {}
	with open(in_block, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			qchr = data[0]
			schr = data[1]
			x1 = int(data[2])
			y1 = int(data[3])
			x2 = int(data[4])
			y2 = int(data[5])
			if qchr not in block_db:
				block_db[qchr] = {}
			if schr not in block_db[qchr]:
				block_db[qchr][schr] = []
			block_db[qchr][schr].append([x1, x2, y1, y2])
	return block_db


def get_ctg_pos(region, pos):
	s = 0
	e = len(region)-1
	while s<=e:
		mid = int((s+e)/2)
		if region[mid][0] > pos:
			e = mid-1
		elif region[mid][0] < pos:
			s = mid+1
		else:
			return mid
	if region[e][1] >= pos:
		return e
	elif e == len(region)-1:
		return e
	else:
		return -1


def draw_dot_plot(in_link, in_block, in_agp, block_s, out_pic):
	print("Reading data")
	chr_list_x, chr_list_y, data_db = read_table(in_link)
	block_db = read_block(in_block)
	agp_db = read_agp(in_agp)
	chr_len_db = {}
	
	print("Calculating chromosomes length")
	for chrx in chr_list_x:
		if chrx not in chr_len_db:
			chr_len_db[chrx] = 0
		for chry in chr_list_y:
			if chry not in chr_len_db:
				chr_len_db[chry] = 0
			if chry not in data_db[chrx]:
				continue
			for x,y in data_db[chrx][chry]:
				if x > chr_len_db[chrx]:
					chr_len_db[chrx] = x
				if y > chr_len_db[chry]:
					chr_len_db[chry] = y
	base_x = 0
	base_y = 0
	offset_db = {}
	for chrx in chr_list_x:
		offset_db[chrx] = base_x
		base_x += chr_len_db[chrx]
	for chry in chr_list_y:
		offset_db[chry] = base_y
		base_y += chr_len_db[chry]				

	print("Converting data")
	data_x = []
	data_y = []
	for chrx in chr_list_x:
		for chry in chr_list_y:
			if chry not in data_db[chrx]:
				continue
			for x,y in data_db[chrx][chry]:
				data_x.append(x+offset_db[chrx])
				data_y.append(y+offset_db[chry])
	
	block_x = []
	block_y = []
	label_x = []
	
	print("Writing contig infos")
	with open("ctg_blocks.txt", 'w') as fout:
		ctg_block_db = {}
		for chrx in chr_list_x:
			for chry in chr_list_y:
				if chrx not in block_db or chry not in block_db[chrx]:
					continue
				for x1, x2, y1, y2 in block_db[chrx][chry]:
					if Euc_dist([x1, y1], [x2, y2])*block_s/Euc_dist([0, 0], [chr_len_db[chrx], chr_len_db[chry]]) < 1:
						continue
					block_x.append([x1+offset_db[chrx], x2+offset_db[chrx]])
					rstart = get_ctg_pos(agp_db[chrx], x1)
					rend = get_ctg_pos(agp_db[chrx], x2)
					if rstart == -1 or rend == -1:
						print(rstart, rend)
					ctg1 = agp_db[chrx][rstart][2]
					ctg2 = agp_db[chrx][rend][2]
					ctg_list = []
					for i in range(rstart, rend+1):
						ctg_list.append(agp_db[chrx][i][2]+agp_db[chrx][i][3])
					label_x.append([ctg1, ctg2])
					if chrx not in ctg_block_db:
						ctg_block_db[chrx] = []
					ctg_block_db[chrx].append([x1, x2, ctg_list])
					block_y.append([y1+offset_db[chry], y2+offset_db[chry]])
		for chrx in chr_list_x:
			if chrx not in ctg_block_db:
				continue
			for x1, x2, ctg_list in sorted(ctg_block_db[chrx]):
				fout.write(">BLOCK_%s_%d_%d\n%s\n"%(chrx, x1, x2, ' '.join(ctg_list)))
	
	max_x = 0
	max_y = 0
	for chrx in chr_list_x:
		max_x += chr_len_db[chrx]
	
	for chry in chr_list_y:
		max_y += chr_len_db[chry]
	
	print("Plotting")
	plt.figure(figsize=(10, 10), dpi=300)
	x_ticks = []
	x_labels = []
	base_x = 0
	for chrx in chr_list_x:
		plt.plot([chr_len_db[chrx]+base_x, chr_len_db[chrx]+base_x], [0, max_y], linestyle='-', color='green', linewidth=0.5, markersize=0)
		x_ticks.append(base_x+int(chr_len_db[chrx]/2))
		x_labels.append(chrx)
		base_x += chr_len_db[chrx]
	
	y_ticks = []
	y_labels = []
	base_y = 0
	for chry in chr_list_y:
		plt.plot([0, max_x], [chr_len_db[chry]+base_y, chr_len_db[chry]+base_y], linestyle='-', color='green', linewidth=0.5, markersize=0)
		y_ticks.append(base_y+int(chr_len_db[chry]/2))
		y_labels.append(chry)
		base_y += chr_len_db[chry]
	plt.plot(data_x, data_y, linestyle ='', color='black', marker='o', markersize=0.5)
	
	for i in range(0, len(block_x)):
		plt.plot(block_x[i], block_y[i], linestyle='-', color='red', linewidth=0.5, markersize=0)
		plt.annotate("%s"%label_x[i][0], xy=(block_x[i][0], block_y[i][0]), fontsize=2.5, color='red')
		plt.annotate("%s"%label_x[i][1], xy=(block_x[i][1], block_y[i][1]), fontsize=2.5, color='red')
	plt.xlim([0, max_x])
	plt.ylim([0, max_y])
	plt.xticks(x_ticks)
	plt.yticks(y_ticks)
	ax = plt.gca()
	ax.set_xticklabels(x_labels, rotation=45)
	ax.set_yticklabels(y_labels, rotation=0)
	ax.xaxis.set_ticks_position('top')
	ax.yaxis.set_ticks_position('right')
	ax.invert_yaxis()
	ax.tick_params(top=False, right=False)
	print("Saving picture")
	plt.savefig(out_pic, filetype=out_pic.split('.')[-1], bbox_inches='tight')
	print("Finished")


if __name__ == "__main__":
	if len(sys.argv) < 6:
		print("Usage: python "+sys.argv[0]+" <in_link> <in_block> <in_agp> <block_size> <out_pic>")
		print("\t<block_size> means 1/block_size of chromosome\' length")
	else:
		in_link, in_block, in_agp, block_s, out_pic = sys.argv[1:]
		block_s = float(block_s)*1.0
		draw_dot_plot(in_link, in_block, in_agp, block_s, out_pic)
		
