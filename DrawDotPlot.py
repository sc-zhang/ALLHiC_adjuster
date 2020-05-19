#!/usr/bin/env python
import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import re


def read_table(in_link):
	data_db = {}
	chr_list_x =[]
	chr_list_y =[]
	with open(in_link, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			chr_x = data[0]
			chr_y = data[2]
			chr_pre_x, index_x = re.findall(r'(\S.*?)(\d.*)', chr_x)[0]
			chr_x = chr_pre_x+"%02d"%(int(index_x))
			chr_pre_y, index_y = re.findall(r'(\S.*?)(\d.*)', chr_y)[0]
			if len(index_x) > 2 or len(index_y) > 2:
				continue
			chr_y = chr_pre_y+"%02d"%(int(index_y))
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


def draw_dot_plot(in_link, out_pic):
	print("Reading data")
	chr_list_x, chr_list_y, data_db = read_table(in_link)
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
	if len(sys.argv) < 3:
		print("Usage: python "+sys.argv[0]+" <in_link> <out_pic>")
	else:
		in_link, out_pic = sys.argv[1:]
		draw_dot_plot(in_link, out_pic)
		
