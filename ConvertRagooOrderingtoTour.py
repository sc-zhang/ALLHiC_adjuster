#!/usr/bin/env python
import sys


def convert_ord_to_tour(in_ord, out_tour):
	with open(in_ord, 'r') as fin:
		with open(out_tour, 'w') as fout:
			tour_list = []
			for line in fin:
				data = line.strip().split()
				tour_list.append(data[0]+data[1])
			fout.write("%s"%(' '.join(tour_list)))


if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Usage: python "+sys.argv[0]+" <in_ordering> <out_tour>")
	else:
		in_ord, out_tour = sys.argv[1:]
		convert_ord_to_tour(in_ord, out_tour)

