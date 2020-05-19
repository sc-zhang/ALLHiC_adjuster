#!/usr/bin/env python
import sys


def merge_group(tour_files, out_tour):
	tour_list = []
	for tour in tour_files.split(','):
		with open(tour, 'r') as fin:
			for line in fin:
				if line.strip() != '':
					last_line = line.strip()
		tour_list.extend(last_line.split())
	with open(out_tour, 'w') as fout:
		fout.write(' '.join(tour_list))


if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Usage: python "+sys.argv[0]+" <tour_files> <out_tour>")
	else:
		tour_files, out_tour = sys.argv[1:]
		merge_group(tour_files, out_tour)
