#!/usr/bin/env python
import sys


def convert_tour_to_txt(in_tour, in_cntREs, out_txt):
	print("Loading countREs")
	cnt_RE_db = {}
	with open(in_cntREs, 'r') as fin:
		for line in fin:
			if line[0] == '#':
				header = line
			else:
				data = line.strip().split()
				cnt_RE_db[data[0]] = line

	print("Converting")
	with open(in_tour, 'r') as fin:
		with open(out_txt, 'w') as fout:
			for line in fin:
				continue
			data = line.strip().split()

			fout.write(header)
			for tig in data:
				fout.write(cnt_RE_db[tig[:-1]])

	print("Finished")


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Usage: python "+sys.argv[0]+" <in_tour> <in_countREs> <out_txt>")
	else:
		in_tour, in_cntREs, out_txt = sys.argv[1:]
		convert_tour_to_txt(in_tour, in_cntREs, out_txt)
