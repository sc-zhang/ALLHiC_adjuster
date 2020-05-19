#!/usr/bin/env python
import sys 
import shutil


def reverse_tour(tour_file):
    tour_lines = []
    with open(tour_file, 'r') as fin:
        for line in fin:
            if line.strip() != '': 
                tour_lines.append(line.strip())
    last_data = tour_lines[-1].split()
    new_last = []
    for i in range(len(last_data)-1, -1, -1):
        ctg = last_data[i][:-1]
        dir = '+' if last_data[i][-1] == '-' else '-' 
        new_last.append(ctg+dir)
    tour_lines[-1] = ' '.join(new_last)
    shutil.copy(tour_file, tour_file+'.bak')
    with open(tour_file, 'w') as fout:
        fout.write('\n'.join(tour_lines))


if __name__ == "__main__":
	if len(sys.argv) < 2:
		print("Usage: python "+sys.argv[0]+" <in_tour>")
	else:
		in_tour = sys.argv[1]
		reverse_tour(in_tour)
