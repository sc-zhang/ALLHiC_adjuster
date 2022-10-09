from shutil import copy


def merge_tours(args):
    tour_files = args.input
    out_tour = args.output
    tour_list = []
    for tour in tour_files.split(','):
        with open(tour, 'r') as fin:
            for line in fin:
                if line.strip() != '':
                    last_line = line.strip()
        tour_list.extend(last_line.split())
    with open(out_tour, 'w') as fout:
        fout.write(' '.join(tour_list))


def reverse_tour(args):
    tour_file = args.input
    tour_lines = []
    with open(tour_file, 'r') as fin:
        for line in fin:
            if line.strip() != '':
                tour_lines.append(line.strip())
    last_data = tour_lines[-1].split()
    new_last = []
    for i in range(len(last_data) - 1, -1, -1):
        ctg = last_data[i][:-1]
        direction = '+' if last_data[i][-1] == '-' else '-'
        new_last.append(ctg + direction)
    tour_lines[-1] = ' '.join(new_last)
    copy(tour_file, tour_file + '.bak')
    with open(tour_file, 'w') as fout:
        fout.write('\n'.join(tour_lines))


def split_group(grp_name, brk_ctgs):
    with open(grp_name + '.tour', 'r') as fin:
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
    with open(grp_name + '.txt', 'r') as fin:
        for line in fin:
            if line[0] == "#":
                header = line
            else:
                data = line.strip().split()
                ctg_db[data[0]] = line
    for i in range(0, len(ctg_grp)):
        with open(grp_name + "_" + str(i + 1) + ".txt", 'w') as fout:
            fout.write(header)
            for ctg in ctg_grp[i]:
                fout.write(ctg_db[ctg])


def split_tour(tour_file, brk_ctgs):
    with open(tour_file, 'r') as fin:
        for line in fin:
            if line.strip() != '':
                last_line = line.strip()
    data = last_line.split()
    brk_ctgs = brk_ctgs.split(',')
    ctg_grp = [[]]
    for ctg in data:
        if ctg[:-1] in brk_ctgs:
            ctg_grp.append([])
        ctg_grp[-1].append(ctg)

    for i in range(0, len(ctg_grp)):
        with open(tour_file[:-5] + "_" + str(i + 1) + ".tour", 'w') as fout:
            fout.write(' '.join(ctg_grp[i]))


def split_file(args):
    in_file = args.input
    brk_ctgs = args.contigs
    if in_file.endswith('.tour'):
        split_tour(in_file, brk_ctgs)
    else:
        split_group(in_file, brk_ctgs)
