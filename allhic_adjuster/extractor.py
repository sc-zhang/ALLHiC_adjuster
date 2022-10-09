from os import path, makedirs, listdir


def extract_ctg_with_tour(in_tour, in_ctg_fasta, out_fasta):
    print("Loading contig fasta")
    ctg_db = {}
    with open(in_ctg_fasta, 'r') as fin:
        for line in fin:
            if line[0] == '>':
                sid = line.strip().split()[0][1:]
                ctg_db[sid] = []
            else:
                ctg_db[sid].append(line.strip())
    for sid in ctg_db:
        ctg_db[sid] = ''.join(ctg_db[sid])

    print("Extracting")
    with open(in_tour, 'r') as fin:
        with open(out_fasta, 'w') as fout:
            for line in fin:
                continue
            data = line.strip().split()

            for tig in data:
                fout.write(">%s\n%s\n" % (tig[:-1], ctg_db[tig[:-1]]))

    print("Finished")


def extract_ctg_with_tours(in_tour_dir, in_ctg_fasta, out_dir):
    if not path.exists(out_dir):
        makedirs(out_dir)
    print("Loading contig fasta")
    ctg_db = {}
    with open(in_ctg_fasta, 'r') as fin:
        for line in fin:
            if line[0] == '>':
                sid = line.strip().split()[0][1:]
                ctg_db[sid] = []
            else:
                ctg_db[sid].append(line.strip())
    for sid in ctg_db:
        ctg_db[sid] = ''.join(ctg_db[sid])
    used_ctg_id = {}

    print("Extracting")
    for fn in listdir(in_tour_dir):
        fn_part = fn.split('.')
        if fn_part[-1] != 'tour':
            continue
        print("\t%s" % fn)
        in_tour = '.'.join(fn_part)
        fn_part[-1] = 'fasta'
        out_fasta = '.'.join(fn_part)
        with open(path.join(in_tour_dir, in_tour), 'r') as fin:
            with open(path.join(out_dir, out_fasta), 'w') as fout:
                for line in fin:
                    continue
                data = line.strip().split()

                for tig in data:
                    used_ctg_id[tig] = 1
                    fout.write(">%s\n%s\n" % (tig[:-1], ctg_db[tig[:-1]]))
    with open(path.join(out_dir, 'unanchored.fasta'), 'w') as fout:
        for tig in sorted(ctg_db):
            if tig not in used_ctg_id:
                fout.write(">%s\n%s\n" % (tig, ctg_db[tig]))

    print("Finished")


def extract_seq_with_list(in_fa, in_list, out_fa):
    fa_db = {}
    with open(in_fa, 'r') as fin:
        for line in fin:
            if line[0] == '>':
                sid = line.strip()[1:]
                fa_db[sid] = []
            else:
                fa_db[sid].append(line.strip())

    for sid in fa_db:
        fa_db[sid] = ''.join(fa_db[sid])

    with open(in_list, 'r') as fin:
        with open(out_fa, 'w') as fout:
            for line in fin:
                if line[0] == '#':
                    continue
                else:
                    tig = line.strip().split()[0]
                    fout.write(">%s\n%s\n" % (tig, fa_db[tig]))


def extract_seq_without_list(in_fa, in_list, out_fa):
    fa_db = {}
    with open(in_fa, 'r') as fin:
        for line in fin:
            if line[0] == '>':
                sid = line.strip()[1:]
                fa_db[sid] = []
            else:
                fa_db[sid].append(line.strip())

    for sid in fa_db:
        fa_db[sid] = ''.join(fa_db[sid])

    with open(in_list, 'r') as fin:
        list_db = {}
        for line in fin:
            if line[0] == '#':
                continue
            list_db[line.strip().split()[0]] = 1

    with open(out_fa, 'w') as fout:
        for tig in fa_db:
            if tig not in list_db:
                fout.write(">%s\n%s\n" % (tig, fa_db[tig]))


def extract_by_tour(args):
    in_file = args.input
    seq_file = args.fasta
    out_file = args.output
    if path.isdir(in_file):
        extract_ctg_with_tours(in_file, seq_file, out_file)
    else:
        extract_ctg_with_tour(in_file, seq_file, out_file)


def extract_by_list(args):
    in_file = args.input
    in_list = args.list
    without = args.without
    out_file = args.output
    if not without:
        extract_seq_with_list(in_file, in_list, out_file)
    else:
        extract_seq_without_list(in_file, in_list, out_file)
