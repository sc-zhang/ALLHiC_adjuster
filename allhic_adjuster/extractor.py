from os import path, makedirs, listdir
from allhic_adjuster.base.file_reader import read_fasta


def extract_ctg_with_tour(in_tour, in_ctg_fasta, out_fasta):
    """
    This function is used for extracting contig sequences with tour file
    """
    print("Loading contig fasta")
    ctg_db = read_fasta(in_ctg_fasta)

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
    """
    This function is used for extracting contig sequences with tour files
    """
    if not path.exists(out_dir):
        makedirs(out_dir)
    print("Loading contig fasta")
    ctg_db = read_fasta(in_ctg_fasta)
    used_ctg_id = set()

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
                    used_ctg_id.add(tig)
                    fout.write(">%s\n%s\n" % (tig[:-1], ctg_db[tig[:-1]]))

    with open(path.join(out_dir, 'unanchored.fasta'), 'w') as fout:
        for tig in sorted(ctg_db):
            if tig not in used_ctg_id:
                fout.write(">%s\n%s\n" % (tig, ctg_db[tig]))

    print("Finished")


def extract_seq_in_list(in_fa, in_list, out_fa):
    """
    This function is used for extracting contig sequences which id in list file
    """
    print("Loading fasta")
    fa_db = read_fasta(in_fa)

    print("Extracting")
    with open(in_list, 'r') as fin:
        with open(out_fa, 'w') as fout:
            for line in fin:
                if line[0] == '#':
                    continue
                else:
                    tig = line.strip().split()[0]
                    fout.write(">%s\n%s\n" % (tig, fa_db[tig]))

    print("Finished")


def extract_seq_not_in_list(in_fa, in_list, out_fa):
    """
    This function is used for extracting contig sequences which id not in list file
    """
    print("Loading fasta")
    fa_db = read_fasta(in_fa)

    print("Loading list")
    id_set = set()
    with open(in_list, 'r') as fin:
        for line in fin:
            if line[0] == '#':
                continue
            id_set.add(line.strip().split()[0])

    print("Extracting")
    with open(out_fa, 'w') as fout:
        for tig in fa_db:
            if tig not in id_set:
                fout.write(">%s\n%s\n" % (tig, fa_db[tig]))

    print("Finished")


def extract_by_tour(args):
    """
    This function is used for extracting contig sequences with tour file or tour files
    """
    in_file = args.input
    seq_file = args.fasta
    out_file = args.output
    if path.isdir(in_file):
        extract_ctg_with_tours(in_file, seq_file, out_file)
    else:
        extract_ctg_with_tour(in_file, seq_file, out_file)


def extract_by_list(args):
    """
    This function is used for extracting contig sequences which id in or not in list file
    """
    in_file = args.input
    in_list = args.list
    is_not_in = args.not_in
    out_file = args.output
    if not is_not_in:
        extract_seq_in_list(in_file, in_list, out_file)
    else:
        extract_seq_not_in_list(in_file, in_list, out_file)
