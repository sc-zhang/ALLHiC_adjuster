from os import path, makedirs, listdir
from allhic_adjuster.base.file_reader import read_bed


def agp2tour(args):
    in_agp = args.agp
    out_dir = args.outdir

    agp_db = {}
    print("Reading agp")
    with open(in_agp, 'r') as fin:
        for line in fin:
            data = line.strip().split()
            if len(line.strip()) == 0 or line[0] == '#' or data[4] != 'W':
                continue
            chrn = data[0]
            tig = data[5]
            dir = data[-1]
            if chrn.startswith("tig"):
                continue
            if chrn not in agp_db:
                agp_db[chrn] = []
            agp_db[chrn].append([tig, dir])

    if not path.exists(out_dir):
        makedirs(out_dir)

    print("Writing tours")
    for chrn in agp_db:
        with open(path.join(out_dir, chrn + ".tour"), 'w') as ftour:
            tigs = []
            for tig, dir in agp_db[chrn]:
                tigs.append(tig + dir)

            ftour.write(" ".join(tigs))

    print("Finished")


def tour2txt(args):
    in_tour = args.input
    in_cntREs = args.countREs
    out_txt = args.output

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


def tours2cluster(args):
    in_tour_dir = args.input
    out_cluster = args.output

    clu_db = {}
    for fn in listdir(in_tour_dir):
        if fn.split('.')[-1] != 'tour':
            continue
        chrn = fn.split('.')[0]
        clu_db[chrn] = []
        with open(path.join(in_tour_dir, fn), 'r') as fin:
            for line in fin:
                continue
            for ctg in line.strip().split():
                clu_db[chrn].append(ctg[:-1])
    with open(out_cluster, 'w') as fout:
        fout.write("#Group\tnContigs\tContigs\n")
        for chrn in sorted(clu_db):
            fout.write("%s\t%d\t%s\n" % (chrn, len(clu_db[chrn]), ' '.join(clu_db[chrn])))


def txt2cluster(args):
    in_txt_dir = args.input
    out_cluster = args.output

    clu_db = {}
    for fn in listdir(in_txt_dir):
        if fn.split('.')[-1] != 'txt':
            continue
        chrn = fn.split('.')[0]
        clu_db[chrn] = []
        with open(path.join(in_txt_dir, fn), 'r') as fin:
            for line in fin:
                if line.strip() == "" or line[0] == '#':
                    continue
                clu_db[chrn].append(line.strip().split()[0])

    with open(out_cluster, 'w') as fout:
        fout.write("#Group\tnContigs\tContigs\n")
        for chrn in sorted(clu_db):
            fout.write("%s\t%d\t%s\n" % (chrn, len(clu_db[chrn]), ' '.join(clu_db[chrn])))


def anchors2circos(args):
    qry_bed = args.query
    ref_bed = args.reference
    anc_file = args.anchors
    out_link = args.output

    qdb = read_bed(qry_bed)
    sdb = read_bed(ref_bed)
    link_db = {}
    with open(anc_file, 'r') as fin:
        for line in fin:
            data = line.strip().split()
            if data[0][0] == '#':
                continue
            qg = data[0]
            sg = data[1]
            if qg not in qdb or sg not in sdb:
                continue
            qchr, qsp, qep = qdb[qg]
            schr, ssp, sep = sdb[sg]
            if qchr not in link_db:
                link_db[qchr] = {}
            if schr not in link_db[qchr]:
                link_db[qchr][schr] = []
            link_db[qchr][schr].append([qsp, qep, ssp, sep])
    with open(out_link, 'w') as fout:
        i = 0
        for chrx in sorted(link_db):
            for chry in sorted(link_db[chrx]):
                for xsp, xep, ysp, yep in sorted(link_db[chrx][chry]):
                    fout.write("link%d\t%s\t%d\t%d\n" % (i, chrx, xsp, xep))
                    fout.write("link%d\t%s\t%d\t%d\n" % (i, chry, ysp, yep))
                    i += 1


def ragoo2agp(args):
    in_ord_dir = args.input
    in_ctg = args.fasta
    out_agp = args.output

    chr_len_db = {}
    print("Getting contig length")
    with open(in_ctg, 'r') as fin:
        for line in fin:
            if line[0] == '>':
                id = line.strip().split()[0][1:]
                chr_len_db[id] = 0
            else:
                chr_len_db[id] += len(line.strip())

    print("Reading ordering file and writing agp file")
    used_ctg = {}
    with open(out_agp, 'w') as fout:
        for ord_file in sorted(os.listdir(in_ord_dir)):
            if 'orderings.txt' not in ord_file:
                continue
            grp = ord_file.replace('_orderings.txt', '')
            if 'Chr0' in ord_file:
                continue
            sbase = 1
            with open(os.path.join(in_ord_dir, ord_file), 'r') as fin:
                cnt = 1
                for line in fin:
                    data = line.strip().split()
                    ctg = data[0]
                    used_ctg[ctg] = 1
                    direct = data[1]
                    ebase = sbase + chr_len_db[ctg] - 1
                    fout.write(
                        "%s\t%d\t%d\t%d\tW\t%s\t1\t%d\t%s\n" % (grp, sbase, ebase, cnt, ctg, chr_len_db[ctg], direct))
                    sbase += chr_len_db[ctg]
                    ebase = sbase + 99
                    cnt += 1
                    fout.write("%s\t%d\t%d\t%d\tU\t100\tcontig\tyes\tmap\n" % (grp, sbase, ebase, cnt))
                    sbase += 100
                    cnt += 1
        for ctg in sorted(chr_len_db):
            if ctg in used_ctg:
                continue
            fout.write("%s\t1\t%d\t1\tW\t%s\t1\t%d\t+\n" % (ctg, chr_len_db[ctg], ctg, chr_len_db[ctg]))

    print("Finished")


def ragoo2tour(args):
    in_ord = args.input
    out_tour = args.output
    with open(in_ord, 'r') as fin:
        with open(out_tour, 'w') as fout:
            tour_list = []
            for line in fin:
                data = line.strip().split()
                tour_list.append(data[0] + data[1])
            fout.write("%s" % (' '.join(tour_list)))
