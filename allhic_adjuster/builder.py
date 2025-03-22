from allhic_adjuster.base.file_reader import read_fasta
import os


def rev_seq(seq):
    base_db = {"A": "T", "T": "A", "G": "C", "C": "G",
               "a": "t", "t": "a", "g": "c", "c": "g"}
    return ''.join([base_db[_] if _ in base_db else _ for _ in seq[::-1]])


def builder(args):
    contig_fa = args.ref
    tours_dir = args.input
    out_dir = args.output

    print("Loading contigs")
    ctg_db = read_fasta(contig_fa)

    print("Building genome and writing agp file")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out_genome = os.path.join(out_dir, "groups.asm.fasta")
    out_agp = os.path.join(out_dir, "groups.agp")

    f_genome = open(out_genome, "w")
    f_agp = open(out_agp, "w")
    anchored_ctgs = set()
    gap = "N" * 100

    for tour in sorted(os.listdir(tours_dir)):
        if tour.startswith("."):
            continue
        if not tour.endswith(".tour"):
            continue

        last_line = ""
        with open(os.path.join(tours_dir, tour), "r") as fin:
            for line in fin:
                if line:
                    last_line = line
        cur_tour = last_line.split()

        idx = 1
        sp = 1
        cur_chrn = ".".join(tour.split(".")[:-1])
        cur_seqs = []
        print("\tBuilding %s" % cur_chrn)

        for _ in range(len(cur_tour)):
            ctg = cur_tour[_][:-1]
            direct = cur_tour[_][-1]
            ctg_len = len(ctg_db[ctg])
            f_agp.write(
                "%s\t%d\t%d\t%d\tW\t%s\t1\t%d\t%s\n" % (cur_chrn, sp, sp + ctg_len - 1, idx, ctg, ctg_len, direct))
            sp += ctg_len
            idx += 1

            if _ != len(cur_tour) - 1:
                f_agp.write("%s\t%d\t%d\t%d\tU\t100\tcontig\tyes\tmap\n" % (cur_chrn, sp, sp + 99, idx))
                idx += 1
                sp += 100
            cur_seqs.append(ctg_db[ctg] if direct == '+' else rev_seq(ctg_db[ctg]))
            anchored_ctgs.add(ctg)

        f_genome.write(">%s\n%s\n" % (cur_chrn, gap.join(cur_seqs)))

    print("Building retain contigs")

    for ctg in ctg_db:
        if ctg in anchored_ctgs:
            continue
        f_agp.write("%s\t1\t%d\t1\tW\t%s\t1\t%d\t+\n" % (ctg, len(ctg_db[ctg]), ctg, len(ctg_db[ctg])))
        f_genome.write(">%s\n%s\n" % (ctg, ctg_db[ctg]))

    f_genome.close()
    f_agp.close()

    print("Finished")
