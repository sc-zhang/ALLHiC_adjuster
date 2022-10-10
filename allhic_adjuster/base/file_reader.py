"""
This file contain several functions of file operate
"""


def read_bed(bed_file):
    """
    This function is used for reading bed file, and return a dictionary
    key -> gene id
    value -> a list contain: chromosome name, start position, end position
    """
    bed_db = {}
    with open(bed_file, 'r') as fin:
        for line in fin:
            data = line.strip().split()
            chrn = data[0]
            sp = int(data[1])
            ep = int(data[2])
            if sp > ep:
                sp, ep = ep, sp
            gene = data[3]
            bed_db[gene] = [chrn, sp, ep]
    return bed_db


def read_agp(in_agp):
    """
    This function is used for reading AGP file, and return a dictionary
    key -> chromosome id
    value -> two dimensions list:
            [[start position, end position, contig id, direction]
            ...]
    """
    agp_db = {}
    with open(in_agp, 'r') as fin:
        for line in fin:
            data = line.strip().split()
            if len(line.strip()) == 0 or line[0] == '#' or data[4] == 'U':
                continue
            chr_x = data[0]
            sp = int(float(data[1]))
            ep = int(float(data[2]))
            ctg = data[5]
            direct = data[-1]
            if chr_x not in agp_db:
                agp_db[chr_x] = []
            agp_db[chr_x].append([sp, ep, ctg, direct])
    return agp_db


def read_fasta(in_fasta):
    """
    This function is used for reading fasta file, and return a dictionary
    key -> chromosome id
    value -> sequence
    """
    seq_db = {}
    with open(in_fasta, 'r') as fin:
        for line in fin:
            if line[0] == '>':
                sid = line.strip().split()[0][1:]
                seq_db[sid] = []
            else:
                seq_db[sid].append(line.strip())
    for sid in seq_db:
        seq_db[sid] = ''.join(seq_db[sid])

    return seq_db
