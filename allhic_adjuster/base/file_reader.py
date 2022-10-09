def read_bed(bed_file):
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
