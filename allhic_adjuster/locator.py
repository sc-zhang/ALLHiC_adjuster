from allhic_adjuster.base.file_reader import read_bed
from allhic_adjuster.base.file_reader import read_agp
from math import sqrt
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use('Agg')


def convert_anchors(qry_bed, ref_bed, anc_file):
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
            qchr, qp = qdb[qg]
            schr, sp = sdb[sg]
            if qchr not in link_db:
                link_db[qchr] = {}
            if schr not in link_db[qchr]:
                link_db[qchr][schr] = []
            link_db[qchr][schr].append([qp, sp])

    link_list = []
    for chrx in sorted(link_db):
        for chry in sorted(link_db[chrx]):
            for x, y in sorted(link_db[chrx][chry]):
                link_list.append([chrx, x, chry, y])

    return link_list


def euc_dist(a, b):
    return sqrt((b[0] - a[0]) ** 2 + (b[1] - a[1]) ** 2)


def get_break_blocks(link_list, resolution):
    link_db = {}
    chr_len_db = {}
    for qchr, qp, schr, sp in link_list:
        if qchr not in link_db:
            link_db[qchr] = {}
        if schr not in link_db[qchr]:
            link_db[qchr][schr] = []
        link_db[qchr][schr].append([qp, sp])
        if qchr not in chr_len_db:
            chr_len_db[qchr] = {}
        if schr not in chr_len_db[qchr]:
            chr_len_db[qchr][schr] = [0, 0]
        if qp > chr_len_db[qchr][schr][0]:
            chr_len_db[qchr][schr][0] = qp
        if sp > chr_len_db[qchr][schr][1]:
            chr_len_db[qchr][schr][1] = sp

    block_db = {}
    for qchr in link_db:
        block_db[qchr] = {}
        for schr in link_db[qchr]:
            groups = []
            chr_len_merge = euc_dist([0, 0], chr_len_db[qchr][schr])
            for x, y in link_db[qchr][schr]:
                if len(groups) == 0:
                    groups.append([[x, y]])
                else:
                    is_add = False
                    for i in range(0, len(groups)):
                        tail_x, tail_y = groups[i][-1]
                        if euc_dist([x, y], [tail_x, tail_y]) * resolution / chr_len_merge < 1:
                            groups[i].append([x, y])
                            is_add = True
                            break
                    if is_add == False:
                        groups.append([[x, y]])

            block_db[qchr][schr] = []
            for group in groups:
                x = []
                y = []
                for i in range(0, len(group)):
                    x.append(group[i][0])
                    y.append(group[i][1])
                min_y = min(y)
                max_y = max(y)
                min_index = y.index(min_y)
                max_index = y.index(max_y)
                min_x = x[min_index]
                max_x = x[max_index]
                sx, sy = group[0]
                ex, ey = group[-1]
                if min_x > max_x:
                    tmp = min_x
                    min_x = max_x
                    max_x = tmp
                    tmp = min_y
                    min_y = max_y
                    max_y = tmp
                tmp_list = []
                tmp_list.extend([sx, sy])
                if sx < min_x < ex:
                    tmp_list.extend([min_x, min_y])
                if sx < max_x < ex:
                    tmp_list.extend([max_x, max_y])
                tmp_list.extend([ex, ey])
                for i in range(0, len(tmp_list) - 2, 2):
                    block_db[qchr][schr].append([tmp_list[i], tmp_list[i + 1], tmp_list[i + 2], tmp_list[i + 3]])

    return block_db


def convert_link(link_list):
    data_db = {}
    chr_list_x = []
    chr_list_y = []
    for chr_x, pos_x, chr_y, pos_y in link_list:
        if chr_x not in chr_list_x:
            chr_list_x.append(chr_x)
        if chr_y not in chr_list_y:
            chr_list_y.append(chr_y)
        if chr_x not in data_db:
            data_db[chr_x] = {}
        if chr_y not in data_db[chr_x]:
            data_db[chr_x][chr_y] = []
        data_db[chr_x][chr_y].append([pos_x, pos_y])

    chr_list_x = sorted(chr_list_x)
    chr_list_y = sorted(chr_list_y)
    return chr_list_x, chr_list_y, data_db


def get_ctg_pos(region, pos):
    s = 0
    e = len(region) - 1
    while s <= e:
        mid = int((s + e) / 2)
        if region[mid][0] > pos:
            e = mid - 1
        elif region[mid][0] < pos:
            s = mid + 1
        else:
            return mid
    if region[e][1] >= pos:
        return e
    elif e == len(region) - 1:
        return e
    else:
        return -1


def draw_dot_plot(link_list, block_db, in_agp, resolution, out_pic):
    print("Reading data")
    chr_list_x, chr_list_y, data_db = convert_link(link_list)
    agp_db = read_agp(in_agp)
    chr_len_db = {}

    print("Calculating chromosomes length")
    for chrx in chr_list_x:
        if chrx not in chr_len_db:
            chr_len_db[chrx] = 0
        for chry in chr_list_y:
            if chry not in chr_len_db:
                chr_len_db[chry] = 0
            if chry not in data_db[chrx]:
                continue
            for x, y in data_db[chrx][chry]:
                if x > chr_len_db[chrx]:
                    chr_len_db[chrx] = x
                if y > chr_len_db[chry]:
                    chr_len_db[chry] = y
    base_x = 0
    base_y = 0
    offset_db = {}
    for chrx in chr_list_x:
        offset_db[chrx] = base_x
        base_x += chr_len_db[chrx]
    for chry in chr_list_y:
        offset_db[chry] = base_y
        base_y += chr_len_db[chry]

    print("Converting data")
    data_x = []
    data_y = []
    for chrx in chr_list_x:
        for chry in chr_list_y:
            if chry not in data_db[chrx]:
                continue
            for x, y in data_db[chrx][chry]:
                data_x.append(x + offset_db[chrx])
                data_y.append(y + offset_db[chry])

    block_x = []
    block_y = []
    label_x = []

    print("Writing contig blocks")
    out_block = out_pic.split('.')
    out_block[-1] = 'block.txt'
    out_block = '.'.join(out_block)
    with open(out_block, 'w') as fout:
        ctg_block_db = {}
        for chrx in chr_list_x:
            for chry in chr_list_y:
                if chrx not in block_db or chry not in block_db[chrx]:
                    continue
                for x1, x2, y1, y2 in block_db[chrx][chry]:
                    if euc_dist([x1, y1], [x2, y2]) * resolution / euc_dist([0, 0],
                                                                            [chr_len_db[chrx], chr_len_db[chry]]) < 1:
                        continue
                    block_x.append([x1 + offset_db[chrx], x2 + offset_db[chrx]])
                    rstart = get_ctg_pos(agp_db[chrx], x1)
                    rend = get_ctg_pos(agp_db[chrx], x2)
                    if rstart == -1 or rend == -1:
                        print(rstart, rend)
                    ctg1 = agp_db[chrx][rstart][2]
                    ctg2 = agp_db[chrx][rend][2]
                    ctg_list = []
                    for i in range(rstart, rend + 1):
                        ctg_list.append(agp_db[chrx][i][2] + agp_db[chrx][i][3])
                    label_x.append([ctg1, ctg2])
                    if chrx not in ctg_block_db:
                        ctg_block_db[chrx] = []
                    ctg_block_db[chrx].append([x1, x2, ctg_list])
                    block_y.append([y1 + offset_db[chry], y2 + offset_db[chry]])
        for chrx in chr_list_x:
            if chrx not in ctg_block_db:
                continue
            for x1, x2, ctg_list in sorted(ctg_block_db[chrx]):
                fout.write(">BLOCK_%s_%d_%d\n%s\n" % (chrx, x1, x2, ' '.join(ctg_list)))

    max_x = 0
    max_y = 0
    for chrx in chr_list_x:
        max_x += chr_len_db[chrx]

    for chry in chr_list_y:
        max_y += chr_len_db[chry]

    print("Plotting")
    plt.figure(figsize=(10, 10), dpi=300)
    x_ticks = []
    x_labels = []
    base_x = 0
    for chrx in chr_list_x:
        plt.plot([chr_len_db[chrx] + base_x, chr_len_db[chrx] + base_x], [0, max_y], linestyle='-', color='green',
                 linewidth=0.5, markersize=0)
        x_ticks.append(base_x + int(chr_len_db[chrx] / 2))
        x_labels.append(chrx)
        base_x += chr_len_db[chrx]

    y_ticks = []
    y_labels = []
    base_y = 0
    for chry in chr_list_y:
        plt.plot([0, max_x], [chr_len_db[chry] + base_y, chr_len_db[chry] + base_y], linestyle='-', color='green',
                 linewidth=0.5, markersize=0)
        y_ticks.append(base_y + int(chr_len_db[chry] / 2))
        y_labels.append(chry)
        base_y += chr_len_db[chry]
    plt.plot(data_x, data_y, linestyle='', color='black', marker='o', markersize=0.5)

    for i in range(0, len(block_x)):
        plt.plot(block_x[i], block_y[i], linestyle='-', color='red', linewidth=0.5, markersize=0)
        plt.annotate("%s" % label_x[i][0], xy=(block_x[i][0], block_y[i][0]), fontsize=2.5, color='red')
        plt.annotate("%s" % label_x[i][1], xy=(block_x[i][1], block_y[i][1]), fontsize=2.5, color='red')
    plt.xlim([0, max_x])
    plt.ylim([0, max_y])
    plt.xticks(x_ticks)
    plt.yticks(y_ticks)
    ax = plt.gca()
    ax.set_xticklabels(x_labels, rotation=45)
    ax.set_yticklabels(y_labels, rotation=0)
    ax.xaxis.set_ticks_position('top')
    ax.yaxis.set_ticks_position('right')
    ax.invert_yaxis()
    ax.tick_params(top=False, right=False)
    print("Saving picture")
    plt.savefig(out_pic, filetype=out_pic.split('.')[-1], bbox_inches='tight')
    print("Finished")


def locator(args):
    qry_bed = args.query
    ref_bed = args.reference
    anc_file = args.anchors
    agp_file = args.agp
    resolution = args.resolution
    out_pic = args.outpic
    link_list = convert_anchors(qry_bed, ref_bed, anc_file)
    block_db = get_break_blocks(link_list, resolution)
    draw_dot_plot(link_list, block_db, agp_file, resolution, out_pic)
