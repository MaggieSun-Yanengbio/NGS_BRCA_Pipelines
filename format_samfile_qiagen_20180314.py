import sys
import time

# get overlap of R1_fastq and R2_fastq
def r1r2_overlap(r1_fastq,r2_fastq,collective):
    dic_f1 = {}
    with open(r1_fastq,'r') as f1:
        while True:
                l1 = f1.readline()
                if not l1:
                    break
                l2 = f1.readline()
                l3 = f1.readline()
                l4 = f1.readline()
                dic_f1[l1] = [l1,l2,l3,l4]
    f1.close()
    dic_f2 = {}
    with open(r2_fastq,'r') as f2:
        while True:
                l1 = f2.readline()
                if not l1:
                    break
                l2 = f2.readline()
                l3 = f2.readline()
                l4 = f2.readline()
                dic_f2[l1] = [l1,l2,l3,l4]
    f2.close()
    dic_overlap = {}
    for key1 in dic_f1:
        if key1 in dic_f2:
            dic_overlap[key1] = dic_f1[key1]
    out = open(collective,'w')
    for values in dic_overlap.values():
        if len(values) != 4:
            continue
        else:
            for line in values:
                out.write(line)
    out.close()

# make dictionary of barcode clustered
def make_dic_barcode(all_barcode_file,overlap_file):
    all_barcode = open(all_barcode_file,'r')
    dic_barcode_all = {}
    dic_qname = {}
    for line1 in all_barcode.readlines():
        key = 'N'
        m = 0
        for bar in line1.split(','):
            barcode_seq = bar.split(' ')[-1].split(';')[0]
            qname = bar.split(' ')[0][1:]
            dic_qname[qname] = barcode_seq
            if m == 0:
                key = barcode_seq
            dic_barcode_all[barcode_seq] = key
            m += 1
    print ('total barcode is %s'%len(dic_barcode_all))
    all_barcode.close()

    overlap = open(overlap_file,'r')
    dic_barcode_overlap = {}
    n = 0
    for line2 in overlap.readlines():
        if n%4 == 0:
            barcode_cluster = line2.split(';')[1]
            dic_barcode_overlap[barcode_cluster] = '0'
        n += 1
    print ('filtered barcode is %s'%len(dic_barcode_overlap))
    overlap.close()
    return dic_qname,dic_barcode_all,dic_barcode_overlap

def format_sam(filter_file,dic_qname,dic_barcode_all,dic_barcode_overlap,output):
    number_success = 0
    out = open(output, 'w')
    with open(filter_file,'r') as samfile:
        for sam in samfile:
            name, tag, chromo, pos, a1, a2, a3, pnext, a4, a5, a6 = sam.split('\t')[:11]
            rest = sam.split('\t')[11:]
            chromo, start, end = chromo.split('_')
            pos = int(start) + int(pos)
            pnext = int(start) + int(pnext)
            rest = '\t'.join(rest)
            before_cluster = dic_qname[name]
            after_cluster = dic_barcode_all[before_cluster]
            if after_cluster in dic_barcode_overlap:
                number_success += 1
            else:
                after_cluster = 'NNNNNNNNNNNN'
            if int(tag) < 128:
                name += ':' + chromo + '-0-' + str(pos) + '-' + after_cluster + ':' + before_cluster
            else:
                name += ':' + chromo + '-1-' + str(pos) + '-' + after_cluster + ':' + before_cluster
            new = name + '\t' + tag + '\t' + chromo + '\t' + str(pos) + '\t' + a1 + '\t' + a2 + '\t' + a3 + '\t' + str(
                pnext) + '\t' + a4 + '\t' + a5 + '\t' + a6 + '\t' + rest
            out.write(new)
    samfile.close()
    out.close()
    print('%s reads have barcode.'%number_success)

def main():
    print ('start time is %s.'%time.ctime())
    (consolidate_read1, consolidate_read2, overlap_file, all_barcode_file, filter_file, output)=sys.argv[1:]
    r1r2_overlap(consolidate_read1,consolidate_read2,overlap_file)
    dic_qname, dic_barcode_all, dic_barcode_overlap = make_dic_barcode(all_barcode_file,overlap_file)
    format_sam(filter_file,dic_qname, dic_barcode_all, dic_barcode_overlap,output)
    print ('end time is %s.'%time.ctime())

if __name__ == '__main__':
    main()

