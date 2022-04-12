import os
import argparse
import pysam
import array

def revcom(seq):
    tab = str.maketrans("ACGT", "TGCA")
    return seq.translate(tab)[::-1]

def mm_generator(seq, pos_list):
    mm_list = ["C+m"]
    for i, end in enumerate(pos_list):
        start = pos_list[i-1] + 2 if i > 0 else 0
        mm_list.append(str(seq[start:end].upper().count('C')))
    return mm_list

def attach_tags(bam_file, tsv_file, out_file):
    hash = {}
    with open(tsv_file, "r") as rf:
        for line in rf:
            words = line.strip().split("\t")
            query_name = str(words[3]) + "/" + str(words[0])
            pos = int(words[1])
            ml = int(float(words[6]) * 256)
            try:
                hash[query_name]['pos_list'].append(pos)
                hash[query_name]['ml_list'].append(ml)
            except:
                hash[query_name] = {'pos_list': [pos], 'ml_list': [ml]}
    bam = pysam.AlignmentFile(bam_file, threads = 8)
    out = pysam.AlignmentFile(out_file, "wb", template=bam, threads = 8)
    for read in bam.fetch():
        query_name = read.query_name
        if query_name in hash:
            seq = revcom(read.query_sequence) if read.flag & 0X10 else read.query_sequence
            # seq = read.query_sequence
            mm_list = mm_generator(seq, hash[query_name]['pos_list'])
           # ml_tag = ','.join(hash[query_name]['ml_list'])
            mm_tag = ','.join(mm_list) + ";"
            read.set_tag('Mm', mm_tag, 'Z')
            read.set_tag('Ml', array.array('B', hash[query_name]['ml_list']))
        out.write(read)

    out.close()
    bam.close()

def main():
    parser = argparse.ArgumentParser(description='attach methylation results from tsv to Mm and Ml tags in the bam file')
    parser.add_argument('-b', '--bam', type=str, required=True,
                        help='input bam file')
    parser.add_argument('-t', '--tsv', type=str, required=True,
                        help='ccsmeth result tsv file')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output bam file with Ml and Mm tags attached')

    args = parser.parse_args()
    tsv_file = os.path.abspath(args.tsv)
    bam_file = os.path.abspath(args.bam)
    if not os.path.exists(tsv_file):
        raise ValueError("--tsv file does not exist!")
    if not os.path.exists(bam_file):
        raise ValueError("--bam file does not exist!")
    # bamfile = args.bam
    # tsvfile = args.tsv
    outfile = args.out
    attach_tags(bam_file, tsv_file, outfile)


if __name__ == '__main__':
    main()
