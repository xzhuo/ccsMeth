import os
import argparse
import sys
import time
import numpy as np
from statsmodels import robust
from subprocess import Popen, PIPE
import multiprocessing as mp
from multiprocessing import Queue
import re
import random
# from collections import Counter

from utils.process_utils import display_args
from utils.process_utils import codecv1_to_frame
from utils.process_utils import generate_samtools_view_cmd
from utils.process_utils import get_refloc_of_methysite_in_motif
from utils.process_utils import get_motif_seqs
from utils.process_utils import complement_seq

code2frames = codecv1_to_frame()
queen_size_border = 1000
time_wait = 1

exceptval = 1000
subreads_value_default = "-"


def check_input_file(inputfile):
    if not (inputfile.endswith(".bam") or inputfile.endswith(".sam")):
        raise ValueError("--input/-i must be in bam/sam format!")
    inputpath = os.path.abspath(inputfile)
    return inputpath


def check_output_file(outputfile, inputfile):
    if outputfile is None:
        fname, fext = os.path.splitext(inputfile)
        output_path = fname + ".features.tsv"
    else:
        output_path = os.path.abspath(outputfile)
    return output_path


def cmd_get_stdout_of_input(inputpath, path_to_samtools):
    if inputpath.endswith(".bam"):
        samtools_view = generate_samtools_view_cmd(path_to_samtools)
        cmd = samtools_view + " " + inputpath
    elif inputpath.endswith(".sam"):
        cmd = "cat " + inputpath
    else:
        raise ValueError()
    return cmd


def _get_holeid(subread_id):
    words = subread_id.strip().split("/")
    holeid = words[0] + "/" + words[1]
    return holeid


def _ccs_words_to_feature(words, args):
    query, seq, tags = [words[i] for i in (0, 9, 11)]
    tags_dict = dict(re.split(':\S:', x) for x in tags.split("\t"))
    lib, id, ccs = query.split("/")
    holeid = lib + "/" + id
    fi, fp, ri, rp = [], [], [], []
    feature_strs = []

    try:
        fi = [int(x) for x in tags_dict['fi'].split(",")[1:]]
        ri = [int(x) for x in tags_dict['ri'].split(",")[1:]]
        fp = [int(x) for x in tags_dict['fp'].split(",")[1:]]
        rp = [int(x) for x in tags_dict['rp'].split(",")[1:]]
        fn = int(tags_dict['fn'])
        rn = int(tags_dict['rn'])
        if fn < args.depth or rn < args.depth:
            return feature_strs

        if not args.no_decode:
            fi = [code2frames[ipdval] for ipdval in fi]
            ri = [code2frames[ipval] for ipval in ri]
            fp = [code2frames[pwval] for pwval in fp]
            rp = [code2frames[pwval] for pwval in rp]
            fi = _normalize_signals(fi, args.norm)
            ri = _normalize_signals(ri, args.norm)
            fp = _normalize_signals(fp, args.norm)
            rp = _normalize_signals(rp, args.norm)
            ri = ri[::-1]
            rp = rp[::-1]

        kmer_width = args.seq_len
        half_width = kmer_width // 2
        motifs = get_motif_seqs(args.motifs)
        len_seq = len(seq)
        for regex in motifs:
            for match in re.finditer(regex, seq):
                if match.start() < half_width or len_seq - match.start() < half_width + 1:
                    continue
                cg_index = match.start()
                kmer_start = cg_index - half_width
                kmer_end = kmer_start + kmer_width
                kmer = seq[kmer_start:kmer_end]
                pos_kmer = kmer[:-1]
                neg_kmer = revcom(kmer[1:])
                repeatfn = str(fn) * kmer_width
                repeatrn = str(rn) * kmer_width
                repeatstd = "0" * kmer_width

                kmer_feature = (ccs, cg_index, "+", holeid, max(fn, rn), 
                                pos_kmer, list(repeatfn), fi[kmer_start:kmer_end], list(repeatstd), fp[kmer_start:kmer_end], list(repeatstd), "-", "-", 
                                neg_kmer, list(repeatrn), ri[kmer_start:kmer_end], list(repeatstd), rp[kmer_start:kmer_end], list(repeatstd), "-", "-",
                                args.methy_label)
                feature_str = _features_to_str_combedfeatures(kmer_feature)
                feature_strs.append(feature_str)

        return feature_strs
    except Exception:
        return feature_strs


def revcom(seq):
    tab = str.maketrans("ACGT", "TGCA")
    return seq.translate(tab)[::-1]


def worker_read(inputfile, readline_q, args):
    sys.stderr.write("read_input process-{} starts\n".format(os.getpid()))
    cmd_view_input = cmd_get_stdout_of_input(inputfile, args.path_to_samtools)
    sys.stderr.write("cmd to view input: {}\n".format(cmd_view_input))
    proc_read = Popen(cmd_view_input, shell=True, stdout=PIPE)
    cnt_holes = 0
    readline_list = []
    while True:
        output = str(proc_read.stdout.readline(), 'utf-8')
        if output != "":
            cnt_holes += 1
            readline_list.append(output)
            if len(readline_list) >= args.holes_batch:
                readline_q.put(readline_list)
                readline_list = []
                while readline_q.qsize() > queen_size_border:
                    time.sleep(time_wait)
        elif proc_read.poll() is not None:
            break
        else:
            continue

    readline_q.put(readline_list)
    while readline_q.qsize() > queen_size_border:
        time.sleep(time_wait)

    readline_q.put("kill")
    rc_read = proc_read.poll()
    sys.stderr.write("read_input process-{} ending, read {} holes, with return_code-{}\n".format(os.getpid(),
                                                                                                 cnt_holes,
                                                                                                 rc_read))


def _ccs_extract(readline_q, featurestr_q, args, holeids_e=None, holeids_ne=None):
    sys.stderr.write("extrac_features process-{} starts\n".format(os.getpid()))
    cnt_linebatch = 0
    while True:
        # print("hole_align_q size:", hole_align_q.qsize(), "; pid:", os.getpid())
        if readline_q.empty():
            time.sleep(time_wait)
            continue
        readline_list = readline_q.get()
        if readline_list == "kill":
            readline_q.put("kill")
            break
        feature_strs = []

        for output in readline_list:
            try:
                if output.startswith("#") or output.startswith("@"):
                    continue
                words = output.strip().split('\t', 11)
                holeid = _get_holeid(words[0])
                if holeids_e is not None and holeid not in holeids_e:
                    continue
                if holeids_ne is not None and holeid in holeids_ne:
                    continue

                # flag = int(words[1])
                # mapq = int(words[4])
                # if not (flag == 0 or flag == 16):  # skip segment alignment
                #     continue
                # if mapq < args.mapq:  # skip low mapq alignment
                #     continue
                # now working on words:
                kmer_feature_list = _ccs_words_to_feature(words, args)
                for kmer_feature in kmer_feature_list:
                    feature_strs.append(kmer_feature)
            #     if holeid != holeid_curr:
            #         if len(hole_align_tmp) > 0:
            #             cnt_holes += 1
            #             holes_align_tmp.append((holeid_curr, hole_align_tmp))
            #             if len(holes_align_tmp) >= args.holes_batch:
            #                 hole_align_q.put(holes_align_tmp)
            #                 holes_align_tmp = []
            #                 while hole_align_q.qsize() > queen_size_border:
            #                     time.sleep(time_wait)
            #         hole_align_tmp = []
            #         holeid_curr = holeid
            #     hole_align_tmp.append(words)
            except Exception:
                raise ValueError("error in extracting lines of input!")
                continue

        featurestr_q.put(feature_strs)
        while featurestr_q.qsize() > queen_size_border:
            time.sleep(time_wait)
        cnt_linebatch += 1
        if cnt_linebatch % 200 == 0:
            sys.stderr.write("extrac_features process-{}, {} hole_batches({}) "
                             "proceed\n".format(os.getpid(), cnt_linebatch, args.holes_batch))
            sys.stderr.flush()

    sys.stderr.write("extrac_features process-{} ending, proceed {} "
                     "hole_batches({})\n".format(os.getpid(), cnt_linebatch, args.holes_batch))


def _normalize_signals(signals, normalize_method="zscore"):
    if normalize_method == 'zscore':
        sshift, sscale = np.mean(signals), np.std(signals)
    elif normalize_method == 'min-max':
        sshift, sscale = np.min(signals), np.max(signals) - np.min(signals)
    elif normalize_method == 'min-mean':
        sshift, sscale = np.min(signals), np.mean(signals)
    elif normalize_method == 'mad':
        sshift, sscale = np.median(signals), np.float(robust.mad(signals))
    else:
        raise ValueError("")
    if sscale == 0.0:
        norm_signals = signals
    else:
        norm_signals = (signals - sshift) / sscale
    return np.around(norm_signals, decimals=6)


def check_excpval(myarray):
    if exceptval in set(myarray):
        return True
    return False


def _features_to_str_combedfeatures(features):
    """

    :param features: a tuple
    :return:
    """
    chrom, abs_loc, strand, holeid, depth_all, \
        kmer_seq, kmer_depth, kmer_ipdm, kmer_ipds, kmer_pwm, kmer_pws, kmer_subr_ipds, kmer_subr_pws, \
        kmer_seq2, kmer_depth2, kmer_ipdm2, kmer_ipds2, kmer_pwm2, kmer_pws2, kmer_subr_ipds2, kmer_subr_pws2, \
        label = features
    kmer_depth_str = ",".join([str(x) for x in kmer_depth])
    kmer_ipdm_str = ",".join([str(x) for x in kmer_ipdm])
    kmer_ipds_str = ",".join([str(x) for x in kmer_ipds])
    kmer_pwm_str = ",".join([str(x) for x in kmer_pwm])
    kmer_pws_str = ",".join([str(x) for x in kmer_pws])
    if kmer_subr_ipds != subreads_value_default:
        kmer_subr_ipds_str = ";".join([",".join([str(x) for x in y]) for y in kmer_subr_ipds])
        kmer_subr_pws_str = ";".join([",".join([str(x) for x in y]) for y in kmer_subr_pws])
    else:
        kmer_subr_ipds_str = kmer_subr_ipds
        kmer_subr_pws_str = kmer_subr_pws

    kmer_depth_str2 = ",".join([str(x) for x in kmer_depth2])
    kmer_ipdm_str2 = ",".join([str(x) for x in kmer_ipdm2])
    kmer_ipds_str2 = ",".join([str(x) for x in kmer_ipds2])
    kmer_pwm_str2 = ",".join([str(x) for x in kmer_pwm2])
    kmer_pws_str2 = ",".join([str(x) for x in kmer_pws2])
    if kmer_subr_ipds2 != subreads_value_default:
        kmer_subr_ipds_str2 = ";".join([",".join([str(x) for x in y]) for y in kmer_subr_ipds2])
        kmer_subr_pws_str2 = ";".join([",".join([str(x) for x in y]) for y in kmer_subr_pws2])
    else:
        kmer_subr_ipds_str2 = kmer_subr_ipds2
        kmer_subr_pws_str2 = kmer_subr_pws2

    return "\t".join([chrom, str(abs_loc), strand, str(holeid), str(depth_all),
                      kmer_seq, kmer_depth_str, kmer_ipdm_str, kmer_ipds_str, kmer_pwm_str, kmer_pws_str,
                      kmer_subr_ipds_str, kmer_subr_pws_str,
                      kmer_seq2, kmer_depth_str2, kmer_ipdm_str2, kmer_ipds_str2, kmer_pwm_str2, kmer_pws_str2,
                      kmer_subr_ipds_str2, kmer_subr_pws_str2,
                      str(label)])


def _write_featurestr_to_file(write_fp, featurestr_q):
    sys.stderr.write('write_process-{} started\n'.format(os.getpid()))
    with open(write_fp, 'w') as wf:
        while True:
            # during test, it's ok without the sleep(time_wait)
            if featurestr_q.empty():
                time.sleep(time_wait)
                continue
            features_str = featurestr_q.get()
            if features_str == "kill":
                sys.stderr.write('write_process-{} finished\n'.format(os.getpid()))
                break
            for one_features_str in features_str:
                wf.write(one_features_str + "\n")
            wf.flush()


def _get_holes(holeidfile):
    holes = set()
    with open(holeidfile, "r") as rf:
        for line in rf:
            words = line.strip().split("\t")
            holeid = words[0]
            holes.add(holeid)
    sys.stderr.write("get {} holeids from {}\n".format(len(holes), holeidfile))
    return holes


def extract_ccs_features(args):
    sys.stderr.write("[extract_features]start..\n")
    start = time.time()

    inputpath = check_input_file(args.input)
    outputpath = check_output_file(args.output, inputpath)

    if not os.path.exists(inputpath):
        raise IOError("input file does not exist!")

    if args.seq_len % 2 == 0:
        raise ValueError("seq_len must be odd")

    holeids_e = None if args.holeids_e is None else _get_holes(args.holeids_e)
    holeids_ne = None if args.holeids_ne is None else _get_holes(args.holeids_ne)

    readline_q = Queue()
    featurestr_q = Queue()

    p_read = mp.Process(target=worker_read, args=(inputpath, readline_q, args))
    p_read.daemon = True
    p_read.start()

    ps_extract = []
    nproc = args.threads
    if nproc == 2:
        nproc -= 1
    if nproc > 2:
        nproc -= 2
    for _ in range(nproc):
        p = mp.Process(target=_ccs_extract, args=(readline_q, featurestr_q, args, holeids_e, holeids_ne))
        p.daemon = True
        p.start()
        ps_extract.append(p)

    # print("write_process started..")
    p_w = mp.Process(target=_write_featurestr_to_file, args=(outputpath, featurestr_q))
    p_w.daemon = True
    p_w.start()

    while True:
        # print("killing _worker_extract process")
        running = any(p.is_alive() for p in ps_extract)
        if not running:
            break

    for p in ps_extract:
        p.join()
    p_read.join()

    # sys.stderr.write("finishing the write_process..\n")
    featurestr_q.put("kill")
    p_w.join()

    endtime = time.time()
    sys.stderr.write("[extract_features]costs {:.1f} seconds\n".format(endtime - start))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--threads", type=int, default=5, required=False,
                        help="number of threads, default 5")
    p_input = parser.add_argument_group("INPUT")
    p_input.add_argument("--input", "-i", type=str, required=True,
                         help="alignment results in bam/sam format. "
                              "We assume that all items/reads are sorted by hole_ids "
                              "in aligned.bam, which generated by align_subreads.py from subreads.bam.")
    p_input.add_argument("--holeids_e", type=str, default=None, required=False,
                         help="file contains holeids to be extracted, default None")
    p_input.add_argument("--holeids_ne", type=str, default=None, required=False,
                         help="file contains holeids not to be extracted, default None")

    p_output = parser.add_argument_group("OUTPUT")
    p_output.add_argument("--output", "-o", type=str, required=False,
                          help="output file path to save the extracted features. "
                               "If not specified, use input_prefix.tsv as default.")

    p_extract = parser.add_argument_group("EXTRACT")
    p_extract.add_argument("--seq_len", type=int, default=21, required=False,
                           help="len of kmer. default 21")
    p_extract.add_argument("--motifs", action="store", type=str,
                           required=False, default='CG',
                           help='motif seq to be extracted, default: CG. '
                                'can be multi motifs splited by comma '
                                '(no space allowed in the input str), '
                                'or use IUPAC alphabet, '
                                'the mod_loc of all motifs must be '
                                'the same')
    p_extract.add_argument("--mod_loc", action="store", type=int, required=False, default=0,
                           help='0-based location of the targeted base in the motif, default 0')
    p_extract.add_argument("--methy_label", action="store", type=int,
                           choices=[1, 0], required=False, default=1,
                           help="the label of the interested modified bases, this is for training."
                                " 0 or 1, default 1")
    p_extract.add_argument("--mapq", type=int, default=20, required=False,
                           help="MAPping Quality cutoff for selecting alignment items, default 20")
    p_extract.add_argument("--identity", type=float, default=0.8, required=False,
                           help="identity cutoff for selecting alignment items, default 0.8")
    p_extract.add_argument("--depth", type=int, default=1, required=False,
                           help="(mean) depth (number of subreads) cutoff for "
                                "selecting high-quality aligned reads/kmers "
                                "per strand of a CCS, default 1.")
    p_extract.add_argument("--norm", action="store", type=str, choices=["zscore", "min-mean", "min-max", "mad"],
                           default="zscore", required=False,
                           help="method for normalizing ipd/pw in subread level. "
                                "zscore, min-mean, min-max or mad, default zscore")
    p_extract.add_argument("--no_decode", action="store_true", default=False, required=False,
                           help="not use CodecV1 to decode ipd/pw")
    p_extract.add_argument("--path_to_samtools", type=str, default=None, required=False,
                           help="full path to the executable binary samtools file. "
                                "If not specified, it is assumed that samtools is in "
                                "the PATH.")
    p_extract.add_argument("--holes_batch", type=int, default=50, required=False,
                           help="number of holes in an batch to get/put in queues")
    p_extract.add_argument("--seed", type=int, default=1234, required=False,
                           help="seed for randomly selecting subreads, default 1234")

    args = parser.parse_args()

    display_args(args, True)
    extract_ccs_features(args)


if __name__ == '__main__':
    main()
