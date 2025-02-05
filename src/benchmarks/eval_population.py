import sys
import argparse
import logging
import time


def hwe_exchet(population_vcf):
    output = open("hwe_exchet.csv", "w")
    output.write("MissingRate\tHWE\tExcHet\n")
    hwe_list = [0 for i in range(22)]
    exchet_list = [0 for i in range(22)]
    with open(population_vcf) as f:
        for line in f:
            if line[0] == "#":
                continue
            else:
                seq = line.strip().split("\t")
                ac = int(seq[7].split(";AC=")[1].split(";")[0])
                hwe = float(seq[7].split(";HWE=")[1].split(";")[0])
                exchet = float(seq[7].split(";ExcHet=")[1])
                if ac == 0:
                    hwe_list[20] += 1
                    exchet_list[20] += 1
                    continue
                missing_cnt = 0
                for i in range(9, 109, 1):
                    if seq[i][0] == ".":
                        missing_cnt += 1
                    if seq[i][2] == ".":
                        missing_cnt += 1
                if missing_cnt > 10:
                    hwe_list[21] += 1
                    exchet_list[21] += 1
                for i in range(20):
                    if i / 20 < hwe <= (i + 1) / 20:
                        hwe_list[i] += 1
                    if i / 20 < exchet <= (i + 1) / 20:
                        exchet_list[i] += 1
                output.write("%f\t%f\t%f\n" % (missing_cnt / 200, hwe, exchet))


def hete(population_vcf):
    output = open("hete.csv", "w")
    output.write("AF\tHete\n")
    with open(population_vcf) as f:
        for line in f:
            if line[0] == "#":
                continue
            else:
                seq = line.strip().split("\t")
                try:
                    af = float(seq[7].split(";AF=")[1].split(";")[0])
                except Exception:
                    ac = 0
                    an = 0
                    for i in range(9, 109, 1):
                        if seq[i][0] == "1":
                            ac += 1
                            an += 1
                        if seq[i][2] == "1":
                            ac += 1
                            an += 1
                        if seq[i][0] == "0":
                            an += 1
                        if seq[i][2] == "0":
                            an += 1
                    if an == 0:
                        af = 0.0
                    else:
                        af = ac / an
                if af <= 0.0001:
                    continue
                hete = 0
                for i in range(9, 109, 1):
                    if seq[i][0] == "0" and seq[i][2] == "1":
                        hete += 1
                output.write("%f\t%f\n" % (af, hete / 100))


def compare(base_file, comp_file):
    def parse(file):
        svs_dict = dict()
        with open(file, "r") as f:
            for line in f:
                if line[0] == "#":
                    continue
                seq = line.strip().split("\t")
                chrom = seq[0]
                pos = int(seq[1])
                svtype = seq[7].split("SVTYPE=")[1].split(";")[0]
                if svtype == "DEL" or svtype == "INS":
                    svlen = abs(int(seq[7].split("SVLEN=")[1].split(";")[0]))
                    af = float(seq[7].split(";AF=")[1].split(";")[0])
                    if chrom not in svs_dict:
                        svs_dict[chrom] = list()
                    svs_dict[chrom].append([pos, svtype, svlen, af])
        return svs_dict

    base = parse(base_file)
    comp = parse(comp_file)
    ans = dict()
    output = open("difference.csv", "w")
    result_5 = []
    pos_bias = 1000
    length_ratio = 0.7
    for chrom in base:
        if chrom in comp:
            ans[chrom] = list()
            for basesv in base[chrom]:
                for compsv in comp[chrom]:
                    if (
                        basesv[1] == compsv[1]
                        and abs(basesv[0] - compsv[0]) <= pos_bias
                        and min(basesv[2], compsv[2]) / max(basesv[2], compsv[2])
                        > length_ratio
                    ):
                        ans[chrom].append(
                            (basesv[1], basesv[0], compsv[0], basesv[3], compsv[3])
                        )
                        output.write(
                            "%s\t%f\t%f\t%f\n"
                            % (basesv[1], basesv[3], compsv[3], basesv[3] - compsv[3])
                        )
                        if abs(basesv[3] - compsv[3]) > 0.5:
                            result_5.append([basesv[1], chrom, basesv[0]])
                        break
    return result_5


def eval_known_callset(sniffles, sniffles2, cutesv, base):
    sniffles = compare(base, sniffles)
    sniffles2 = compare(base, sniffles2)
    cutesv = compare(base, cutesv)
    ss2 = 0
    sc = 0
    s2c = 0
    ss2c = 0
    for item1 in sniffles:
        for item2 in sniffles2:
            if item1[1] == item2[1] and item1[0] == item2[0] and item1[2] == item2[2]:
                ss2 += 1
        for item2 in cutesv:
            if item1[1] == item2[1] and item1[0] == item2[0] and item1[2] == item2[2]:
                sc += 1
    for item1 in cutesv:
        for item2 in sniffles2:
            if item1[1] == item2[1] and item1[0] == item2[0] and item1[2] == item2[2]:
                s2c += 1
                for item3 in sniffles:
                    if (
                        item1[1] == item3[1]
                        and item1[0] == item3[0]
                        and item1[2] == item3[2]
                    ):
                        ss2c += 1
    logging.info(
        "Sniffles1: %d\tSniffles2: %d\t cuteSV2: %d"
        % (len(sniffles), len(sniffles2), len(cutesv))
    )
    logging.info(
        "Sniffles1&Sniffes2: %d\tSniffles1&cuteSV2: %d\tSniffles2&cuteSV2: %d"
        % (ss2, sc, s2c)
    )
    logging.info("Sniffles1&Sniffes2&cuteSV2: %d" % (ss2c))


def main(argv):
    args = parseArgs(argv)
    setupLogging(False)
    starttime = time.time()
    if args.handle == "HWE" or args.handle == "ExcHet":
        hwe_exchet(args.comp)
    if args.handle == "Hete":
        hete(args.comp)
    if args.handle == "Known":
        eval_known_callset(args.sniffles, args.sniffles2, args.cutesv, args.base)
    logging.info("Finished in %0.2f seconds." % (time.time() - starttime))


USAGE = """\
	Evaluate SV callset generated by simulations.
	Author: Shuqi Cao
	Email: sqcao@stu.hit.edu.cn
"""


def parseArgs(argv):
    parser = argparse.ArgumentParser(
        prog="evaluation on HWE/ExcHet/Hete/Known worldwide cohort",
        description=USAGE,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "handle",
        type=str,
        help="The aspect of evaluation, contains HWE/ExcHet/Hete/Known.",
    )
    parser.add_argument("--comp", type=str, help="BND callsets to be benched.")
    parser.add_argument("--base", type=str, help="Ground truth of BNDs.")
    parser.add_argument("--sniffles", type=str, help="Callsets of Sniffles1.")
    parser.add_argument("--sniffles2", type=str, help="Callsets of Sniffles2.")
    parser.add_argument("--cutesv", type=str, help="Callsets of cuteSV2.")
    parser.add_argument(
        "-o",
        "--offect",
        help="Offect of translocation overlaping.[%(default)s]",
        default=1000,
        type=int,
    )
    args = parser.parse_args(argv)
    return args


def setupLogging(debug=False):
    logLevel = logging.DEBUG if debug else logging.INFO
    logFormat = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig(stream=sys.stderr, level=logLevel, format=logFormat)
    logging.info("Running %s" % " ".join(sys.argv))


if __name__ == "__main__":
    main(sys.argv[1:])
