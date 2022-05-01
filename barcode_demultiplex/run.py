import os
import click
import pandas as pd
import subprocess
import glob
import sys
import logging
import numpy as np
import dreem.run
import shutil
import pickle

import rna_library as rl
from seq_tools.sequence import get_reverse_complement

np.warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)
APP_LOGGER_NAME = "barcode-demultiplex"


def get_default_dreem_args():
    args = {
        "fasta": "",
        "fastq1": "",
        "fastq2": None,
        "dot_bracket": None,
        "param_file": None,
        "overwrite": False,
        "log_level": "INFO",
        "restore_org_behavior": False,
        "map_overwrite": False,
        "skip": False,
        "skip_fastqc": False,
        "skip_trim_galore": False,
        "bt2_alignment_args": None,
        "bv_overwrite": False,
        "qscore_cutoff": None,
        "num_of_surbases": None,
        "map_score_cutoff": None,
        "mutation_count_cutoff": None,
        "percent_length_cutoff": None,
        "summary_output_only": False,
        "plot_sequence": False,
    }
    return args


def setup_applevel_logger(
    logger_name=APP_LOGGER_NAME, is_debug=True, file_name=None
):
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG if is_debug else logging.INFO)

    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(formatter)
    logger.handlers.clear()
    logger.addHandler(sh)

    if file_name:
        fh = logging.FileHandler(file_name)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger


log = setup_applevel_logger()


def get_read_length(fastq):
    f = open(fastq)
    lines = f.readlines()
    f.close()
    min = 10000
    for i in range(0, len(lines), 4):
        if len(lines[i + 1]) < min:
            min = len(lines[i + 1])
    return min


class Grepdemultiplexer(object):
    def __init__(self):
        pass

    def run(self, df, hpos, fastq1, fastq2):
        for i, row in df.iterrows():
            s = rl.SecStruct(
                row["structure"], row["sequence"].replace("U", "T")
            )
            helices = []
            for m in s:
                if m.is_helix():
                    helices.append(m)
            path = f"data/{row['name']}/"
            for hp in hpos:
                minpos = helices[hp].strands()[0][0]
                maxpos = helices[hp].strands()[0][-1]
                spl = helices[hp].sequence().split("&")
                cmd = f'cut -c {minpos - 10}-{maxpos + 10} {fastq1} | grep -n "{spl[0]}" > {path}/hp_{hp}.out'
                # print(cmd)
                # subprocess.call(cmd, shell=True)
                nums = self.__get_line_nums(glob.glob(path + "/*.out"))
                # self.__generate_sub_fastqs(nums, fastq1, fastq2, f"data/{row['name']}")
                os.chdir(path)
                subprocess.call(
                    "dreem -fa test.fasta -fq1 test_mate1.fastq -fq2 test_mate2.fastq --dot_bracket test.csv",
                    shell=True,
                )
                os.chdir("../..")

    def __get_nums_from_file(self, file):
        nums = []
        f = open(file)
        lines = f.readlines()
        for l in lines:
            spl = l.split(":")
            nums.append(int(spl[0]))
        f.close()
        return set(nums)

    def __get_line_nums(self, files):
        sets = []
        for file in files:
            s = self.__get_nums_from_file(file)
            sets.append(s)
        s = sets[0]
        for i in range(1, len(sets)):
            s = s.intersection(sets[i])
        return list(s)

    def __generate_sub_fastqs(self, nums, fastq1, fastq2, outputdir):
        os.makedirs(outputdir, exist_ok=True)
        f = open(fastq1)
        lines_r1 = f.readlines()
        f.close()
        f = open(fastq2)
        lines_r2 = f.readlines()
        f.close()
        f1 = open(f"{outputdir}/test_mate1.fastq", "w")
        f2 = open(f"{outputdir}/test_mate2.fastq", "w")
        for i in range(0, len(lines_r1), 4):
            if i + 2 in nums:
                f1.writelines(lines_r1[i : i + 4])
                f2.writelines(lines_r2[i : i + 4])
        f1.close()
        f2.close()


class Seqkitdemultiplexer(object):
    def __init__(self):
        self.check_fwd = True
        self.check_rev = True
        self.r1_length = -1
        self.r2_length = -1
        self.r1_lines = []
        self.r2_lines = []
        self.all_ids = []
        self.max = 9999

    def run(self, df, fastq1, fastq2):
        self.r1_length = get_read_length(fastq1)
        self.r2_length = get_read_length(fastq2)
        f = open(fastq1)
        self.r1_lines = f.readlines()
        f.close()
        f = open(fastq2)
        self.r2_lines = f.readlines()
        f.close()
        log.info(f"total number of reads: {len(self.r1_lines) / 4}")
        self.reads_used = 0
        count = 0
        for i, row in df.iterrows():
            if count >= self.max:
                log.info(f"max number of barcodes processed: {self.max}")
                break
            if self.check_fwd:
                self.__run_fwd_demultiplex(row, fastq1)
            if self.check_rev:
                self.__run_rev_demultiplex(row, fastq2)
            self.__filter_reads(row["path"])
            count += 1
        log.info(f"total reads found: {self.reads_used}")
        # self.__output_unused_reads(fastq1)

    def __run_fwd_demultiplex(self, row, fastq1):
        cmds = []
        count = 0
        for seqs, bounds in zip(row["barcodes"], row["barcode_bounds"]):
            for seq, bound in zip(seqs, bounds):
                b1, b2 = bound
                if b2 + 10 >= self.r1_length:
                    continue
                count += 1
                if count == 1:
                    cmds.append(
                        f'seqkit grep -s -p "{seq}" -m 1 -P -R {b1 - 10}:{b2 + 10} {fastq1} '
                    )
                else:
                    cmds.append(
                        f'seqkit grep -s -p "{seq}" -m 1 -P -R {b1 - 10}:{b2 + 10} '
                    )
        cmd = "| ".join(cmds) + f"> {row['path']}/test_fwd.fastq"
        subprocess.call(cmd, shell=True)

    def __run_rev_demultiplex(self, row, fastq2):
        cmds = []
        count = 0
        for seqs, bounds in zip(row["barcodes"], row["barcode_bounds"]):
            for seq, bound in zip(seqs, bounds):
                b1, b2 = bound
                b3 = row["len"] - b1
                b4 = row["len"] - b2
                if b3 + 10 >= self.r2_length:
                    continue
                count += 1
                rev_seq = get_reverse_complement(seq)
                if count == 1:
                    cmds.append(
                        f'seqkit grep -s -p "{rev_seq}" -m 1 -P -R {b4 - 10}:{b3 + 10} {fastq2}'
                    )
                else:
                    cmds.append(
                        f'seqkit grep -s -p "{rev_seq}" -m 1 -P -R {b4 - 10}:{b3 + 10} '
                    )

        cmd = "| ".join(cmds) + f"> {row['path']}/test_rev.fastq"
        subprocess.call(cmd, shell=True)

    def __get_read_ids(self, fpath):
        f = open(fpath)
        lines = f.readlines()[::4]
        f.close()
        ids = []
        for l in lines:
            spl = l.split()
            ids.append(spl[0])
        return set(ids)

    def __filter_reads(self, path):
        ids = set()
        if self.check_fwd and self.check_rev:
            fwd_ids = self.__get_read_ids(f"{path}/test_fwd.fastq")
            rev_ids = self.__get_read_ids(f"{path}/test_rev.fastq")
            ids = fwd_ids.intersection(rev_ids)
        elif self.check_fwd:
            ids = self.__get_read_ids(f"{path}/test_fwd.fastq")
        elif self.check_rev:
            ids = self.__get_read_ids(f"{path}/test_rev.fastq")
        self.reads_used += len(ids)
        for id in ids:
            if id not in self.all_ids:
                self.all_ids.append(id)
        log.info(f"{path}: {len(ids)}")
        f = open(f"{path}/test_mate1.fastq", "w")
        for i in range(0, len(self.r1_lines), 4):
            spl = self.r1_lines[i].split()
            if spl[0] not in ids:
                continue
            f.writelines(self.r1_lines[i : i + 4])
        f.close()
        f = open(f"{path}/test_mate2.fastq", "w")
        for i in range(0, len(self.r2_lines), 4):
            spl = self.r2_lines[i].split()
            if spl[0] not in ids:
                continue
            f.writelines(self.r2_lines[i : i + 4])
        f.close()

    def __output_unused_reads(self, fastq1):
        pos_ids = self.__get_read_ids(fastq1)
        ids = []
        for id in pos_ids:
            if id not in self.all_ids:
                ids.append(id)
        f = open(f"test_mate1.fastq", "w")
        for i in range(0, len(self.r1_lines), 4):
            spl = self.r1_lines[i].split()
            if spl[0] not in ids:
                continue
            f.writelines(self.r1_lines[i : i + 4])
        f.close()
        f = open(f"test_mate2.fastq", "w")
        for i in range(0, len(self.r2_lines), 4):
            spl = self.r2_lines[i].split()
            if spl[0] not in ids:
                continue
            f.writelines(self.r2_lines[i : i + 4])
        f.close()


def setup_directories(df, dirname="data"):
    os.makedirs(dirname, exist_ok=True)
    bc = 0
    data = []
    for i, g in df.groupby("full_barcode"):
        bc_dir = f"{dirname}/bc-{bc:03d}"
        os.makedirs(bc_dir, exist_ok=True)
        g.to_json(f"{bc_dir}/constructs.json", orient="records")
        f = open(f"{bc_dir}/test.fasta", "w")
        for j, row in g.iterrows():
            f.write(f">{row['name']}\n")
            f.write(row["sequence"].replace("U", "T") + "\n")
        f.close()
        db_file_df = g[["name", "sequence", "structure"]]
        db_file_df.to_csv(f"{bc_dir}/test.csv", index=False)
        row = g.iloc[0]
        data_row = [
            bc_dir,
            len(g),
            row["barcodes"],
            row["barcode_bounds"],
            row["full_barcode"],
            len(row["sequence"]),
        ]
        data.append(data_row)
        bc += 1
    return pd.DataFrame(
        data,
        columns="path,num,barcodes,barcode_bounds,full_barcode,len".split(","),
    )


def find_barcodes(df, helices):
    def __get_subsection(h1, h2, pos1, pos2):
        h1_new = h1[pos1 : pos2 + 1]
        h2_new = h2[::-1][pos1 : pos2 + 1][::-1]
        return [h1_new, h2_new]

    df["barcodes"] = [[] for _ in range(len(df))]
    df["barcode_bounds"] = [[] for _ in range(len(df))]
    df["full_barcode"] = ""
    for i, row in df.iterrows():
        s = rl.SecStruct(row["structure"], row["sequence"].replace("U", "T"))
        row_helices = list(s.helix())
        all_barcodes = []
        all_bounds = []
        for j, h in enumerate(helices):
            row_h = row_helices[h[0]]
            seqs = row_h.sequence().split("&")
            strands = row_h.strands()
            b_seq = __get_subsection(seqs[0], seqs[1], h[1], h[2])
            b_strands = __get_subsection(strands[0], strands[1], h[1], h[2])
            b_bounds = [
                [min(b_strands[0]), max(b_strands[0])],
                [min(b_strands[1]), max(b_strands[1])],
            ]
            all_barcodes.append(b_seq)
            all_bounds.append(b_bounds)
        full_barcode = "_".join(np.concatenate(all_barcodes).flat)
        df.at[i, ["barcodes", "barcode_bounds", "full_barcode"]] = [
            all_barcodes,
            all_bounds,
            full_barcode,
        ]
    df.to_json("test.json", orient="records")
    return df


def run_dreem_prog(df, max):
    count = 0
    all_mhs = {}
    for i, row in df.iterrows():
        if count >= max:
            break
        args = get_default_dreem_args()
        path = row["path"]
        args["fasta"] = path + "/test.fasta"
        args["fastq1"] = path + "/test_mate1.fastq"
        args["fastq2"] = path + "/test_mate2.fastq"
        args["dot_bracket"] = path + "/test.csv"
        args["summary_output_only"] = True
        dreem.run.run(args)
        pickle_path = "output/BitVector_Files/mutation_histos.p"
        mhs = pickle.load(open(pickle_path, "rb"))
        for name, mh in mhs.items():
            if name in all_mhs:
                raise ValueError(
                    f"something went very wrong should not be duplicate entries"
                    f" in MutationalHistograms! {name} is duplicated"
                )
            all_mhs[name] = mh
        log.info(f"currently {len(all_mhs)} mutational histograms!")
        shutil.move(pickle_path, path)
        shutil.rmtree("input")
        shutil.rmtree("log")
        shutil.rmtree("output")
        count += 1
    os.makedirs("output/BitVector_Files/", exist_ok=True)
    pickle.dump(all_mhs, open("output/BitVector_Files/mutation_histos.p", "wb"))



# cli #########################################################################


@click.command()
@click.option("-csv", "--rna-csv", required=True)
@click.option("-fq1", "--fastq1", required=True)
@click.option("-fq2", "--fastq2", required=False)
@click.option("-rd", "--run-dreem", is_flag=True)
@click.option("-m", "--max", default=99999)
@click.option(
    "-helix",
    "--helix",
    "helices",
    type=click.Tuple([int, int, int]),
    required=True,
    multiple=True,
)
def main(rna_csv, fastq1, fastq2, helices, run_dreem, max):
    df = pd.read_csv(rna_csv)
    df = find_barcodes(df, helices)
    df_sum = setup_directories(df, "data")
    log.info(f"{len(df_sum)} unique barcodes found!")
    df_sum.to_csv("summary.csv", index=False)
    demult = Seqkitdemultiplexer()
    demult.max = max
    demult.run(df_sum, fastq1, fastq2)
    if run_dreem:
        run_dreem_prog(df_sum, max)


if __name__ == "__main__":
    main()
