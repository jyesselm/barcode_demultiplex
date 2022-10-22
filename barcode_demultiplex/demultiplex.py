import os
import click
import pandas as pd
import subprocess
import numpy as np
import dreem.run
import shutil
import pickle

import rna_library as rl
from seq_tools.sequence import get_reverse_complement
from dreem.run import get_default_run_args

from barcode_demultiplex.logger import *

np.warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)


def get_read_length(fastq):
    f = open(fastq)
    lines = f.readlines()
    f.close()
    min_len = 10000
    for i in range(0, len(lines), 4):
        if len(lines[i + 1]) < min_len:
            min_len = len(lines[i + 1])
    return min_len


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
        self.reads_used = 0

    def run(self, df, fastq1, fastq2):
        log = get_logger("")
        self.r1_length = get_read_length(fastq1)
        f = open(fastq1)
        self.r1_lines = f.readlines()
        f.close()
        if self.check_rev:
            self.r2_length = get_read_length(fastq2)
            f = open(fastq2)
            self.r2_lines = f.readlines()
            f.close()
        log.info(f"total number of reads: {len(self.r1_lines) / 4}")
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
                        f'seqkit grep -s -p "{seq}" -m 1 -P -R {b1 - 10}:{b2 + 10} '
                        f"{fastq1} "
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
                        f'seqkit grep -s -p "{rev_seq}" -m 1 -P -R {b4 - 10}:{b3 + 10} '
                        f"{fastq2}"
                    )
                else:
                    cmds.append(
                        f'seqkit grep -s -p "{rev_seq}" -m 1 -P -R {b4 - 10}:{b3 + 10} '
                    )

        cmd = "| ".join(cmds) + f"> {row['path']}/test_rev.fastq"
        subprocess.call(cmd, shell=True)

    @staticmethod
    def __get_read_ids(fpath):
        f = open(fpath)
        lines = f.readlines()[::4]
        f.close()
        ids = []
        for line in lines:
            spl = line.split()
            ids.append(spl[0])
        return set(ids)

    def __filter_reads(self, path):
        log = get_logger("")
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
        for read_id in ids:
            if id not in self.all_ids:
                self.all_ids.append(read_id)
        log.info(f"{path}: {len(ids)}")
        f = open(f"{path}/test_mate1.fastq", "w")
        for i in range(0, len(self.r1_lines), 4):
            spl = self.r1_lines[i].split()
            if spl[0] not in ids:
                continue
            f.writelines(self.r1_lines[i : i + 4])
        f.close()
        if len(self.r2_lines) == 0:
            return
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
        for read_id in pos_ids:
            if read_id not in self.all_ids:
                ids.append(read_id)
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
        bc_dir = f"{dirname}/bc-{bc:04d}"
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
    return df


def run_dreem_prog(df, max_barcodes, include_all_dreem_outputs):
    # need another file setup for keeping all the output files
    if include_all_dreem_outputs:
        return run_dreem_prog_multi(df, max_barcodes)
    log = get_logger("RUN_DREEM")
    count = 0
    all_mhs = {}
    for i, row in df.iterrows():
        if count >= max_barcodes:
            break
        args = get_default_run_args()
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

def run_dreem_prog_multi(df, max_barcodes):
    log = get_logger("RUN_DREEM")
    count = 0
    for i, row in df.iterrows():
        if count >= max_barcodes:
            break
        args = get_default_run_args()
        path = row["path"]
        args["fasta"] = path + "/test.fasta"
        args["fastq1"] = path + "/test_mate1.fastq"
        args["fastq2"] = path + "/test_mate2.fastq"
        args["dot_bracket"] = path + "/test.csv"
        dir_name = row["construct"] + "_" + row["code"] + "_" + row["data_type"]
        print(dir_name)
        #dreem.run.run(args)
        #count += 1


