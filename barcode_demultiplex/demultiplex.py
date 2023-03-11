import os
import click
import pandas as pd
import subprocess
import numpy as np
import shutil
import pickle

from rna_secstruct import SecStruct
from barcode_demultiplex.external_cmd import (
    does_program_exist,
    SeqkitGrepOpts,
    run_seqkit_grep_fwd,
    run_seqkit_grep_rev,
    run_seqkit_common,
)
from barcode_demultiplex.logger import get_logger


log = get_logger("DEMULTIPLEX")


def get_read_length(fastq):
    f = open(fastq)
    lines = f.readlines()
    f.close()
    min_len = 10000
    for i in range(0, len(lines), 4):
        if len(lines[i + 1]) < min_len:
            min_len = len(lines[i + 1])
    return min_len


class Demultiplexer:
    def setup(self, df_barcode, outdir):
        # check if seqkit is installed
        if not does_program_exist("seqkit"):
            log.error("seqkit not found")
            raise ValueError("seqkit not found")
        setup_directories(df_barcode, outdir)

        self.df_barcode = df_barcode
        self.outdir = outdir
        setup_directories(df_barcode, outdir)

    def run(self, fastq1, fastq2):
        bc = -1
        for i, row in self.df_barcode.iterrows():
            bc += 1
            bc_dir = f"{self.outdir}/bc-{bc:04d}/"
            barcode_seqs = row["barcodes"]
            barcode_bounds = row["barcode_bounds"]
            fwd_seqs = [seqs[0] for seqs in barcode_seqs]
            fwd_bounds = [bounds[0] for bounds in barcode_bounds]
            rev_seqs = [seqs[1] for seqs in barcode_seqs]
            rev_bounds = [bounds[1] for bounds in barcode_bounds]
            opts = SeqkitGrepOpts(
                mut_num=0,
                read_len=151,
                seq_len=len(row["sequence"]),
                buffer=5,
            )
            run_seqkit_grep_fwd(fwd_seqs, fwd_bounds,
                                fastq1, "fwd.fastq.gz", opts)
            run_seqkit_grep_rev(rev_seqs, rev_bounds,
                                fastq2, "rev.fastq.gz", opts)
            run_seqkit_common(
                "fwd.fastq.gz", "rev.fastq.gz", f"{bc_dir}/test_mate1.fastq.gz"
            )
            run_seqkit_common(
                "rev.fastq.gz",
                f"{bc_dir}/test_mate1.fastq.gz",
                f"{bc_dir}/test_mate2.fastq.gz",
            )
            os.remove("fwd.fastq.gz")
            os.remove("rev.fastq.gz")

            # TODO run rna-map here


def setup_directories(df, dirname="data"):
    os.makedirs(dirname, exist_ok=True)
    bc = 0
    data = []
    for _, g in df.groupby("full_barcode"):
        bc_dir = f"{dirname}/bc-{bc:04d}"
        os.makedirs(bc_dir, exist_ok=True)
        g.to_json(f"{bc_dir}/constructs.json", orient="records")
        f = open(f"{bc_dir}/test.fasta", "w")
        for _, row in g.iterrows():
            f.write(f">{row['name']}\n")
            f.write(row["sequence"].replace("U", "T") + "\n")
        f.close()
        db_file_df = g[["name", "sequence", "structure"]]
        db_file_df.to_csv(f"{bc_dir}/test.csv", index=False)
        row = g.iloc[0]
        data_row = [
            row["name"],
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
        columns="name,path,num,barcodes,barcode_bounds,full_barcode,len".split(
            ","),
    )


def find_helix_barcodes(df, helices):
    """
    Finds the sequence and bounds of helix barcodes in a dataframe of sequences and
    and structures
    :param df: A dtataframe with columns "sequence" and "structure"
    :param helices: A list of tuples of the form (helix_index, start_pos, end_pos)
    :return: A dataframe with the same columns as the input, plus the columns barocodes,
    barcode_bounds, and full_barcode
    """

    def __get_subsection(h1, h2, pos1, pos2):
        h1_new = h1[pos1: pos2 + 1]
        h2_new = h2[::-1][pos1: pos2 + 1][::-1]
        return [h1_new, h2_new]

    df["barcodes"] = [[] for _ in range(len(df))]
    df["barcode_bounds"] = [[] for _ in range(len(df))]
    df["full_barcode"] = ""
    for i, row in df.iterrows():
        s = SecStruct(row["sequence"].replace("U", "T"), row["structure"])
        row_helices = list(s.get_helices())
        all_barcodes = []
        all_bounds = []
        for j, h in enumerate(helices):
            row_h = row_helices[h[0]]
            seqs = row_h.sequence.split("&")
            strands = row_h.strands
            b_seq = __get_subsection(seqs[0], seqs[1], h[1], h[2])
            b_strands = __get_subsection(strands[0], strands[1], h[1], h[2])
            b_bounds = [
                [min(b_strands[0]), max(b_strands[0])],
                [min(b_strands[1]), max(b_strands[1])],
            ]
            all_barcodes.append(b_seq)
            all_bounds.append(b_bounds)
        full_barcode = "_".join(np.concatenate(all_barcodes).flat)
        df.at[i, "barcodes"] = all_barcodes
        df.at[i, "barcode_bounds"] = all_bounds
        df.at[i, "full_barcode"] = full_barcode
    return df


def run_dreem_prog(df, max_barcodes, include_all_dreem_outputs):
    # need another file setup for keeping all the output files
    if include_all_dreem_outputs:
        return run_dreem_prog_multi(df, max_barcodes)
    count = 0
    all_mhs = {}
    for i, row in df.iterrows():
        if count == max_barcodes:
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
    pickle.dump(all_mhs, open(
        "output/BitVector_Files/mutation_histos.p", "wb"))


def run_dreem_prog_multi(df, max_barcodes):
    count = 0
    for i, row in df.iterrows():
        if count == max_barcodes:
            break
        args = get_default_run_args()
        path = os.path.abspath(row["path"])
        args["fasta"] = path + "/test.fasta"
        args["fastq1"] = path + "/test_mate1.fastq"
        args["fastq2"] = path + "/test_mate2.fastq"
        args["dot_bracket"] = path + "/test.csv"
        os.makedirs(row["full_barcode"], exist_ok=True)
        os.chdir(row["full_barcode"])
        try:
            dreem.run.run(args)
        except:
            pass
        os.chdir("..")
        count += 1
