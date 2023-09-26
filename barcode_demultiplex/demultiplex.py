import os
import pandas as pd
from pathlib import Path
import numpy as np
import yaml
import shutil
import pickle
import gzip
import json

from rna_secstruct import SecStruct
from seq_tools.dataframe import calc_edit_distance

from barcode_demultiplex.external_cmd import (
    does_program_exist,
    SeqkitGrepOpts,
    run_seqkit_grep_fwd,
    run_seqkit_grep_rev,
    run_seqkit_common,
    run_seqkit_stats,
)
from barcode_demultiplex.logger import get_logger

from rna_map.run import run as run_rna_map


log = get_logger("DEMULTIPLEX")


LIB_DIR = Path(os.path.dirname(os.path.realpath(__file__)))


# TODO test what happens when there is an empty fastq file for rna_map
# TODO if is not summary true need to recombine the bitvector files


def get_read_length(fastq_file: Path):
    is_gzipped = fastq_file.suffix == ".gz"
    seq_len = -1
    count = 0
    if not is_gzipped:
        with open(fastq_file) as f:
            for line_number, line in enumerate(f):
                line = line.rstrip()
                if line_number % 4 == 1:
                    count += 1
                    if seq_len < len(line):
                        seq_len = len(line)
    else:
        with gzip.open(fastq_file, "rt") as f:
            for line_number, line in enumerate(f):
                line = line.rstrip()
                if line_number % 4 == 1:
                    count += 1
                    if seq_len < len(line):
                        seq_len = len(line)
    log.info(f"read length is {seq_len} for {fastq_file} (count={count})")
    return seq_len


def min_max_for_sublist(sublists):
    return [min([s[0] for s in sublists]), max([s[1] for s in sublists])]


def get_barcode_bounds_for_construct_group(group):
    num_structures = len(group["barcode_bounds"].iloc[0])
    # For each structure
    results = []
    for i in range(num_structures):
        # Extract the sublists for the current structure across all rows
        current_structure_sublists = [row[i] for row in group["barcode_bounds"]]
        # Transpose to group corresponding sublists
        transposed_lists = list(zip(*current_structure_sublists))
        # Find the min and max values for each group of corresponding sublists
        result_for_current_structure = [
            min_max_for_sublist(sublist) for sublist in transposed_lists
        ]
        results.append(result_for_current_structure)
    return results


class Demultiplexer:
    def setup(self, df_barcode, outdir, params):
        # check if seqkit is installed
        if not does_program_exist("seqkit"):
            log.error("seqkit not found")
            raise ValueError("seqkit not found")
        self.df_dirs = setup_directories(
            df_barcode, params["general"]["directory_style"], outdir
        )
        self.df_barcode = df_barcode
        self.outdir = outdir
        self.params = params
        self.fwd_read_len = 100
        self.rev_read_len = 100

    def __run_seqkit(self, fastq1, fastq2, df_constructs, bc_dir):
        full_barcode = df_constructs.iloc[0]["full_barcode"]
        barcode_seqs = df_constructs.iloc[0]["barcodes"]
        barcode_bounds = get_barcode_bounds_for_construct_group(df_constructs)
        # this assumes only 1 barcode ...
        fwd_seqs = [seqs[0] for seqs in barcode_seqs]
        fwd_bounds = [bounds[0] for bounds in barcode_bounds]
        rev_seqs = [seqs[1] for seqs in barcode_seqs]
        rev_bounds = [bounds[1] for bounds in barcode_bounds]
        fwd_args = self.params["fwd"].copy()
        del fwd_args["check"]
        rev_args = self.params["rev"].copy()
        del rev_args["check"]
        avg_seq_len = np.mean(
            [len(row["sequence"]) for _, row in df_constructs.iterrows()]
        )
        avg_seq_len = int(avg_seq_len)
        fwd_opts = SeqkitGrepOpts(
            read_len=self.fwd_read_len, seq_len=avg_seq_len, **fwd_args
        )
        rev_opts = SeqkitGrepOpts(
            read_len=self.rev_read_len, seq_len=avg_seq_len, **rev_args
        )
        # TODO probably break this up into seperate functions this is too long
        if self.params["general"]["directory_style"] == "individual":
            if self.params["fwd"]["check"]:
                run_seqkit_grep_fwd(
                    fwd_seqs, fwd_bounds, fastq1, f"{bc_dir}/fwd.fastq.gz", fwd_opts
                )
            else:
                shutil.copy(fastq1, f"{bc_dir}/fwd.fastq.gz")
            if self.params["rev"]["check"]:
                run_seqkit_grep_rev(
                    rev_seqs, rev_bounds, fastq2, f"{bc_dir}/rev.fastq.gz", rev_opts
                )
            else:
                shutil.copy(fastq2, f"{bc_dir}/rev.fastq.gz")
            run_seqkit_common(
                f"{bc_dir}/fwd.fastq.gz",
                f"{bc_dir}/rev.fastq.gz",
                f"{bc_dir}/test_mate1.fastq.gz",
            )
            run_seqkit_common(
                f"{bc_dir}/rev.fastq.gz",
                f"{bc_dir}/test_mate1.fastq.gz",
                f"{bc_dir}/test_mate2.fastq.gz",
            )
        else:
            if self.params["fwd"]["check"]:
                run_seqkit_grep_fwd(
                    fwd_seqs,
                    fwd_bounds,
                    fastq1,
                    f"{self.outdir}/fwd.fastq.gz",
                    fwd_opts,
                )
            else:
                shutil.copy(fastq1, f"{self.outdir}/fwd.fastq.gz")
            if self.params["rev"]["check"]:
                run_seqkit_grep_rev(
                    rev_seqs,
                    rev_bounds,
                    fastq2,
                    f"{self.outdir}/rev.fastq.gz",
                    rev_opts,
                )
            else:
                shutil.copy(fastq2, f"{self.outdir}/rev.fastq.gz")
            run_seqkit_common(
                f"{self.outdir}/fwd.fastq.gz",
                f"{self.outdir}/rev.fastq.gz",
                f"{self.outdir}/{full_barcode}_mate1.fastq.gz",
            )
            run_seqkit_common(
                f"{self.outdir}/rev.fastq.gz",
                f"{self.outdir}/{full_barcode}_mate1.fastq.gz",
                f"{self.outdir}/{full_barcode}_mate2.fastq.gz",
            )
            os.remove(f"{self.outdir}/fwd.fastq.gz")
            os.remove(f"{self.outdir}/rev.fastq.gz")

    def __run_rna_map(self, bc_dir, all_mhs, rna_map_params):
        fastq1_path = f"{bc_dir}/test_mate1.fastq.gz"
        fastq2_path = f"{bc_dir}/test_mate2.fastq.gz"
        fa_path = f"{bc_dir}/test.fasta"
        csv_path = f"{bc_dir}/test.csv"
        rna_map_params = rna_map_params.copy()
        # map rna-map directories to the barcode directory
        rna_map_params["dirs"]["output"] = f"{bc_dir}/output"
        rna_map_params["dirs"]["input"] = f"{bc_dir}/input"
        rna_map_params["dirs"]["log"] = f"{bc_dir}/log"
        try:
            run_rna_map(fa_path, fastq1_path, fastq2_path, csv_path, rna_map_params)
        except:
            log.warning(f"rna-map failed for {bc_dir}")
            return
        json_file = f"{bc_dir}/output/BitVector_Files/mutation_histos.json"
        if not os.path.exists(json_file):
            log.warning(f"rna-map did not produce a json file for {bc_dir}")
            return
        with open(json_file) as f:
            mhs = json.load(f)
            num_reads = sum([mh["num_aligned"] for mh in mhs.values()])
            log.info(f"{bc_dir} has {num_reads} aligned reads")
        for name, mh in mhs.items():
            if name in all_mhs:
                raise ValueError(
                    f"something went very wrong should not be duplicate entries"
                    f" in MutationalHistograms! {name} is duplicated"
                )
            all_mhs[name] = mh

    def run(self, fastq1: Path, fastq2: Path):
        # check to make sure files actually exist
        if not fastq1.parts and not fastq1.is_file():
            log.error(f"{fastq1} not found")
            raise ValueError(f"{fastq1} not found")
        if not fastq2.parts and not fastq2.is_file():
            log.error(f"{fastq2} not found")
            raise ValueError(f"{fastq2} not found")

        self.fwd_read_len = get_read_length(fastq1)
        self.rev_read_len = get_read_length(fastq2)
        if self.params["rna-map"]["run"]:
            param_path = LIB_DIR / "resources/rna-map.yml"
            log.info(f"using rna-map params from {param_path}")
            rna_map_params = yaml.safe_load(open(param_path))
            log.info("running rna-map on all barcodes")
            log.info("data will be found in output/BitVector_Files/")
            os.makedirs("output/BitVector_Files/", exist_ok=True)
        all_mhs = {}
        count = 0
        for barcode, group in self.df_barcode.groupby("full_barcode"):
            if len(group) == 1:
                continue
            bc_dir = f"bc-{count:04d}"
            log.info(f"on barcode {count+1} - {barcode}")
            self.__run_seqkit(fastq1, fastq2, group, bc_dir)
            if self.params["rna-map"]["run"]:
                self.__run_rna_map(bc_dir, all_mhs, rna_map_params)
            count += 1
            if count >= self.params["general"]["max"]:
                log.info("max barcodes reached specified in params")
                break
        if self.params["rna-map"]["run"]:
            json.dump(all_mhs, open("output/BitVector_Files/mutation_histos.json", "w"))


def setup_directories(df, style, dirname="data"):
    os.makedirs(dirname, exist_ok=True)
    bc = 0
    data = []
    for _, g in df.groupby("full_barcode"):
        bc_dir = f"{dirname}/bc-{bc:04d}"
        if style == "individual":
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
            row["sequence"],
            row["structure"],
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
        columns="name,sequence,structure,path,num,barcodes,barcode_bounds,full_barcode,len".split(
            ","
        ),
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
    if len(helices) > 1:
        raise ValueError("only supporting one helix now! Will be changed soon!!")

    def __get_subsection(h1, h2, pos1, pos2):
        h1_new = h1[pos1 : pos2 + 1]
        h2_new = h2[::-1][pos1 : pos2 + 1][::-1]
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


def demultiplex(df, fastq1, fastq2, helices, outdir="data", params=None):
    if params is None:
        params = yaml.safe_load(open(LIB_DIR / "resources/default.yml"))
    # TODO need to validate params
    df_barcodes = find_helix_barcodes(df, helices)
    df_barcodes.to_json(f"barcodes.json", orient="records")
    dmulter = Demultiplexer()
    dmulter.setup(df_barcodes, outdir, params)
    dmulter.run(fastq1, fastq2)
