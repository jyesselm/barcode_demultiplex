import os
from pathlib import Path

from barcode_demultiplex.external_cmd import (
    run_seqkit_grep_fwd,
    run_seqkit_grep_rev,
    run_seqkit_common,
    SeqkitGrepOpts,
)

TEST_DIR = Path(os.path.dirname(os.path.realpath(__file__)))


def test_run_seqkit_grep_fwd():
    opts = SeqkitGrepOpts(
        mut_num=0,
        read_len=151,
        seq_len=169,
        buffer=5,
    )
    barcode_seqs = ["GTCGTCGA"]
    barcode_bounds = [[47, 54]]
    in_fastq_file = TEST_DIR / "resources/test_cases/C0098/R2.sub.fastq.gz"
    out_fastq_file = "test.fastq.gz"
    run_seqkit_grep_fwd(
        barcode_seqs, barcode_bounds, in_fastq_file, out_fastq_file, opts
    )
    assert Path(out_fastq_file).exists()
    os.remove(out_fastq_file)


def test_run_seqkit_grep_fwd_2():
    opts = SeqkitGrepOpts(
        mut_num=1,
        read_len=151,
        seq_len=169,
        buffer=3,
    )
    barcode_seqs = ["GTCGTCGA", "TCGATGAC"]
    barcode_bounds = [[47, 54], [139, 146]]
    in_fastq_file = TEST_DIR / "resources/test_cases/C0098/R2.sub.fastq.gz"
    out_fastq_file = "test.fastq.gz"
    run_seqkit_grep_fwd(
        barcode_seqs, barcode_bounds, in_fastq_file, out_fastq_file, opts
    )
    assert Path(out_fastq_file).exists()
    os.remove(out_fastq_file)


def test_run_seqkit_grep_rev():
    opts = SeqkitGrepOpts(
        mut_num=0,
        read_len=151,
        seq_len=169,
        buffer=5,
    )
    barcode_seqs = ["TCGATGAC"]
    barcode_bounds = [[139, 146]]
    in_fastq_file = TEST_DIR / "resources/test_cases/C0098/R1.sub.fastq.gz"
    out_fastq_file = "test.fastq.gz"
    run_seqkit_grep_rev(
        barcode_seqs, barcode_bounds, in_fastq_file, out_fastq_file, opts
    )
    assert Path(out_fastq_file).exists()
    os.remove(out_fastq_file)


def test_run_seqkit_grep_rev_2():
    opts = SeqkitGrepOpts(
        mut_num=1,
        read_len=151,
        seq_len=169,
        buffer=3,
    )
    barcode_seqs = ["GTCGTCGA", "TCGATGAC"]
    barcode_bounds = [[47, 54], [139, 146]]
    in_fastq_file = TEST_DIR / "resources/test_cases/C0098/R1.sub.fastq.gz"
    out_fastq_file = "test.fastq.gz"
    run_seqkit_grep_fwd(
        barcode_seqs, barcode_bounds, in_fastq_file, out_fastq_file, opts
    )
    assert Path(out_fastq_file).exists()
    os.remove(out_fastq_file)


def test_run_seqkit_common():
    """
    test run_seqkit_common command
    :return:
    """
    in_fastq_file_1 = TEST_DIR / "resources/test_cases/C0098/R1.sub.fastq.gz"
    in_fastq_file_2 = TEST_DIR / "resources/test_cases/C0098/R2.sub.fastq.gz"
    out_fastq_file = "test.fastq.gz"
    run_seqkit_common(in_fastq_file_1, in_fastq_file_2, out_fastq_file)
    assert Path(out_fastq_file).exists()
    os.remove(out_fastq_file)

