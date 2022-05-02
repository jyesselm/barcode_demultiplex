from barcode_demultiplex.run import *
from click.testing import CliRunner
from pathlib import Path
import numpy as np
import os
import shutil

PATH = str(Path(__file__).parents[1])

np.warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

def test_run():
    example_path = PATH + "/examples/ires_minittrs/"
    args = [
        "-csv",
        f"{example_path}/test.csv",
        "-fq1",
        f"{example_path}/test_mate1.fastq",
        "-fq2",
        f"{example_path}/test_mate2.fastq",
        "-helix",
        "1",
        "0",
        "8",
        "--max",
        1,
    ]
    runner = CliRunner()
    result = runner.invoke(cli, args, prog_name="barcode_demultiplex")
    assert result.exit_code == 0
    assert os.path.isfile("data/bc-0000/test_mate1.fastq")
    assert os.path.isfile("data/bc-0000/test_mate2.fastq")
    shutil.rmtree("data")

def test_run_only_fwd():
    example_path = PATH + "/examples/ires_minittrs/"
    args = [
        "-csv",
        f"{example_path}/test.csv",
        "-fq1",
        f"{example_path}/test_mate1.fastq",
        "-helix",
        "1",
        "0",
        "8",
        "--max",
        1,
    ]
    runner = CliRunner()
    result = runner.invoke(cli, args, prog_name="barcode_demultiplex")
    assert result.exit_code == 0
    assert os.path.isfile("data/bc-0000/test_mate1.fastq")
    shutil.rmtree("data")
