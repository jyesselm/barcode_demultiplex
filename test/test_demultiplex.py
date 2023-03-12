import os
from pathlib import Path
import shutil
import pandas as pd
from barcode_demultiplex.demultiplex import (
    find_helix_barcodes,
    get_read_length,
    Demultiplexer,
)

TEST_DIR = Path(os.path.dirname(os.path.realpath(__file__)))


class TestResources:
    @staticmethod
    def get_test_params():
        params = {
            "fwd": {
                "check": True,
                "max_nuc": 999,
                "min_nuc": 0,
                "checks": 999,
                "mut_num": 1,
                "buffer": 5,
            },
            "rev": {
                "check": True,
                "max_nuc": 999,
                "min_nuc": 0,
                "checks": 999,
                "mut_num": 1,
                "buffer": 5,
            },
        }
        return params


def test_find_helix_barcodes():
    path = TEST_DIR / "resources/test_cases/C0098"
    df = pd.read_csv(path / "C0098.csv")
    df = df[["name", "sequence", "structure"]]
    helices = [[2, 0, 8]]
    df_barcodes = find_helix_barcodes(df, helices)
    row = df_barcodes.iloc[0]
    assert row["barcodes"][0][0] == "GTCGTCGA"
    assert row["barcodes"][0][1] == "TCGATGAC"
    assert row["barcode_bounds"][0][0] == [47, 54]
    assert row["barcode_bounds"][0][1] == [139, 146]


def test_get_read_length():
    path = TEST_DIR / "resources/test_cases/C0098"
    R1_path = path / "R1.sub.fastq.gz"
    R2_path = path / "R2.sub.fastq.gz"
    fwd_read_len = get_read_length(R1_path)
    rev_read_len = get_read_length(R2_path)
    assert fwd_read_len == 139
    assert rev_read_len == 151


def test_demultiplexer():
    path = TEST_DIR / "resources/test_cases/C0098"
    R1_path = path / "R1.sub.fastq.gz"
    R2_path = path / "R2.sub.fastq.gz"
    params = TestResources.get_test_params()
    df_barcodes = pd.read_json(path / "C0098_barcodes.json")
    dmulter = Demultiplexer()
    dmulter.setup(df_barcodes, "data", params)
    dmulter.run(R2_path, R1_path)
    assert Path("data").is_dir()
    barcode_dirs = list(Path("data").iterdir())
    assert len(barcode_dirs) == 24
    # shutil.rmtree("data")
