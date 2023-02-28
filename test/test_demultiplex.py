import os
from pathlib import Path
import pandas as pd
from barcode_demultiplex.demultiplex import find_helix_barcodes, Demultiplexer

TEST_DIR = Path(os.path.dirname(os.path.realpath(__file__)))


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


def test_demultiplexer():
    path = TEST_DIR / "resources/test_cases/C0098"
    R1_path = path / "R1.sub.fastq.gz"
    R2_path = path / "R2.sub.fastq.gz"
    df_barcodes = pd.read_json(path / "C0098_barcodes.json")
    dmulter = Demultiplexer()
    dmulter.setup(df_barcodes, "data")
    dmulter.run(R2_path, R1_path)




