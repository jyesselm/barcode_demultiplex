import click
import pandas as pd

from barcode_demultiplex.demultiplex import *
from barcode_demultiplex.logger import get_logger, setup_applevel_logger

log = get_logger("CLI")


# cli #########################################################################


@click.command(
    help="a lightweight program to demultiplex rna sequencing data from internal"
    " barcodes sequestered in helices. Required arguments are a csv file thats"
    " contains the name, sequence and structure of each construct "
)
@click.option(
    "-csv",
    "--rna-csv",
    required=True,
    help="a csv containing the name, sequence and structure of all RNAs that "
    "have reads",
)
@click.option(
    "-fq1", "--fastq1", required=True, help="the path to the forward fastq file"
)
@click.option(
    "-fq2",
    "--fastq2",
    required=False,
    default=None,
    help="the path to the reverse fastq file",
)
@click.option(
    "-rd", "--run-dreem", is_flag=True, help="should we run dreem after?"
)
@click.option(
    "--include-all-dreem-outputs",
    is_flag=True,
    help="usually want only summary outputs but some cases want everythign from dreem",
)
@click.option(
    "-dp",
    "--data-path",
    default="data",
    help="the location of where to store all the demultiplexed data default='data'",
)
@click.option(
    "-m",
    "--max",
    "max_barcodes",
    default=99999,
    help="used for debugging purposes to see if a few barcodes work",
)
@click.option(
    "-helix",
    "--helix",
    "helices",
    type=click.Tuple([int, int, int]),
    required=True,
    multiple=True,
    help="which helices should we use as barcodes? In the format helix number as it "
    "appears 5' to 3' and what are the bounds of barcode in the helix. If the "
    "the barcode is from the beginning of the helix to the 8th position this can "
    "be specified as -helix 1 0 8 assuming this is the second helix as 0 is the first",
)
def cli(
    rna_csv,
    fastq1,
    fastq2,
    helices,
    run_dreem,
    data_path,
    max_barcodes,
    include_all_dreem_outputs,
):
    setup_applevel_logger()
    df = pd.read_csv(rna_csv)
    log.info(f"{rna_csv} contains {len(df)} unique sequences")
    required_cols = "name,sequence,structure".split(",")
    for c in required_cols:
        if c not in df:
            raise ValueError(
                f"{c} is a required column for the input csv file!"
            )
    df = find_helix_barcodes(df, helices)
    df_sum = setup_directories(df, data_path)
    log.info(f"{len(df_sum)} unique barcodes found!")
    demult = Seqkitdemultiplexer()
    demult.max = max_barcodes
    if fastq2 is None:
        demult.check_rev = False
    demult.run(df_sum, fastq1, fastq2)

    if run_dreem:
        run_dreem_prog(df_sum, max_barcodes, include_all_dreem_outputs)


if __name__ == "__main__":
    cli()
