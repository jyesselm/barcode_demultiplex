from barcode_demultiplex.external_cmd import (
    run_seqkit_grep_fwd,
    SeqkitGrepOpts,
)


def test_run_seqkit_grep_fwd():
    opts = SeqkitGrepOpts(
        mut_num=0,
        read_len=150,
        seq_len=169,
        buffer=5,
    )
