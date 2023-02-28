import shutil
import subprocess
from typing import Optional
from dataclasses import dataclass

from seq_tools.sequence import get_reverse_complement


@dataclass(frozen=True, order=True)
class ProgOutput:
    """
    Class to store the output of an external program
    """

    output: Optional[str]
    error: Optional[str]


@dataclass(frozen=True, order=True)
class SeqkitGrepOpts:
    mut_num: int
    read_len: int
    seq_len: int
    buffer: int = 5


def does_program_exist(prog_name: str) -> bool:
    """
    Check if a program exists
    :param prog_name: name of the program
    """
    if shutil.which(prog_name) is None:
        return False
    else:
        return True


def run_command(cmd: str) -> ProgOutput:
    """
    Run a command and return the output
    :param cmd: command to run
    """
    output, error_msg = None, None
    try:
        output = subprocess.check_output(
            cmd, shell=True, stderr=subprocess.STDOUT
        ).decode("utf8")
    except subprocess.CalledProcessError as exc:
        error_msg = exc.output.decode("utf8")
    return ProgOutput(output, error_msg)


def run_seqkit_grep_fwd(
    barcode_seqs,
    barcode_bounds,
    in_fastq_file,
    out_fastq_file,
    opts: SeqkitGrepOpts,
) -> ProgOutput:
    count = 0
    cmds = []
    max_len = opts.seq_len
    if max_len > opts.read_len:
        max_len = opts.read_len
    for seq, bound in zip(barcode_seqs, barcode_bounds):
        b1, b2 = bound[0] - opts.buffer, bound[1] + opts.buffer
        if b2 >= max_len:
            continue
        count += 1
        if count == 1:
            cmds.append(
                f'seqkit grep -s -p "{seq}" -m 1 -P -R {b1}:{b2} '
                f"{in_fastq_file} "
            )
        else:
            cmds.append(f'seqkit grep -s -p "{seq}" -m 1 -P -R {b1}:{b2} ')
    cmd = "| ".join(cmds) + f" -o {out_fastq_file}"
    return run_command(cmd)


def run_seqkit_grep_rev(
    barcode_seqs,
    barcode_bounds,
    in_fastq_file,
    out_fastq_file,
    opts: SeqkitGrepOpts,
) -> ProgOutput:
    # TODO hardcoded for now but you lost 12 nt from the RTB barcodes on R1
    rtb_length = 12
    count = 0
    cmds = []
    max_len = opts.seq_len
    if max_len > opts.read_len - rtb_length:
        max_len = opts.read_len - rtb_length
    for seq, bound in zip(barcode_seqs, barcode_bounds):
        # need rev complement
        seq = get_reverse_complement(seq)
        b1, b2 = bound[0], bound[1]
        # remap direction
        b3 = opts.seq_len - b1 + opts.buffer
        b4 = opts.seq_len - b2 - opts.buffer
        if b3 >= max_len:
            continue
        count += 1
        if count == 1:
            cmds.append(
                f'seqkit grep -s -p "{seq}" -m 1 -P -R {b4}:{b3} '
                f"{in_fastq_file} "
            )
        else:
            cmds.append(f'seqkit grep -s -p "{seq}" -m 1 -P -R {b4}:{b3} ')
    cmd = "| ".join(cmds) + f" -o {out_fastq_file}"
    return run_command(cmd)


def run_seqkit_common(in_fastq_file_1, in_fastq_file_2, out_fastq_file):
    """
    Finds reads in common betwen in_fastq_file_1 and in_fastq_file_2 and
    outputs them to out_fastq_file
    :param in_fastq_file_1: first fastq file. These are the reads that will be
    moved into out_fastq_file
    :param in_fastq_file_2: second fastq file.
    :param out_fastq_file: Outputed fastq files that combine the read ids
    contained between the two files
    :return: the stdout/stderr of runs
    """
    cmd = (
        f"seqkit common {in_fastq_file_1} {in_fastq_file_2} -o {out_fastq_file}"
    )
    run_command(cmd)
