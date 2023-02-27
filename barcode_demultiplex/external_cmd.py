import shutil
import subprocess
from typing import Optional, List
from dataclasses import dataclass


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
    for seqs, bounds in zip(barcode_seqs, barcode_bounds):
        max_len = opts.seq_len
        if max_len > opts.read_len:
            max_len = opts.read_len
        for seq, bound in zip(seqs, bounds):
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
                cmds.append(
                    f'seqkit grep -s -p "{seq}" -m 1 -P -R {b1}:{b2} '
                )
    cmd = "| ".join(cmds) + f"-o {out_fastq_file}"
    return run_command(cmd)
