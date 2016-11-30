#!/usr/bin/env python

"""
Convert eukaryotic gff3 files to gff files and gtf files.

required for downstream analsysis, will use gffread

"""
# --- standard library imports
import argparse # noqa
import os # noqa
# --- third-party imports

# --- project specific imports


def cmdline_parser():
    """
    Create an argparse instance.

    Combination of different options for this script.
    """
    parser = argparse.ArgumentParser(description="""converts eukaryotic gff3
        files to gff and gtf files""")
    parser.add_argument("-f", "--FASTA", help="""fasta file corresponding to
        the gff gile""")
    parser.add_argument("-g", "--GFF", help="""input gff file""")
    parser.add_argument("-o", "--OUT_TABLE", help="""table from the GOTTCHA2
        sequence header""")
    return parser


def main():
    """
    Main function.

    All functions are called here.
    """
    parser = cmdline_parser()
    args = parser.parse_args()


def grep_convert(fasta, out_table):
    """
    Convert fasta header to table.

    from fasta file generate a table file from sequence header
    """

if __name__ == '__main__':
    main()
