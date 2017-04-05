"""Script that merges all abundance file produced from StringTie."""

from __future__ import print_function
import pandas as pd
import argparse
import os


def cmdline_parser():
    """Create an argparse instance."""
    parser = argparse.ArgumentParser(description=""""to merge abundance table
        containing FPKM, TPM, from each sample""")
    parser.add_argument("-W", "--WORK_DIR", help="""working directory of the
        project""")
    return parser


def main():
    """The main function."""
    parser = cmdline_parser()
    args = parser.parse_args()

    abun_files = get_file_list(args.WORK_DIR)

    merged_df = merge_df(abun_files)

    out_file = '/'.join([args.WORK_DIR, "all_metrics.csv"])
    merged_df.to_csv(out_file, index=False)


def merge_df(file_list):
    """Merge list of dataframes."""
    df_list = []
    for file in file_list:
        samp_name = file.split("/")[-1].split("_")[0]
        df = pd.read_csv(file, index_col=None, header=0, sep="\t")
        df = df[~df['Gene ID'].str.contains(samp_name)]
        df = df.rename(columns={'FPKM': 'FPKM' + '_' + samp_name,
                                'TPM': 'TPM' + '_' + samp_name,
                                'Coverage': 'Coverage' + '_' + samp_name})
        df_list.append(df)
    all_df = pd.DataFrame(columns=['Gene ID', 'Gene Name', 'Reference',
                                   'Strand', 'Start', 'End', 'Coverage',
                                   'FPKM', 'TPM'])
    for dframe in df_list:
        all_df = pd.merge(all_df, dframe, on=['Gene ID', 'Gene Name',
                                              'Reference', 'Strand',
                                              'Start', 'End', 'Strand'],
                          how='outer')
    all_butdf = all_df.drop(['Coverage', 'FPKM', 'TPM'], axis=1)
    return all_butdf


def get_file_list(workdir):
    """Get list of all tab files."""
    file_list = []
    for root, dirs, files in os.walk(workdir):
        for file in files:
            if file.endswith(".tab"):
                full_file = os.path.abspath(os.path.join(root, file))
                file_list.append(full_file)
    return file_list

if __name__ == '__main__':
    main()
