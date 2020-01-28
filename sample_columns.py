# Program to sample cells from paired tumor and non-tumor samples; in the input files,
# rows are genes and columns are cells, so sampling cells in effect means sampling columns.
# Author: Pattara Sukprasert
import argparse
import random
import os
import errno
from typing import Optional

def parse() -> argparse.Namespace:
    """ Parse arguments into argparse.Namespace to be used later.
    
    Returns:
        A namespace containing all parameters.
    """
    parser = argparse.ArgumentParser(description='''
    Given a dataset generate smaller datasets with specific number of cells randomly.
    This tool can be used to generate many replicates in one run when the replication flag is
    specified. Also, the output directory can be different from the directory where the tool is run.
    A key assumption is that tumor files and non-tumor files are paired by having similar names,
    if there are any non-tumor files; one can also use the program with only tumor files.
    ''')

    # Note: When nargs='+' that means the argument is a actually a non-empty list, which may contain as few as one item
    # When nargs='*' that means the argument is a actually a list, which may be empty
    # In our usages, the values for args.tumor_files and args.non_tumor_files may have the * character in them
    # When the * character is used, the string is expanaded to a list by bash; the user must choose the strings carefully
    # to ensure that the pairs of files that are supposed to be matched occur in the same positions (ranks or indices) in
    # the two expanded lists
    parser.add_argument('--tumor_files', metavar='FILE', type=str, nargs='+',
                        help='Tumor file names in space-separated-value format.')
    parser.add_argument('--non_tumor_files', metavar='MEAN_FILE', type=str, nargs='*',
                        help='Non-tumor file names used for calculating average in'
                        + 'space-separated-value format.')

    # This flag allows some flexibility in the data format.
    parser.add_argument('--data_column', default=1, type=int,
                        help='The first column that store data. [default=1]')

    
    parser.add_argument('--num_cells', type=int, nargs='+',
                        help='Number of cells to be sampled.')

    parser.add_argument('--output_dir', type=str,
                        help='Output directory of the subset of cells.')

    parser.add_argument('--replication', type=int, default=1,
                        help='Number of replication for each dataset.')
    
    # This flag allows dataset replication across machines.
    parser.add_argument('--seed', type=int, default=-1,
                        help=(
                        'Seed to be used for pseudo random number '
                        ' generator. If -1 (default), then do not specify any seed'
                        ))


    return parser.parse_args()

def sample(args: argparse.Namespace, tumor_file: str, non_tumor_file: Optional[str] = None) -> None:
    """ sample is the central procedure for sampling from one pair of files.

    Args:
        args: The namespace containing command line arguments.
        tumor_file: A path containing tumor cells of an individual.
        non_tumor_file: A path containing non-tumor cells of the same individual.
        This file, if specified, should have the same named rows )representing genes)
        as the file containing tumor cells. The named columns (representing cells) will differ
        in the two files, but should have values in the same format.
    """

    tumor_data = open(tumor_file, 'r')
    tumor_lines = list(line.rstrip().split('\t') for line in tumor_data)
    # Assign to cell_set the set of columns to the left of the data_column.
    cell_set = set(tumor_lines[0][args.data_column:])
    if non_tumor_file is not None:
        non_tumor_data = open(non_tumor_file, 'r')
        non_tumor_lines = list(line.rstrip().split('\t') for line in non_tumor_data)
        non_tumor_cell_set = set(non_tumor_lines[0][args.data_column:])
        # Re-assign to cell_set the union of the tumor cells and the non-tumor cells
        cell_set = cell_set.union(non_tumor_cell_set)

    cell_set = list(cell_set)
    cell_set.sort() # sort to make sure that the pseudo rng seed works.

    for numcell in args.num_cells:
        for i in range(args.replication):
            # Decide on the number of cells to sample, which cannot be more than the number available
            actual_numcell = min(len(cell_set), numcell)
            selected_cell_set = set(random.sample(list(cell_set), actual_numcell))

            with open(args.output_dir+'/num_cell_{}_replication_{}_{}_actual_{}'.format(
                numcell, i+1, tumor_file[tumor_file.rfind('/') + 1:], actual_numcell), 'w') as file_to_write:

                for line in tumor_lines:
                    # Select non data-point column.
                    left = line[:args.data_column]

                    # To the right of the data column, include only those items in the selected cell set
                    right = list(line[k] for k in range(args.data_column, len(line))
                                 if tumor_lines[0][k] in selected_cell_set)
                    file_to_write.write('\t'.join(left+right) + '\n')


            if non_tumor_file is not None:
                with open(args.output_dir+'/num_cell_{}_replication_{}_{}'.format(
                    numcell, i+1, non_tumor_file[non_tumor_file.rfind('/') + 1:]), 'w') \
                as file_to_write:
                    for line in non_tumor_lines:
                        left = line[:args.data_column]
                        right = list(line[k] for k in range(args.data_column, len(line))
                                     if non_tumor_lines[0][k] in selected_cell_set)
                        file_to_write.write('\t'.join(left+right) + '\n')

def main():
    args = parse()

    if args.seed != -1:
        random.seed(args.seed)

    if not os.path.exists(os.path.dirname(args.output_dir + '/')):
        try:
            os.makedirs(os.path.dirname(args.output_dir + '/'))
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
    # Run sampling for each (pair of) file independently.
    for idx, _ in enumerate(args.tumor_files):
        if args.non_tumor_files:
            sample(args, args.tumor_files[idx], args.non_tumor_files[idx])
        else:
            sample(args, args.tumor_files[idx])

# Call main().
if __name__ == '__main__':
    main() 
