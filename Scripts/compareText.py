"""
compareText.py

This script compares two different directories of text files to see if one differs from the other.
It will throw an error if the datasets differ by more than a small tolerance.

This script has not been tested on TimeSeries-format files.

Usage:
    python3 compareText.py [a_directory] [b_directory]
"""
import sys

TOL = 10e-15

time_format_folder="Final"
vars_folders = ["Conserved", "Auxiliary", "Primitive"]
vars_files = ["cons", "aux", "prims"]
extension = ".dat"

for index in range(len(vars_folders)):
    a_filename = "/".join([sys.argv[1], time_format_folder, vars_folders[index], vars_files[index]])
    b_filename = "/".join([sys.argv[2], time_format_folder, vars_folders[index], vars_files[index]])
    a_filename = a_filename+extension
    b_filename = b_filename+extension
    print("Processing: " + a_filename + ", " + b_filename)

    try:
        with open(a_filename, 'r') as a_dat_file:
            with open(b_filename, 'r') as b_dat_file:
                skip_header = 1
                line_number = 0
                for a_line, b_line in zip(a_dat_file, b_dat_file):
                    if skip_header:
                        skip_header = 0
                        continue
                    a_val = float(a_line)
                    b_val = float(b_line)
                    line_number = line_number + 1
                    if abs(a_val-b_val) > TOL:
                        print("\n\n!! Error in {} (val={}, line={}), {}, (val={})\n\n".format(
                            a_filename, a_val, line_number, b_filename, b_val)
                        )
                        break

    except IOError:
        print("Could not read file:", a_filename, b_filename)


