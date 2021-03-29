"""
compareHDF5.py

This script compares two different HDF5 files to see if one differs from the other.
It will throw an error if they have different groups, attributes or datasets,
or warn if the attributes or datasets have different values (beyond a small tolerance).

This script has not been tested on TimeSeries-format files, but should work.

Usage:
    python3 compareHDF5.py [a_file.hdf5] [b_file.hdf5]

Authors:
    Sam Mangham (github:smangham)
"""
import sys
import numpy as np
import h5py


def compare(file1, file2):
    # These attributes can be absent or different- it doesn't matter!
    whitelist_attributes: list = [
        'Nx', 'Ny', 'Nz'
    ]
    # These datasets can be absent or different- it doesn't matter!
    whitelist_datasets: list = [
        'x', 'y', 'z'
    ]
    
    warnings_found = False
    
    print(
        "Ignoring differences in:\n"
        "- Attributes: "+', '.join(whitelist_attributes)+"\n"
        "- Datasets: "+', '.join(whitelist_datasets)
    )
    
    try:
        a_file = h5py.File(file1, 'r')
    except IOError:
        print("Error: Could not read file:", file1)
        exit(1)
    
    try:
        b_file = h5py.File(file2, 'r')
    except IOError:
        print("Error: Could not read file:", file2)
        exit(1)
    
    # First, we look at the top-level groups: Are they the same?
    group_difference: set = set(a_file.keys()) - set(b_file.keys())
    if len(group_difference):
        print(
            "Error: Root groups differ!\n"
            " - "+file1+": "+', '.join(a_file.keys())+"\n"
            " - "+file2+": "+', '.join(b_file.keys())
        )
        exit(1)
    
    
    # Then, we compare the attributes of the two files: Do they have the same attributes, and do those
    # attributes have the same values?
    attribute_difference: set = set(a_file.attrs.keys()) - set(b_file.attrs.keys()) - set(whitelist_attributes)
    if len(attribute_difference):
        print(
            "Error: Root attributes differ!\n"
            " - "+file1+": "+a_file.attrs.keys()+"\n"
            " - "+file2+": "+b_file.attrs.keys()
        )
        exit(1)
    
    for attribute_name, a_attribute in a_file.attrs.items():
        if attribute_name not in whitelist_attributes:
            b_attribute = b_file.attrs[attribute_name]
            if not np.allclose(a_attribute, b_attribute):
                warnings_found = True
                print(
                    "Warning: root attribute '"+attribute_name+"' values differ!\n",
                    " - "+file1+": "+a_attribute+"\n"
                    " - "+file2+": "+b_attribute
                )
    
    # For each group, compare the attributes and dataset values
    for group_name, a_group in a_file.items():
        b_group = b_file[group_name]
        attribute_difference: set = set(a_group.attrs.keys()) - set(b_group.attrs.keys()) - set(whitelist_attributes)
        if len(attribute_difference):
            print(
                "Error: "+group_name+" attributes differ!\n"
                " - "+file1+": "+', '.join(a_group.attrs.keys())+"\n"
                " - "+file2+": "+', '.join(b_group.attrs.keys())
            )
            exit(1)
    
        for attribute_name, a_attribute in a_group.attrs.items():
            if attribute_name not in whitelist_attributes:
                b_attribute = b_group.attrs[attribute_name]
                if a_attribute.dtype.char == 'S':
                    if not a_attribute == b_attribute:
                        warnings_found = True
                        print(
                            "Warning: "+group_name+" attribute '"+attribute_name+"' values differ!\n"
                            " - "+file1+": "+a_attribute+"\n"
                            " - "+file2+": "+b_attribute
                        )
                elif not np.allclose(a_attribute, b_attribute):
                    warnings_found = True
                    print(
                        "Warning: "+group_name+" attribute '"+attribute_name+"' values differ!\n"
                        " - "+file1+": "+a_attribute+"\n"
                        " - "+file2+": "+b_attribute
                    )
    
        dataset_difference: set = set(a_group.keys()) - set(b_group.keys()) - set(whitelist_datasets)
        if len(attribute_difference):
            print(
                "Error: "+group_name+" datasets differ!\n"
                " - "+file1+": "+', '.join(a_group.keys())+"\n"
                " - "+file2+": "+', '.join(b_group.keys())
            )
            exit(1)
    
        for dataset_name, a_dataset in a_group.items():
            if dataset_name not in whitelist_datasets:
                b_dataset: h5py.Dataset = b_file[group_name][dataset_name]
                if not a_dataset.shape == b_dataset.shape:
                    warnings_found = True
                    print(
                        "Error: "+group_name+" datasets "+dataset_name+" shapes differ!\n"
                        " - "+file1+": "+str(a_dataset.shape)+"\n"
                        " - "+file2+": "+str(b_dataset.shape)
                    )
    
                elif not np.allclose(a_dataset, b_dataset, atol=1e-16):
                    warnings_found = True
                    print(
                        "Warning: "+group_name+" dataset '"+dataset_name+"' values differ!"
                    )
    
    if not warnings_found:
        print("Files "+file1+" and "+file2+" are the same to within tolerances")
        return 1


if __name__ == "__main__":
    import sys
    compare(sys.argv[1], sys.argv[2])
