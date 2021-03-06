"""
text2hdf.py

This script converts text-format output to HDF-format output, for comparison.

This script does not yet work on TimeSeries format outputs.

Usage:
    python3 text2hdf.py [directory] [Optional: output_filename.hdf5]

Authors:
    Sam Mangham (github:smangham)
"""
import h5py
import sys
import numpy as np
from pathlib import Path
from typing import List, Union, Tuple
from h5py import File

final_constants: Path = Path('Final/Constants/constants.dat')
final_domain: Path = Path('Final/Domain/domain.dat')
timeseries_constants: Path = Path('TimeSeries/Constants/constants.dat')

# Which constants belong to which sections?
constants_domain: List[str] = [
    'nx', 'ny', 'nz', 'Nx', 'Ny', 'Nz', 'Ng', 'dx', 'dy', 'dz', 'dt',
    'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax', 'endTime'
]
constants_primitive: List[str] = [
    'Nprims'
]
constants_auxiliary: List[str] = [
    'Naux'
]
constants_conserved: List[str] = [
    'Ncons'
]
constants_root: List[str] = [
    't', 'cfl', 'cp', 'gamma', 'sigma'
]

final_directories: dict = {
    'Auxiliary': {
        'file': 'aux.dat',
        'constants': constants_auxiliary
    },
    'Conserved': {
        'file': 'cons.dat',
        'constants': constants_conserved
    },
    'Primitive': {
        'file': 'prims.dat',
        'constants': constants_primitive
    }
}

# Take input directory from the command line, and check it's one!
if len(sys.argv) < 2:
    raise NotADirectoryError("Please provide an input directory")
directory: Path = Path(sys.argv[1])
if not directory.is_dir():
    raise NotADirectoryError("Please provide an input directory.")

# Create an HDF5 file, and create the two groups within it
hdf5: File = h5py.File(
    sys.argv[2] if len(sys.argv) == 3 else directory.with_suffix('.converted.hdf5'), 'w'
)

# Now we want to start with the Final data; let's look up the constants
with open(directory/final_constants) as f:
    # Names are "Constants = a, b, c" so trim off constants =, split by comma, trim whitespace
    names: List[str] = [s.strip() for s in f.readline().split('=')[1].split(',')]
    # Values are space-separated list e.g. "1 2 3 4 5", but may contain 'inf' or 'nan'...
    values: List[Union[float, int]] = [float(s) if '.' in s or 'n' in s else int(s) for s in f.readline().split()]
    # Create dictionary of constants
    constants: dict = dict(zip(names, values))


# Deciding on the data shape is more complex as it may or may be only 1 or 2-d.
# The capital values ('NZ') include ghost cells. Lower case values ('nz') do not.
if constants['nz']:
    data_shape: Tuple[int, int, int] = (constants['nx'], constants['ny'], constants['nz'])
elif constants['ny']:
    data_shape: Tuple[int, int] = (constants['nx'], constants['ny'])
else:
    data_shape: Tuple[int] = (constants['nx'])

data_length: int = np.product(data_shape)

# Now we write the constants as HDF5 attributes
for constant in constants_root:
    hdf5.attrs.create(constant, constants[constant])

# Domain info is stored in its own format. We want to read that out. However, parallel files do not contain it
domain = hdf5.create_group('Domain')
for constant in constants_domain:
    domain.attrs.create(constant, constants[constant])

try:
    domain.create_dataset(
        name='x',
        data=np.genfromtxt(
            directory/final_domain,
            skip_header=0, max_rows=1
        ).T[constants['Ng']:-constants['Ng']]
    )

    # Y and Z may be single element, in which case there is no domain
    if constants['ny']:
        domain.create_dataset(
            name='y',
            data=np.genfromtxt(
                directory/final_domain,
                skip_header=1, max_rows=1
            ).T[constants['Ng']:-constants['Ng']]
        )

    if constants['nz']:
        domain.create_dataset(
            name='z',
            data=np.genfromtxt(
                directory/final_domain,
                skip_header=2, max_rows=1
            ).T[constants['Ng']:-constants['Ng']]
        )

except Exception as e:
    print("No domain information! OK if this is a parallel text output.")


for subdirectory, details in final_directories.items():
    subgroup = hdf5.create_group(subdirectory)
    for constant in details['constants']:
        subgroup.attrs.create(constant, constants[constant])

    with open(directory/'Final'/subdirectory/details['file']) as f:
        # Names are "Constants = a, b, c" so trim off constants =, split by comma, trim whitespace
        datasets = [s.strip() for s in f.readline().split('=')[1].split(',')]
        for index, dataset in enumerate(datasets):
            # Dataset is N shape x*y*z long, varies fastest in Z, then Y, then X.
            data = np.genfromtxt(
                directory/'Final'/subdirectory/details['file'],
                skip_header=1 + index*data_length, max_rows=data_length
            ).reshape(
                data_shape
            )
            subgroup.create_dataset(name=dataset, data=data)
