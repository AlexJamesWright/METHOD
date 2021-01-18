## Tests Instructions

These are instructions to run CPU unit tests as a batch job on Iridis 5

## Setting up python env

In the root METHOD folder, create a python venv using

```
module purge
module load gcc/6.4.0
module load python/3.6.4
module load hdf5/1.10.2/gcc/parallel
python3 -m venv venv
source venv/bin/activate
```

Then install python modules using

```
python -m pip install -r Scripts/IridisEnv/requirements.txt
```

## Runing unit tests as a batch job

From `Tests/CPU` run `sbatch ../../Scripts/IridisEnv/tests_job.sh`

This will run all CPU tests including tests of the hdf5 serial and parallel writers




