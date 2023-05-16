#!/bin/bash
PARTITION="markov"
# make_ext_path outputs the new path required for external partions.
function make_ext_path() {
	local path=$1
	[[ $path == "/datainbackup/"* ]] \
	&& echo $path | \
		awk '{sub("/datainbackup/", "/bayes_datainbackup/"); print;}'
	[[ $path == "/data/"* ]] \
	&& echo $path | \
		awk '{sub("/data/", "/bayes_data/"); print;}'
}
# Determine if partition is on another machine.
[[ $PARTITION == *"shannon" ]] || [[ $PARTITION == *"markov" ]] \
	&& EXTERNAL_PARTITION=1 \
	|| EXTERNAL_PARTITION=0

# Determine SLURM output directory depending on partition.
OUTDIR="$(pwd -P)"
[[ $EXTERNAL_PARTITION -eq 1 ]] && OUTDIR=$(make_ext_path $OUTDIR)

# Modify $SCRIPT location if on an external partition.
# We can do below because we know user is in its home directory and we have 
# the real path of $SCRIPT (i.e. not the symlink).
[[ $EXTERNAL_PARTITION -eq 1 ]] && SCRIPT=$(make_ext_path $SCRIPT)

# Notes about sbatch script:
#
# If we use a partition on another machine we only define SLURM_* variables
#	from the user environment (using --export=NONE); i.e. no user-defined
#	variables will be transferred. This is so that e.g. TMPDIR will
#	point to whatever the other machine point TMPDIR to, instead of 
#	/data/tmp, because that directory doesn't exist on other machines.
#	(NFS mounts bayes:/data --> other_machine:/bayes_data.)

SBATCH_SCRIPT="#!/bin/bash
#SBATCH --partition=$PARTITION
#SBATCH --cpus-per-task=$CPUS_PER_TASK
#SBATCH --output=$OUTDIR/slurm-%j.out
#SBATCH --error=$OUTDIR/slurm-%j.error
#SBATCH --chdir=$OUTDIR # Working directory
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1G
$([[ $EXTERNAL_PARTITION -eq 1 ]] && echo '#SBATCH --export=ALL,TEMP=/scratch,TMP=/scratch,TMPDIR=/scratch')

/data/MS-2023/atjohan/bin/apptainer run --bind ./csv_results:/thesis/csv_results/,./hdf5_results:/thesis/hdf5_results thesis.sif
"