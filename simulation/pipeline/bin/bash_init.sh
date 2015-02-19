set -o pipefail
set -e

#prevents uninitialised variables but this causes an issue with PBS
#we would need to set PBS_ENVIRONMENT to something in each script.
#set -u

function rollback_rm_file() {
	echo "Rollback from premature exit" > /dev/stderr
	echo "Removing output target $1" > /dev/stderr
	rm $1
}

function rollback_rm_dir() {
	echo "Rollback from premature exit" > /dev/stderr
	echo "Removing directory $1"
	rm -r $1
}

function rollback_rm_files() {
	echo "Rollback from premature exit" > /dev/stderr
	echo "Removing list of files ${1[@]}"
	for fn in "${1[@]}"
	do
		echo "Removing $fn"
		rm $fn
	done
}
