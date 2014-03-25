#!/bin/bash

if [ $# -ne 1 ]
then
	echo "Usage: [run dir]"
	exit 1
fi

# Nestly binaries
RUN="$HOME/python/bin/nestrun"

# Concurrency of jobs
OPTIONS="-j 20"

# Make WGS
$RUN $OPTIONS -d $1 --template "$PWD/bin/pbsrun_ART.sh {seed} {wgs_insert_length} {wgs_insert_sd} {wgs_read_length} {wgs_xfold} {community}/{refseq} {wgs_base}"

# Make assemblies
$RUN $OPTIONS -d $1 --template "$PWD/bin/a5submit.sh -mw {wgs_base}1.fq {wgs_base}2.fq {wgs_base}"
