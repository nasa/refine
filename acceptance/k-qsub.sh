#!/usr/bin/env bash

set -e # exit on first error
set -u # Treat unset variables as error
set -x # echo

uname -n

pwd # echo path

testname=$1

# Get the best queue to submit to (removed disabled K3a-debug)
queues="K4-debug K3-debug K2-debug K2-fun3d"
for queue in $queues; do

    if [[ "$queue" == "K4-"* ]]; then
	nprocs=40
	walltime=50
    elif [[ "$queue" == "K3-"* ]]; then
	nprocs=16
	walltime=50
    elif [[ "$queue" == "K3a-"* ]]; then
	nprocs=16
	walltime=50
    elif [[ "$queue" == "K2-"* ]]; then
	nprocs=12
	walltime=50
    else
	echo "unknown queue requested: $queue"
	exit 1
    fi

    output=`pwd`/../log-${testname}-status.txt

    cat << EOF > ${testname}.pbs
#PBS -S /bin/bash -q ${queue} -l select=1:ncpus=${nprocs}:mpiprocs=${nprocs} -l walltime=00:${walltime}:00
#PBS -N ${testname}
#PBS -o ${testname}
#PBS -W umask=027
#PBS -m n

set -o xtrace   # show commands invoked
set -o errexit  # exit on first error
set -o pipefail # ensure exit with pipe

uname -mrn

cd \$PBS_O_WORKDIR

hostname 2>&1 | tee ${output}
pwd 2>&1 | tee ${output}

(./acceptance/${testname}.sh) 2>&1 | tee ${output}

EOF

    printf "Trying K queue: $queue"

    rm -rf ${output}
    (qsub -V -Wblock=true ${testname}.pbs > pbs_id 2>&1) &
    pid=$!

    sleep 20

    job_running="R"
    job_finished="F"
    job_status=$(qstat -x `cat pbs_id` | tail -n1 | sed -r 's/\s+/ /g' | cut -d' ' -f8)
    if [[ "${job_status}" == "${job_running}" ]]; then
	echo "Job Running!"
	break
    elif [[ "${job_status}" == "${job_finished}" ]]; then
	echo "Job Finished!"
	break
    else
	echo "try next queue! Status >${job_status}<"
	qdel -Wforce `cat pbs_id`
    fi
done

echo "===================================================================="
echo "     ------>  Running ${testname} on K queue: ${queue}"
echo "===================================================================="

# This is where the winning PBS job starts,
# due to the $BUILDLOG file being created
set -x
touch $output
tail -n+0 -f --pid $pid $output

# Capture test-suite error code, turn of exit-on-error
set +e
wait ${pid}; exit_code=$?

exit ${exit_code}

