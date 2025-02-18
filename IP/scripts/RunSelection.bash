#!/bin/bash
era=${1}
channel=${2}

cat > condor/job_${era}_${channel}.submit <<EOF
universe = vanilla
executable = condor/job_${era}_${channel}.sh
log = /dev/null
should_transfer_files = YES
transfer_input_files = condor/job_$_{era}_${channel}.sh
when_to_transfer_output = ON_EXIT
output = out.txt
error = err.txt
notification = Never
initialdir = /afs/cern.ch/work/r/rasp/CMSSW_14_1_0_pre4/src/CPHiggs/IP
+MaxRuntime = 300
+RequestRuntime = 300
MY.WantOS = el9
arguments = 123
queue
EOF


cat > condor/job_${era}_${channel}.sh <<EOF1
cd /afs/cern.ch/work/r/rasp/CMSSW_14_1_0_pre4/src
setenv SCRAM_ARCH el9_amd64_gcc10
cmsenv
cd CPHiggs/IP
./scripts/RunSelection.py --era ${era} --channel ${channel} > condor/job_${era}_${channel}.out
EOF1

chmod u+x condor/job_${era}_${channel}.sh
condor_submit condor/job_${era}_${channel}.submit
