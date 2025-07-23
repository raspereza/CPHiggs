#!/bin/bash
era=$1
channel=$2
sample=$3
period=$4

folder=/afs/cern.ch/work/r/rasp/CPHiggs/Analysis

cat > ${folder}/condor/job_${era}_${channel}_${sample}.submit <<EOF
universe = vanilla
executable = /afs/cern.ch/work/r/rasp/CPHiggs/Analysis/condor/job_${era}_${channel}_${sample}.sh
output = /afs/cern.ch/work/r/rasp/CPHiggs/Analysis/condor/job_${era}_${channel}_${sample}.out
error = /afs/cern.ch/work/r/rasp/CPHiggs/Analysis/condor/job_${era}_${channel}_${sample}.err
log = /afs/cern.ch/work/r/rasp/CPHiggs/Analysis/condor/job_${era}_${channel}_${sample}.log
notification = Never
initialdir = /afs/cern.ch/work/r/rasp/CMSSW_14_1_0_pre4/src/CPHiggs/Analysis
+MaxRuntime = ${period}
+RequestRuntime = ${period}
MY.WantOS = el9
queue
EOF


cat > /afs/cern.ch/work/r/rasp/CPHiggs/Analysis/condor/job_${era}_${channel}_${sample}.sh <<EOF1
#!/bin/bash
cd /afs/cern.ch/work/r/rasp/CMSSW_14_1_0_pre4/src
echo $PWD
export SCRAM_ARCH=el9_amd64_gcc10
cmsenv
cd CPHiggs/Analysis
./scripts/RunSelection.py --era ${era} --channel ${channel} --sample ${sample} --applyIPSigLep1Cut --applyIPSigLep2Cut
EOF1

chmod u+x /afs/cern.ch/work/r/rasp/CPHiggs/Analysis/condor/job_${era}_${channel}_${sample}.sh
condor_submit /afs/cern.ch/work/r/rasp/CPHiggs/Analysis/condor/job_${era}_${channel}_${sample}.submit
