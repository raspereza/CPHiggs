#!/bin/bash
era=$1
channel=$2
period=$3

cat > condor/job_${era}_${channel}.submit <<EOF
universe = vanilla
executable = condor/job_${era}_${channel}.sh
output = condor/job_${era}_${channel}.out
error = condor/job_${era}_${channel}.err
log = condor/job_${era}_${channel}.log
notification = Never
initialdir = /afs/cern.ch/work/r/rasp/CMSSW_14_1_0_pre4/src/CPHiggs/IP
+MaxRuntime = ${period}
+RequestRuntime = ${period}
MY.WantOS = el9
queue
EOF


cat > condor/job_${era}_${channel}.sh <<EOF1
#!/bin/bash
cd /afs/cern.ch/work/r/rasp/CMSSW_14_1_0_pre4/src
echo $PWD
export SCRAM_ARCH=el9_amd64_gcc10
cmsenv
cd CPHiggs/IP
./scripts/RunSelection.py --era ${era} --channel ${channel} --applyIPSigPromptLepSF 
EOF1

chmod u+x condor/job_${era}_${channel}.sh
condor_submit condor/job_${era}_${channel}.submit
