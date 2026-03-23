#!/bin/bash

n=$#
if [[ $n -ne 6 ]]; then
    echo usage : RunSelection.bash [ERA] [CHANNEL] [SAMPLE] [ANALYSIS_TYPE] [APPLY_IPCUT] [APPLY_FF]
    echo ERA = [Run3_2022preEE, Run3_2022postEE, Run3_2022preBPix, Run3_2022postBPix]
    echo CHANNEL = [mt,et,mm,ee]
    echo SAMPLE = sample
    echo ANALYSIS_TYPE = [baseline,ipSig,datacardsPhiCP,jetFakes]
    echo APPLY_IPCUT1 = [0,1]
    echo APPLY_FF = [0,1]
    exit
fi

era=${1}
channel=${2}
sample=${3}
analysisType=${4}
ipcut=${5}
ff=${6}

folder=/afs/cern.ch/work/r/rasp/CPHiggs/Analysis
suffix=${era}_${channel}_${sample}_${analysisType}_${ipcut}_${ff}

cat > ${folder}/condor/job_${suffix}.submit <<EOF
universe = vanilla
executable = /afs/cern.ch/work/r/rasp/CPHiggs/Analysis/condor/job_${suffix}.sh
output = /afs/cern.ch/work/r/rasp/CPHiggs/Analysis/condor/job_${suffix}.out
error = /afs/cern.ch/work/r/rasp/CPHiggs/Analysis/condor/job_${suffix}.err
log = /afs/cern.ch/work/r/rasp/CPHiggs/Analysis/condor/job_${suffix}.log
notification = Never
initialdir = /afs/cern.ch/work/r/rasp/CMSSW_14_1_0_pre4/src/CPHiggs/Analysis
+MaxRuntime = 50000
+RequestRuntime = 50000
MY.WantOS = el9
queue
EOF

cat > ${folder}/condor/job_${suffix}.sh <<EOF1
#!/bin/tcsh
cd /afs/cern.ch/work/r/rasp/CMSSW_14_1_0_pre4/src
setenv SCRAM_ARCH el9_amd64_gcc10
export SCRAM_ARCH=el9_amd64_gcc10
cmsenv
cd CPHiggs/Analysis
echo $PWD
./scripts/RunSelection.py --era ${era} --channel ${channel} --sample ${sample} --analysisType ${analysisType} --applyFakeFactor ${ff} --applyIPSigLepCut ${ipcut}
EOF1

chmod u+x /afs/cern.ch/work/r/rasp/CPHiggs/Analysis/condor/job_${suffix}.sh
condor_submit /afs/cern.ch/work/r/rasp/CPHiggs/Analysis/condor/job_${suffix}.submit
