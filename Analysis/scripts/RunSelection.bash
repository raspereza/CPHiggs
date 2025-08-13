#!/bin/bash

n=$#
if [[ $n -ne 9 ]]; then
    echo usage : RunSelection.bash [ERA] [CHANNEL] [SAMPLE] [ANALYSIS_TYPE] [APPLY_CROSS_TRIGGER] [APPLY_IPCUT1] [APPLY_IPCUT2] [APPLY_PROMPT_SF] [APPLY_TAU_SF]
    echo ERA = [Run3_2022preEE, Run3_2022postEE, Run3_2022preBPix, Run3_2022postBPix]
    echo CHANNEL = [mt,et,mm,ee]
    echo SAMPLE = sample
    echo ANALYSIS_TYPE = [baseline,ipSig,phiCP]
    echo APPLY_CROSS_TRIGGER = [0,1]
    echo APPLY_IPCUT1 = [0,1]
    echo APPLY_IPCUT2 = [0,1]
    echo APPLY_PROMPT_SF = [0,1]
    echo APPLY_TAU_SF = [0,1]
    exit
fi

era=$1
channel=$2
sample=$3
analysisType=$4
xtrig=$5
ipcut1=$6
ipcut2=$7
promptSF=$8
tauSF=$9

folder=/afs/cern.ch/work/r/rasp/CPHiggs/Analysis
suffix=${era}_${channel}_${sample}_${analysisType}_${xtrig}_${ipcut1}_${ipcut2}_${promptSF}_${tauSF}

cat > ${folder}/condor/job_${suffix}.submit <<EOF
universe = vanilla
executable = /afs/cern.ch/work/r/rasp/CPHiggs/Analysis/condor/job_${suffix}.sh
output = /afs/cern.ch/work/r/rasp/CPHiggs/Analysis/condor/job_${suffix}.out
error = /afs/cern.ch/work/r/rasp/CPHiggs/Analysis/condor/job_${suffix}.err
log = /afs/cern.ch/work/r/rasp/CPHiggs/Analysis/condor/job_${suffix}.log
notification = Never
initialdir = /afs/cern.ch/work/r/rasp/CMSSW_14_1_0_pre4/src/CPHiggs/Analysis
+MaxRuntime = 20000
+RequestRuntime = 20000
MY.WantOS = el9
queue
EOF

cat > /afs/cern.ch/work/r/rasp/CPHiggs/Analysis/condor/job_${suffix}.sh <<EOF1
#!/bin/tcsh
cd /afs/cern.ch/work/r/rasp/CMSSW_14_1_0_pre4/src
setenv SCRAM_ARCH el9_amd64_gcc10
export SCRAM_ARCH=el9_amd64_gcc10
cmsenv
cd CPHiggs/Analysis
echo $PWD
./scripts/RunSelection.py --era ${era} --channel ${channel} --sample ${sample} --analysisType ${analysisType} --useCrossTrigger ${xtrig} --applyIPSigLep1Cut ${ipcut1} --applyIPSigLep2Cut ${ipcut2} --applyIPSigPromptLepSF ${promptSF} --applyIPSigTauLepSF ${tauSF}
EOF1

chmod u+x /afs/cern.ch/work/r/rasp/CPHiggs/Analysis/condor/job_${suffix}.sh
condor_submit /afs/cern.ch/work/r/rasp/CPHiggs/Analysis/condor/job_${suffix}.submit
