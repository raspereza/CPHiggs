# C++ Wrapper of fastMTT

## Content

The subpackage `CPHiggs/fastMTT` contains the following components:
* `fastmtt_cpp.cpp` - C++ wrapper of fastMTT algorithm
* `compile_fastmtt.bash` - compilation script (creates shared library)
* `check_fastmtt.py` - testing script
* `pybind11` - folder of pybind11 package (C++ binding to python)
* `PlotFastMTT.py` - plotting macro
* `ntuple.h` - content of the IC tuple (produced with the HiggsDNA package)

## Getting code from git

```
cd $CMSSW_BASE/src
git clone https://github.com/raspereza/CPHiggs.git
cd $CMSSW_BASE/src/CPHiggs
scramv1 b -j 4
```

## Compiling shared library

```
cd $CMSSW_BASE/src/CPHiggs/fastMTT
./compile_fastmtt.bash
```

## Running test

Example:
```
./check_fastmtt.py --era Run3_2022 --channel mt --sample higgs --nevts 1000000 
```
The script will create RooT file named `higgs_Run3_2022_mt.root` in the directory `$CMSSW_BASE/src/CPHiggs/fastMTT`.

## Plotting output of testing script

Example:
```
./PlotFastMTT.py --era Run3_2022 --channel mt --sample higgs
```
The script will create the following files
* `higgs_Run3_2022_mt_mass.png` - distributions of `m_vis` and `m_fastmtt`;
* `higgs_Run3_2022_mt_dpt1.png` - distribution of pT(1,reco)/pT(1,gen) : ratio of the first tau pT reconstructed using fastMTT (tau->mu in this particular example) and generated pT;
* `higgs_Run3_2022_nt_dpt2.png` - distribution of pT(2,reco)/pT(2,gen) : ratio of the second tau pT reconstructed using fastMTT (tauh in this particular example) and generated pT;

