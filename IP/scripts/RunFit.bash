#!/bin/bash
# $1 - era
# $2 - chan
era=$1
chan=$2
cd datacards
for pt in 1 2 3 4
do
    for eta in 1 2
    do
	name=${chan}_${era}_binPt${pt}_binEta${eta}
	combineCards.py ${name}_pass.txt ${name}_fail.txt > ${name}.txt
	combineTool.py -M T2W -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO '"map=^.*/*_pass:r_pass[1,0,2]"' --PO '"map=^.*/*_fail:r_fail[1,0,2]"' -o ${name}.root -i ${name}.txt 
	combineTool.py -M FitDiagnostics --saveNormalizations --saveShapes --saveWithUncertainties --saveNLL --redefineSignalPOIs r_pass,r_fail --robustFit 1 -m 91 -d ${name}.root --cminDefaultMinimizerTolerance 0.1 --cminDefaultMinimizerStrategy 1 -v 2
	mv fitDiagnostics.Test.root ${name}_fit.root
	rm higgsCombine*

    done
done
cd -
