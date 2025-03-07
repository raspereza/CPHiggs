#!/bin/bash
# $1 - postfix

#for chan in mm ee
#do
#    rm selection/${chan}_Run3_2022All.root
#    rm selection/${chan}_Run3_2023All.root
#    hadd selection/${chan}_Run3_2022All.root selection/${chan}_Run3_2022.root selection/${chan}_Run3_2022EE.root
#    hadd selection/${chan}_Run3_2023All.root selection/${chan}_Run3_2023.root selection/${chan}_Run3_2023BPix.root    
#done

for chan in mt
do
    rm selection/${chan}_Run3_2022All_${1}.root
    rm selection/${chan}_Run3_2023All_${1}.root
    hadd selection/${chan}_Run3_2022All_${1}.root selection/${chan}_Run3_2022_${1}.root selection/${chan}_Run3_2022EE_${1}.root
    hadd selection/${chan}_Run3_2023All_${1}.root selection/${chan}_Run3_2023_${1}.root selection/${chan}_Run3_2023BPix_${1}.root    
done

#for chan in mt
#do
#    rm selection/${chan}_Run3_2022All.root
#    rm selection/${chan}_Run3_2023All.root
#    hadd selection/${chan}_Run3_2022All.root selection/${chan}_Run3_2022.root selection/${chan}_Run3_2022EE.root
#    hadd selection/${chan}_Run3_2023All.root selection/${chan}_Run3_2023.root selection/${chan}_Run3_2023BPix.root    
#done

