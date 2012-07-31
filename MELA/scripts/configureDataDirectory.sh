#!/bin/sh

mkdir 7TeVplus8TeV_FSR

list="PRODFSR_8TeV PRODFSR"

for dir in $list
do 
for i in $(ls $dir/2mu2e/HZZ4lTree*.root | awk -F _ '{print $2}')
    do mv $dir/2mu2e/HZZ4lTree_$i 7TeVplus8TeV_FSR/HZZ2e2muTree_$i
    done
    
for i in $(ls $dir/4mu/HZZ4lTree*.root | awk -F _ '{print $2}')
    do mv $dir/4mu/HZZ4lTree_$i 7TeVplus8TeV_FSR/HZZ4muTree_$i
    done

for i in $(ls $dir/4e/HZZ4lTree*.root | awk -F _ '{print $2}')
    do mv $dir/4e/HZZ4lTree_$i 7TeVplus8TeV_FSR/HZZ4eTree_$i
    done
done

exit
