#!/bin/sh

dirList="PRODFSR PRODFSR_8TeV JHU JHU_8TeV"

#make directories if they don't exist

if [[ ! -d 7plus8TeV_FSR ]] 
    then
    mkdir 7plus8TeV_FSR
fi
if [[ ! -d JHUsignal ]]
    then
    mkdir JHUsignal
fi

#####################################

for dir in $dirList
  do

  if [[ -d $dir ]] 
      then

      cd $dir

      #check whether files in dir are 7TeV or 8TeV

      if [[ $dir == *8TeV* ]] 
	  then 
	  append=8TeV
	  else
	  append=7TeV
      fi

      ############################################

      #move files for channel sub-directory and change 4l to corresponding channel name

      for i in $(ls 2mu2e/HZZ4lTree*.root | awk -F _ '{print $2}' | awk -F . '{print $1}')
	do 
	mv 2mu2e/HZZ4lTree_${i}.root HZZ2e2muTree_${i}_${append}.root
      done
      
      for i in $(ls 4mu/HZZ4lTree*.root | awk -F _ '{print $2}' | awk -F . '{print $1}')
	do 
	mv 4mu/HZZ4lTree_${i}.root HZZ4muTree_${i}_${append}.root
      done
      
      for i in $(ls 4e/HZZ4lTree*.root | awk -F _ '{print $2}' | awk -F . '{print $1}')
	do 
	mv 4e/HZZ4lTree_${i}.root HZZ4eTree_${i}_${append}.root
      done
      
      rm -r 4mu
      rm -r 4e 
      rm -r 2mu2e

      cd -

      ###################################################################################

      #move files to combine 7 plus 8 TeV directory

      if [[ $dir = *JHU* ]]
	  then
	  mv ${dir}/* JHUsignal/.
	  rm -r $dir
	  else
	  mv ${dir}/* 7plus8TeV_FSR/.
	  rm -r $dir
      fi

      #############################################

      fi
done

exit
