#!/bin/bash

# script to run at atmosphere through synthe

# calling sequence: synthe.sh directory input_atm

echo " "
date

#set input and output directories
kgrids=/Users/cconroy/kurucz/
indir="$kgrids/grids/$1/atm/"
outdir="$kgrids/grids/$1/spec/"
echo " Output Dir: ${outdir}"

#check if the directory exists
if [ ! -d $indir ]; then
    echo $indir
    echo "synthe error: input directory does not exist!"
    exit
fi
if [ ! -d $outdir ]; then
    mkdir $outdir
fi

rundir=${ATLAS12}/workdir/$1_$2
#delete an old version of the tmp dir if one exists
/bin/rm -rf $rundir

mkdir $rundir
cd $rundir

#run synthe
echo " Running synthe...."
$ATLAS12/bin/synthe.exe $indir$2 wlbeg=350 wlend=1500 resolu=300000

#save the molecular number density profiles
#moldir="$kgrids/grids/$1/molnden/"
#mol="${asc/spec/mol}"
#/bin/mv fort.35 ${moldir}/${mol}
#save the line formation depths
#/bin/mv *linform ${moldir}/

echo " "
echo " ...done!"
cd ../
#/bin/rm -rf $rundir

date
