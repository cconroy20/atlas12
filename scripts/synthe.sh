#!/bin/bash

# script to run at atmosphere through synthe

# calling sequence: synthe.sh directory input_atm

echo " "
date

bindir="$ATLAS12/bin/"

kgrids=/Users/cconroy/kurucz/

#set input and output directories
indir="$kgrids/grids/$1/atm/"
outdir="$kgrids/grids/$1/spec/"
moldir="$kgrids/grids/$1/molnden/"

#check if the directory exists
if [ ! -d $indir ]; then
    echo $indir
    echo "synthe error: input directory does not exist!"
    exit
fi
if [ ! -d $outdir ]; then
    mkdir $outdir
fi

echo " Output Dir: ${outdir}"

arr=(`ls $indir`)
len=${#arr[*]}

rundir=${ATLAS12}/workdir/$1_$2

#delete an old version of the tmp dir if one exists
/bin/rm -rf $rundir

mkdir $rundir
cd $rundir

# Set the directory containing the line information
linedir="Lines_RV31new"
echo " Input line dir:" $linedir

#generate synthe-ready input file
/bin/cp $indir$2 $2

#link the input files generated from synthe.setup
ln -s ${ATLAS12}/${linedir}/tfort.12 fort.12
ln -s ${ATLAS12}/${linedir}/tfort.14 fort.14
ln -s ${ATLAS12}/${linedir}/tfort.19 fort.19
ln -s ${ATLAS12}/${linedir}/tfort.20 fort.20
ln -s ${ATLAS12}/${linedir}/tfort.93 fort.93

#run synthe, the main program
echo " Running synthe...."
$bindir/synthe.exe $2

#save the molecular number density profiles
#mol="${asc/spec/mol}"
#/bin/mv fort.35 ${moldir}/${mol}

#save the line formation depths (not currently)
#lin="${asc/spec/linform}"
#/bin/mv fort.33 ${moldir}/${lin}

echo " "
echo " ...done!"
cd ../
#/bin/rm -rf $rundir

date

