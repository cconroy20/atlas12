#!/bin/bash

# script to run a directory full of atmospheres through synthe

# calling sequence: synthe.sh start end output_dir [input_atm]

date

tout=$3
tout="${tout///}"

bindir="$ATLAS12/bin/"

#set input and output directories
indir="${KURUCZ_HOME}/grids/$3/atm/"
outdir="${KURUCZ_HOME}/grids/$3/spec/"

#test if the directory exists
if [ ! -d $indir ]; then
    echo $indir
    echo "synthe error: input directory does not exist!"
    exit
fi
if [ ! -d $outdir ]; then
    mkdir $outdir
fi

moldir="${KURUCZ_HOME}/grids/$3/molnden/"
echo "Output Dir: ${outdir}"

arr=(`ls $indir`)
len=${#arr[*]}

# model(s) to analyze
declare -i i1 i2
if [ ! $4 ]; then
    i1=$1-1 
    i2=$2-1
else
    i1=$1
    i2=$2
fi

rundir=${ATLAS12}/workdir/${tout}_$1

#delete an old version of the tmp dir if one exists
/bin/rm -rf $rundir

while [ $i1 -le $i2 ]; do

echo "Moving into temporary working directory...."
mkdir $rundir
cd $rundir

# Run the specific model passed, if $4 exists
if [ $4 ]; then
    model=$4
else
    model=${arr[$i1]}
fi

# Set the directory containing the line information
linedir="Lines_RV31new"


asc="${model/atm/spec}"
asc="${asc/dat/spec}"

#only run the model if no spec file exists
if [ ! -f "${outdir}/${asc}" ] && [ ! -f "${outdir}/${asc}.gz" ]; then

    echo "$i1: $model -> $asc"
    echo $linedir

    #generate synthe-ready input file
    $bindir/at12tosyn.exe $indir$model $model

    #link the input files generated from synthe.setup
    ln -s ${ATLAS12}/${linedir}/tfort.12 fort.12
    ln -s ${ATLAS12}/${linedir}/tfort.14 fort.14
    ln -s ${ATLAS12}/${linedir}/tfort.19 fort.19
    ln -s ${ATLAS12}/${linedir}/tfort.20 fort.20
    ln -s ${ATLAS12}/${linedir}/tfort.93 fort.93

    ln -s $model fort.5

    #run synthe, the main program
    echo "running synthe...."
    $bindir/synthe_spectrv.exe 

    #save the molecular number density profiles
    #mol="${asc/spec/mol}"
    #/bin/mv fort.35 ${moldir}/${mol}

    #save the line formation depths (not currently)
    #lin="${asc/spec/linform}"
    #/bin/mv fort.33 ${moldir}/${lin}

    #/bin/mv fort.2 $outdir$asc

fi

echo "this model is done"

echo "removing working directory"
cd ../
#/bin/rm -rf $rundir

let i1++

done

date

