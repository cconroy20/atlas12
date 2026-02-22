#!/bin/bash

# calling sequence: 
#   atlas12.sh start end element abundance_change output_dir input_dir [element abundance_change]

# if element=99 -> vary alpha elements (O,Ne,NaMg,Si,S,Ca,Ti)
# if elem=99 is the second set of ab changes, then no Na
# if element=98 -> vary a subset of alpha elements (O,Ne,S)
# if element=1  -> vary Z.  (Units are log wrt Solar; e.g., "-1.0" generates a 1/10 Solar model)
# if element=2  -> vary He. (Units are linear)

#define the mixing length
mlt=2.03

#----------------------------------------------------------#

# Asplund et al. (2009) abundances of first 99 elements
# in ATLAS format (12.04)
#ee=(0.92068 0.07837 -10.99 -10.66 -9.34 -3.61 -4.21 -3.35 -7.48 -4.11 -5.80 -4.44 -5.59 -4.53 -6.63 -4.92 -6.54 -5.64 -7.01 -5.70 -8.89 -7.09 -8.11 -6.40 -6.61 -4.54 -7.05 -5.82 -7.85 -7.48 -9.00 -8.39 -9.74 -8.70 -9.50 -8.79 -9.52 -9.17 -9.83 -9.46 -10.58 -10.16 -20.00 -10.29 -11.13 -10.47 -11.10 -10.33 -11.24 -10.00 -11.03 -9.86 -10.49 -9.80 -10.96 -9.86 -10.94 -10.46 -11.32 -10.62 -20.00 -11.08 -11.52 -10.97 -11.74 -10.94 -11.56 -11.12 -11.94 -11.20 -11.94 -11.19 -12.16 -11.19 -11.78 -10.64 -10.66 -10.42 -11.12 -10.87 -11.14 -10.29 -11.39 -20.00 -20.00 -20.00 -20.00 -20.00 -20.00 -12.02 -20.00 -12.58 -20.00 -20.00 -20.00 -20.00 -20.00 -20.00 -20.00)

#GS98 abundances (from Dotter)
#changed Li to the A09 value (Dotter provided the meteoritic value)
ee=(0.92042 0.07834  -10.99 -10.62  -9.25  -3.52  -4.12  -3.21  -7.56  -3.96  -5.72  -4.46  -5.55  -4.48  -6.48  -4.84  -6.76  -5.64  -6.91  -5.69  -8.94  -7.10  -8.02  -6.35  -6.51  -4.54  -7.13  -5.79  -7.75  -7.37  -8.91  -8.41  -9.67  -8.63  -9.41  -8.73  -9.63  -9.12  -9.81  -9.43 -10.64 -10.07 -32.04 -10.21 -10.94 -10.34 -10.80 -10.28 -11.22  -9.90 -11.01  -9.80 -10.53  -9.87 -10.91  -9.82 -10.82 -10.41 -11.24 -10.55 -32.04 -11.06 -11.49 -10.95 -11.69 -10.87 -11.53 -11.07 -11.89 -11.08 -11.91 -11.29 -12.17 -11.35 -11.76 -10.65 -10.67 -10.35 -11.19 -10.91 -11.21  -9.98 -11.33 -32.04 -32.04 -32.04 -32.04 -32.04 -32.04 -11.95 -32.04 -12.54 -20.00 -20.00 -20.00 -20.00 -20.00 -20.00 -20.00)

# Change CNO to Grevesse & Sauval (1994)
#ee[6-1]=$(echo "scale=2;  0.12 + ${ee[6-1]}"  | bc)
#ee[7-1]=$(echo "scale=2;  $4 + ${ee[7-1]}"  | bc)
#ee[8-1]=$(echo "scale=2;  $4 + ${ee[8-1]}"  | bc)

# Caffau et al. (2011) scale (for Li,C,N,O,P,S,K,Fe,Eu,Hf,Os,Th)
#ee=(0.92068 0.07837 -11.01 -10.66 -9.34 -3.54 -4.18 -3.28 -7.48 -4.11 -5.80 -4.44 -5.59 -4.53 -6.58 -4.88 -6.54 -5.64 -6.93 -5.70 -8.89 -7.09 -8.11 -6.40 -6.61 -4.52 -7.05 -5.82 -7.85 -7.48 -9.00 -8.39 -9.74 -8.70 -9.50 -8.79 -9.52 -9.17 -9.83 -9.46 -10.58 -10.16 -20.00 -10.29 -11.13 -10.47 -11.10 -10.33 -11.24 -10.00 -11.03 -9.86 -10.49 -9.80 -10.96 -9.86 -10.94 -10.46 -11.32 -10.62 -20.00 -11.08 -11.52 -10.97 -11.74 -10.94 -11.56 -11.12 -11.94 -11.20 -11.94 -11.17 -12.16 -11.19 -11.78 -10.68 -10.66 -10.42 -11.12 -10.87 -11.14 -10.29 -11.39 -20.00 -20.00 -20.00 -20.00 -20.00 -20.00 -11.96 -20.00 -12.58 -20.00 -20.00 -20.00 -20.00 -20.00 -20.00 -20.00)

declare -i i1 i2 ie ii molflag dmflag

#----------------------Set Defaults------------------------#

# default metallicity abundance scale in units of Solar
zabnd=1.00000
# default Helium abundance in units of number fraction
heabnd=${ee[1]}
# default He+H abundance
tabnd=$(echo "scale=5; ${ee[0]} + ${ee[1]}" | bc)
# default H abundance
habnd=${ee[0]}
# set the input directory
ind=$6

# modify He abundance by hand
#heabnd=$(echo "scale=5; 0.0 * ${ee[1]}" | bc)
#habnd=$(echo "scale=5; $tabnd - $heabnd" | bc)

#-------------------Setup filenames and directories-----------------#

# output directory
outdir=$KURUCZ_HOME/grids/$5/

# create an output file head, stripping off the possible pre-directories
head=$5
head=${head#*/}
head=${head#*/}
head=${head#*/}
head=${head#*/}
head=${head#*/}

head2=$5
head2=${head2///}
head2=${head2///}
head2=${head2///}

# create the directory if it does not already exist
if [ ! -d $outdir ]; then
    echo "creating new directory:" $outdir
    mkdir -p $outdir
    mkdir $outdir/atm $outdir/molnden $outdir/outfiles $outdir/spec $outdir/taunu
    mkdir $outdir/flux $outdir/fail $outdir/crash $outdir/sed
fi

# input atmospheres
indir=$KURUCZ_HOME/grids/$ind/atm/

# get a list of the input atmospheres
arr=(`ls $indir`)
len=${#arr[*]}

# make sure that the input doesn't contain a plus sign
test=$4
if [ ${test:0:1} = "+" ]; then
    echo "syntax error: remove the +, exiting."
    exit
fi

# which model to analyze?
i1=$1-1
i2=$2-1
ie=1


#--------------Make Changes to the Abundance List-------------------#

# alpha abundance change (incl Na)
if [ $3 -eq 99 ]; then

    ee[8-1]=$(echo "scale=2;  $4 + ${ee[8-1]}"  | bc)
    ee[10-1]=$(echo "scale=2; $4 + ${ee[10-1]}" | bc)
    ee[11-1]=$(echo "scale=2; $4 + ${ee[11-1]}" | bc)
    ee[12-1]=$(echo "scale=2; $4 + ${ee[12-1]}" | bc)
    ee[14-1]=$(echo "scale=2; $4 + ${ee[14-1]}" | bc)
    ee[16-1]=$(echo "scale=2; $4 + ${ee[16-1]}" | bc)
    ee[20-1]=$(echo "scale=2; $4 + ${ee[20-1]}" | bc)
    ee[22-1]=$(echo "scale=2; $4 + ${ee[22-1]}" | bc)
    echo "[a/Fe]=$4"

elif [ $3 -eq 98 ]; then

    ee[8-1]=$(echo "scale=2;  $4 + ${ee[8-1]}"  | bc)
    ee[10-1]=$(echo "scale=2; $4 + ${ee[10-1]}" | bc)
    ee[16-1]=$(echo "scale=2; $4 + ${ee[16-1]}" | bc)
    echo "[O,Ne,S/Fe]=$4"

# change overall abundance scale
elif [ $3 -eq 1 ]; then

    #zabnd=$(echo "scale=5; e($4*l(10))" | bc -l)
    while [ $ie -le 98 ]; do
         ee[$ie]=$(echo "scale=2; $4 + ${ee[$ie]}" | bc)
	 let ie++
    done    

    echo "Overall abundance scaled by a factor of "$4 "with respect to Solar"

# change He abundance (and H)
elif [ $3 -eq 2 ]; then

    heabnd=$(echo "scale=5; $4 * ${ee[1]}" | bc)
    habnd=$(echo "scale=5; $tabnd - $heabnd" | bc)
    echo "He abundance changed to " $heabnd 
    echo "H  abundance changed to " $habnd

# Change abundance of a single element
else
    if [ $3 -gt 0 ]; then
	ee[$3-1]=$(echo "scale=2; ${ee[$3-1]} + $4" | bc)
	echo "Element "$3" changed by a factor of "$4
    fi
fi

# Option to include additional abundance variations
if [ $7 ]; then

    if [ $7 -eq 99 ]; then
	ee[8-1]=$(echo "scale=2;  $8 + ${ee[8-1]}"  | bc)
	ee[10-1]=$(echo "scale=2; $8 + ${ee[10-1]}" | bc)
	ee[12-1]=$(echo "scale=2; $8 + ${ee[12-1]}" | bc)
	ee[14-1]=$(echo "scale=2; $8 + ${ee[14-1]}" | bc)
	ee[16-1]=$(echo "scale=2; $8 + ${ee[16-1]}" | bc)
	ee[20-1]=$(echo "scale=2; $8 + ${ee[20-1]}" | bc)
	ee[22-1]=$(echo "scale=2; $8 + ${ee[22-1]}" | bc)
	echo "[a/Fe]=$8"
    elif [ $7 -eq 98 ]; then
	ee[8-1]=$(echo "scale=2;  $8 + ${ee[8-1]}"  | bc)
	ee[10-1]=$(echo "scale=2; $8 + ${ee[10-1]}" | bc)
	ee[16-1]=$(echo "scale=2; $8 + ${ee[16-1]}" | bc)
	echo "[O,Ne,S/Fe]=$8"
    elif [ $7 -eq 1 ]; then
	echo "no change to abundances"
    else
	ee[$7-1]=$(echo "scale=2; ${ee[$7-1]} + $8" | bc)
	echo "Element "$7" changed by a factor of "$8
    fi

fi

#------------------Now change X so that X+Y+Z=1-------------------#

nz=0.00
ii=3
while [ $ii -le 99 ]; do
    tt=$(echo "scale=4; ${ee[$ii-1]}+0.0" | bc)
    nz=$(echo "scale=12; $nz + e(${tt}*l(10))" | bc -l)
    let ii++
done

tot=0.0
tot=$(echo "scale=5; $habnd + $heabnd" | bc)
tot=$(echo "scale=5; $tot + 0${nz:0:6}" | bc)
tot=$(echo "scale=5; 1 - $tot " | bc)
# change H abund so that nz+ny+nz=1.0
# note that the He abundnace stays constant
habnd=$(echo "scale=5; $habnd + $tot" | bc)

# get the model filename
model=${arr[$i1]}
# we want a filename with all the header info stripped off
tmodel=${model#ap00}
#tmodel=${tmodel#at12_feh+0.00_afe+0.0_}
tmpd=$6
tmpd=${tmpd#c3k_v1.2/}
tmpd=${tmpd#c3k_v1.3/}
tmpd=${tmpd#c3k_v2.1/}
tmpd=${tmpd#c3k_v2.3/}
tmodel=${tmodel#$tmpd}
tmodel=${tmodel#Z+0.0}
tmodel=${tmodel#_}


if [ "${tmodel:6:1}" = "g" ]; then
    teff=${tmodel:1:5}
    logg=${tmodel:7:4}
else
    teff=${tmodel:1:4}
    logg=${tmodel:6:4}
fi

if [ ${teff} -lt 10000 ]; then
    if [ ${teff:0:1} -eq 0 ]; then
	teff=${teff:1:4}
    fi
    outfile=${head}_t0${teff}g${logg}
else
    outfile=${head}_t${teff}g${logg}
fi

logg=${logg}0

echo " "
echo $i1":"
date "+%Y-%m-%d %H:%M:%S"
echo "Input directory: "${indir}
echo "Input model: "${model}
echo "Output directory:"${outdir}
echo "Output model: "${outfile}
echo "Teff="$teff", logg="$logg

# check that the input atm model exists
inm="${indir}/${model}"
if [ ! -f $inm ]; then
    echo "error: input file does not exist"
    exit
fi

# check to see if the spec file already exists
# if it does, then we don't need to run this model
specf="${outdir}/sed/${outfile}.sed"
fluxf="${outdir}/flux/${outfile}.flux"
if [ -f $fluxf ]; then
    echo "flux file already exists, moving on..."
    exit
fi

# Delete an old version of the tmp dir if one exists
/bin/rm -rf tmp_${head2}.${1}
mkdir tmp_${head2}.${1}
echo "Moving into temporary working directory...."
cd tmp_${head2}.${1}

ln -s ${indir}${model} input_model.dat

# Check to see if the atm file already exists
if [ ! -f ${outdir}/atm/${outfile}.atm ]; then

# Run ATLAS12

$ATLAS12/bin/atlas12c.exe ${outfile} <<EOF> ${outfile}.out
MOLECULES ON
READ MOLECULES
READ PUNCH
READ LINES
TITLE ATLAS12 l/H=$mlt
OPACITY ON LINES
OPACITY ON XLINES
CONVECTION OVER $mlt 0 30
ITERATIONS 30
PRINT 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3
PUNCH 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2
SCALE MODEL 80 -6.875 0.125 ${teff}. ${logg}
VTURB 2.00E+05
ABUNDANCE SCALE   $zabnd ABUNDANCE CHANGE 1 $habnd 2 $heabnd
 ABUNDANCE CHANGE  3 ${ee[3-1]}  4 ${ee[4-1]}  5 ${ee[5-1]}  6 ${ee[6-1]}  7 ${ee[7-1]}  8 ${ee[8-1]}
 ABUNDANCE CHANGE  9 ${ee[9-1]} 10 ${ee[10-1]} 11 ${ee[11-1]} 12 ${ee[12-1]} 13 ${ee[13-1]} 14 ${ee[14-1]}
 ABUNDANCE CHANGE  15 ${ee[15-1]} 16 ${ee[16-1]} 17 ${ee[17-1]} 18 ${ee[18-1]} 19 ${ee[19-1]} 20 ${ee[20-1]}
 ABUNDANCE CHANGE  21 ${ee[21-1]} 22 ${ee[22-1]} 23 ${ee[23-1]} 24 ${ee[24-1]} 25 ${ee[25-1]} 26 ${ee[26-1]}
 ABUNDANCE CHANGE  27 ${ee[27-1]} 28 ${ee[28-1]} 29 ${ee[29-1]} 30 ${ee[30-1]} 31 ${ee[31-1]} 32 ${ee[32-1]}
 ABUNDANCE CHANGE  33 ${ee[33-1]} 34 ${ee[34-1]} 35 ${ee[35-1]} 36 ${ee[36-1]} 37 ${ee[37-1]} 38 ${ee[38-1]}
 ABUNDANCE CHANGE  39 ${ee[39-1]} 40 ${ee[40-1]} 41 ${ee[41-1]} 42 ${ee[42-1]} 43 ${ee[43-1]} 44 ${ee[44-1]}
 ABUNDANCE CHANGE  45 ${ee[45-1]} 46 ${ee[46-1]} 47 ${ee[47-1]} 48 ${ee[48-1]} 49 ${ee[49-1]} 50 ${ee[50-1]}
 ABUNDANCE CHANGE  51 ${ee[51-1]} 52 ${ee[52-1]} 53 ${ee[53-1]} 54 ${ee[54-1]} 55 ${ee[55-1]} 56 ${ee[56-1]}
 ABUNDANCE CHANGE  57 ${ee[57-1]} 58 ${ee[58-1]} 59 ${ee[59-1]} 60 ${ee[60-1]} 61 ${ee[61-1]} 62 ${ee[62-1]}
 ABUNDANCE CHANGE  63 ${ee[63-1]} 64 ${ee[64-1]} 65 ${ee[65-1]} 66 ${ee[66-1]} 67 ${ee[67-1]} 68 ${ee[68-1]}
 ABUNDANCE CHANGE  69 ${ee[69-1]} 70 ${ee[70-1]} 71 ${ee[71-1]} 72 ${ee[72-1]} 73 ${ee[73-1]} 74 ${ee[74-1]}
 ABUNDANCE CHANGE  75 ${ee[75-1]} 76 ${ee[76-1]} 77 ${ee[77-1]} 78 ${ee[78-1]} 79 ${ee[79-1]} 80 ${ee[80-1]}
 ABUNDANCE CHANGE  81 ${ee[81-1]} 82 ${ee[82-1]} 83 ${ee[83-1]} 84 ${ee[84-1]} 85 ${ee[85-1]} 86 ${ee[86-1]}
 ABUNDANCE CHANGE  87 ${ee[87-1]} 88 ${ee[88-1]} 89 ${ee[89-1]} 90 ${ee[90-1]} 91 ${ee[91-1]} 92 ${ee[92-1]}
 ABUNDANCE CHANGE  93 ${ee[93-1]} 94 ${ee[94-1]} 95 ${ee[95-1]} 96 ${ee[96-1]} 97 ${ee[85-1]} 98 ${ee[98-1]}
 ABUNDANCE CHANGE  99 ${ee[99-1]}
BEGIN 
END
EOF

#-----------------------------------------------------------------------#

date "+%Y-%m-%d %H:%M:%S"
echo "atlas is finished"

exit

# Check if the atm file was sucessfully created
if [ -f fort.7 ]; then
    #test convergence
    echo "testing output for convergence..."
    /n/home12/cconroy/kurucz/bin/checkconv.exe
    /bin/mv fort.50 $outdir/taunu/$taunu
    /bin/mv fort.66 $outdir/outfiles/$iter
    /bin/mv fort.67 $outdir/outfiles/$tcorr
    if [ -f model.failed ]; then
	echo "atlas12 model did not converge, exiting..."
        #save the atm file
        /bin/mv fort.7 $outdir/fail/$outfile
    else
        #save the atm file
        /bin/mv fort.7 $outdir/atm/$outfile
        #save the flux output
	if [ ! -d $outdir/flux ]; then
	    mkdir $outdir/flux
	fi
	flux="${outfile/atm/flux}"
	/bin/mv fort.8 $outdir/flux/$flux
        #clean up
	rm -f fort.* tmp.bin
    fi 
else
    echo "atlas12 apparently crashed, exiting..."
    /bin/cp fort.3 $outdir/crash/$model
    if [ -f fort.66 ]; then
	/bin/mv fort.66 $outdir/outfiles/$iter
    fi
fi

fi


# Clean up
echo "Removing working directory"
cd ../
date "+%Y-%m-%d %H:%M:%S"
