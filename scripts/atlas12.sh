#!/bin/bash
#set -euo pipefail

# Calling sequence:
#   atlas12.sh model_number element abundance_change output_dir input_dir
#
# If element=99 -> vary alpha elements (O, Ne, Na, Mg, Si, S, Ca, Ti)
# If element=98 -> vary a subset of alpha elements (O, Ne, S)
# If element=1  -> vary Z  (units are log wrt Solar; e.g., "-1.0" generates a 1/10 Solar model)
# If element=2  -> vary He (units are linear)

if [ $# -lt 5 ]; then
    echo "Usage: $(basename "$0") model_number element abundance_change output_dir input_dir"
    exit 1
fi

# Assign named variables from positional parameters
model_number="$1"
element="$2"
abund_change="$3"
output_dir="$4"
input_dir="$5"

# Define the mixing length
mlt=2.03

#----------------------------------------------------------#

# Select abundance scale: GS98, A09, or C11
abund_scale="GS98"

case "$abund_scale" in

    # Grevesse & Sauval (1998), from Dotter
    # Li changed to the A09 value (Dotter provided the meteoritic value)
    GS98)
    ee=(0.92042 0.07834 \
        -10.99 -10.62  -9.25  -3.52  -4.12  -3.21  -7.56  -3.96 \
         -5.72  -4.46  -5.55  -4.48  -6.48  -4.84  -6.76  -5.64 \
         -6.91  -5.69  -8.94  -7.10  -8.02  -6.35  -6.51  -4.54 \
         -7.13  -5.79  -7.75  -7.37  -8.91  -8.41  -9.67  -8.63 \
         -9.41  -8.73  -9.63  -9.12  -9.81  -9.43 -10.64 -10.07 \
        -32.04 -10.21 -10.94 -10.34 -10.80 -10.28 -11.22  -9.90 \
        -11.01  -9.80 -10.53  -9.87 -10.91  -9.82 -10.82 -10.41 \
        -11.24 -10.55 -32.04 -11.06 -11.49 -10.95 -11.69 -10.87 \
        -11.53 -11.07 -11.89 -11.08 -11.91 -11.29 -12.17 -11.35 \
        -11.76 -10.65 -10.67 -10.35 -11.19 -10.91 -11.21  -9.98 \
        -11.33 -32.04 -32.04 -32.04 -32.04 -32.04 -32.04 -11.95 \
        -32.04 -12.54 -20.00 -20.00 -20.00 -20.00 -20.00 -20.00 \
        -20.00)
    ;;

    # Asplund et al. (2009), ATLAS format (12.04)
    A09)
    ee=(0.92068 0.07837 \
        -10.99 -10.66  -9.34  -3.61  -4.21  -3.35  -7.48  -4.11 \
         -5.80  -4.44  -5.59  -4.53  -6.63  -4.92  -6.54  -5.64 \
         -7.01  -5.70  -8.89  -7.09  -8.11  -6.40  -6.61  -4.54 \
         -7.05  -5.82  -7.85  -7.48  -9.00  -8.39  -9.74  -8.70 \
         -9.50  -8.79  -9.52  -9.17  -9.83  -9.46 -10.58 -10.16 \
        -20.00 -10.29 -11.13 -10.47 -11.10 -10.33 -11.24 -10.00 \
        -11.03  -9.86 -10.49  -9.80 -10.96  -9.86 -10.94 -10.46 \
        -11.32 -10.62 -20.00 -11.08 -11.52 -10.97 -11.74 -10.94 \
        -11.56 -11.12 -11.94 -11.20 -11.94 -11.19 -12.16 -11.19 \
        -11.78 -10.64 -10.66 -10.42 -11.12 -10.87 -11.14 -10.29 \
        -11.39 -20.00 -20.00 -20.00 -20.00 -20.00 -20.00 -12.02 \
        -20.00 -12.58 -20.00 -20.00 -20.00 -20.00 -20.00 -20.00 \
        -20.00)
    ;;

    # Caffau et al. (2011) (for Li,C,N,O,P,S,K,Fe,Eu,Hf,Os,Th)
    C11)
    ee=(0.92068 0.07837 \
        -11.01 -10.66  -9.34  -3.54  -4.18  -3.28  -7.48  -4.11 \
         -5.80  -4.44  -5.59  -4.53  -6.58  -4.88  -6.54  -5.64 \
         -6.93  -5.70  -8.89  -7.09  -8.11  -6.40  -6.61  -4.52 \
         -7.05  -5.82  -7.85  -7.48  -9.00  -8.39  -9.74  -8.70 \
         -9.50  -8.79  -9.52  -9.17  -9.83  -9.46 -10.58 -10.16 \
        -20.00 -10.29 -11.13 -10.47 -11.10 -10.33 -11.24 -10.00 \
        -11.03  -9.86 -10.49  -9.80 -10.96  -9.86 -10.94 -10.46 \
        -11.32 -10.62 -20.00 -11.08 -11.52 -10.97 -11.74 -10.94 \
        -11.56 -11.12 -11.94 -11.20 -11.94 -11.17 -12.16 -11.19 \
        -11.78 -10.68 -10.66 -10.42 -11.12 -10.87 -11.14 -10.29 \
        -11.39 -20.00 -20.00 -20.00 -20.00 -20.00 -20.00 -11.96 \
        -20.00 -12.58 -20.00 -20.00 -20.00 -20.00 -20.00 -20.00 \
        -20.00)
    ;;

    *)
    echo "error: unknown abundance scale '$abund_scale' (use GS98, A09, or C11)"
    exit 1
    ;;
esac

echo ""
date "+%Y-%m-%d %H:%M:%S"
echo "Abundance scale: $abund_scale"

#----------------------Set Defaults------------------------#

# Lightweight floating-point helper: one awk process per call.
# Usage: _calc "expression" [decimal_places]
_calc() {
    awk -v PREC="${2:-5}" "BEGIN { printf \"%.*f\n\", PREC, $1 }"
}

# Default metallicity abundance scale in units of Solar
zabnd=1.00000
# Default Helium abundance in units of number fraction
heabnd=${ee[1]}
# Default He+H abundance
tabnd=$(_calc "${ee[0]} + ${ee[1]}" 5)
# Default H abundance
habnd=${ee[0]}

#-------------------Setup filenames and directories-----------------#

outdir="${KURUCZ_HOME}/grids/${output_dir}/"
indir="${KURUCZ_HOME}/grids/${input_dir}/atm/"

# Output filename prefix: basename of output_dir
head="${output_dir##*/}"
# Temp directory name: output_dir with slashes removed, plus model number
tmpdir="tmp_${output_dir//\//}.${model_number}"

# Create the output directory tree
mkdir -p "${outdir}"/{atm,molnden,outfiles,spec,taunu,flux,fail,crash,sed}

# Get a list of the input atmospheres
arr=( $(ls "$indir") )

# Make sure that the abundance change doesn't start with a plus sign
if [ "${abund_change:0:1}" = "+" ]; then
    echo "syntax error: remove the + from abundance_change, exiting."
    exit 1
fi

# Convert 1-based model number to 0-based index
model_idx=$(( model_number - 1 ))

#--------------Make Changes to the Abundance List-------------------#

if [ "$element" -eq 99 ]; then

    # Alpha abundance change (O, Ne, Na, Mg, Si, S, Ca, Ti)
    for idx in 7 9 10 11 13 15 19 21; do
        ee[$idx]=$(_calc "$abund_change + ${ee[$idx]}" 2)
    done
    echo "[a/Fe]=$abund_change"

elif [ "$element" -eq 98 ]; then

    # Subset of alpha elements (O, Ne, S)
    for idx in 7 9 15; do
        ee[$idx]=$(_calc "$abund_change + ${ee[$idx]}" 2)
    done
    echo "[O,Ne,S/Fe]=$abund_change"

elif [ "$element" -eq 1 ]; then

    # Change overall metallicity: add delta to elements 1..98 (0-indexed)
    read -ra ee <<< "$(printf '%s\n' "${ee[@]}" | awk -v delta="$abund_change" '
        NR == 1 { print; next }
        NR >= 2 && NR <= 99 { printf "%.2f\n", $1 + delta; next }
        { print }
    ' | tr '\n' ' ')"
    echo "Overall abundance scaled by a factor of $abund_change with respect to Solar"

elif [ "$element" -eq 2 ]; then

    # Change He abundance (and adjust H)
    heabnd=$(_calc "$abund_change * ${ee[1]}" 5)
    habnd=$(_calc "$tabnd - $heabnd" 5)
    echo "He abundance changed to $heabnd"
    echo "H  abundance changed to $habnd"

elif [ "$element" -gt 0 ]; then

    # Change abundance of a single element
    ee[$element-1]=$(_calc "${ee[$element-1]} + $abund_change" 2)
    echo "Element $element changed by a factor of $abund_change"

fi

#-------------------Extract model parameters-----------------------#

model="${arr[$model_idx]}"

# Extract Teff and logg from the filename (e.g. "...t05500g3.50.atm")
if [[ "$model" =~ t([0-9]+)g([0-9]+\.[0-9]+) ]]; then
    teff="${BASH_REMATCH[1]#0}"   # strip leading zero
    logg="${BASH_REMATCH[2]}0"    # ATLAS12 expects trailing zero (e.g. 3.500)
else
    echo "error: cannot parse Teff/logg from filename: $model"
    exit 1
fi
outfile="${head}_t$(printf '%05d' "$teff")g${logg%0}"

echo "Input:  ${indir}${model}"
echo "Output: ${outdir}atm/${outfile}.atm"
echo "Teff=${teff}, logg=${logg}"

# Check that the input atm model exists
if [ ! -f "${indir}/${model}" ]; then
    echo "error: input file does not exist: ${indir}/${model}"
    exit 1
fi

#-------------------Set up working directory-----------------------#

/bin/rm -rf "$tmpdir"
mkdir "$tmpdir"
cd "$tmpdir"
ln -s "${indir}${model}" input_model.dat

#-------------------Generate input abundance file------------------#

# Write abundance file for elements 3-99 from the ee array.
# Format: Z  log10(number_fraction), one element per line.
abund_file="${outfile}.abund"
{
  for ((i = 3; i <= 99; i++)); do
    printf "%2d  %s\n" "$i" "${ee[$((i - 1))]}"
  done
} > "$abund_file"

#-------------------Run ATLAS12------------------------------------#

if [ ! -f "${outdir}/atm/${outfile}.atm" ]; then

${ATLAS12}/bin/atlas12c.exe ${outfile} numit=30 vturb=2.0 mlt=2.03 \
    teff=${teff} logg=${logg} zscale=${zabnd} heabnd=${heabnd} \
    abund=${abund_file} > ${outfile}.out

date "+%Y-%m-%d %H:%M:%S"
echo "atlas is finished"

exit

# Check if the atm file was successfully created
if [ -f fort.7 ]; then
    echo "testing output for convergence..."
    "${ATLAS12}/bin/checkconv.exe"
    if [ -f model.failed ]; then
        echo "atlas12 model did not converge"
        /bin/mv fort.7 "${outdir}/fail/${outfile}.atm"
    else
        /bin/mv fort.7 "${outdir}/atm/${outfile}.atm"
        /bin/mv fort.8 "${outdir}/flux/${outfile}.flux"
        rm -f fort.* tmp.bin
    fi
else
    echo "atlas12 apparently crashed"
    /bin/cp fort.3 "${outdir}/crash/${model}"
fi

fi

date "+%Y-%m-%d %H:%M:%S"

echo "Removing working directory"
cd ../
/bin/rm -rf "${tmpdir}"
