#!/bin/bash

#   parameterSweeps.sh
#   Author:   Gregory Kimmel
#   Date:     12/12/19
#
# This bash script calls go_or_grow to do parameter sweeps to calculate the
# degree of infiltration. We give two variable names, max, min and range to
# sweep over.
#
# NOTE: go_or_grow is called with usage:
#   ./go_or_grow <tFinal> <R> <R0> <dt> <dr> <u0> <v0> <d(eta)> <a> <xi_u>
#        <ku> <xi_v> <kv> <tfile> <rfile> <Efile> <ufile> <vfile> <saveFiles>
#        <verbosity>
# 
# USAGES (two possible):
#         ./parameterSweeps.sh p1name p1min p1steps p1max p2name
#             p2min p2steps p2max filename
#
#         ./parameterSweeps.sh p1name p1min p1steps p1max filename
#
#
# OUTPUTS
#   To terminal:
#       The output of Gillespie (e.g. mean first passage time)
#   To output files:
#       Saves the distribution of times to the two output files
#


# Check to make sure correct number of arguments are given
if [ "$#" -ne 9 ] && [ "$#" -ne 5 ]; then
  echo "Usage:" 
  echo "./exe p1name p1min p1steps p1max p2name p2min p2steps p2max filename"
  echo "./exe p1name p1min p1steps p1max filename"
  exit 1
fi

## Checking to make sure input arguments are correct ##
# This checks that values are numeric
if [[ -n ${2//[0-9,.]/} || ${3//[0-9,.]/} || ${4//[0-9,.]/} || ${6//[0-9,.]/} \
|| ${7//[0-9,.]/} || ${8//[0-9,.]/} ]]; then
    echo "Starting value should be a number!"
    exit 2
fi

if [ "$#" = 9 ]; then
  outfile=$9
fi
if [ "$#" = 5 ]; then
  outfile=$5
fi

# Default parameter values
tfinal=1500.0
R=10.0
R0=3.0
dr=0.1
dt=0.01
eta=10000.0
a=3.5
xi_u=1000.0
xi_v=1.2
ku=0.05
kv=0.2
u0=1.0
v0=0.0

# Place all variables into an array
allVars=(tfinal R R0 dr dt eta a xi_u xi_v ku kv u0 v0)

# Set empty array to be filled with the variables that are changed
varsVaried=()

# Loop through all variable names to find which ones are selected
for i in "${!allVars[@]}"; do
    if [[ "${allVars[$i]}" = "$1" ]]; then
        # echo "${allVars[$i]}";
        varsVaried+=("${allVars[$i]}")
    fi
done

# Loop through all variable names to find which ones are selected
for i in "${!allVars[@]}"; do
    if [[ "${allVars[$i]}" = "$5" ]]; then
        # echo "${allVars[$i]}";
        varsVaried+=("${allVars[$i]}")
    fi
done

echo "Iterating over ${#varsVaried[@]} variable(s): ${varsVaried[@]}"

# Loop over all possible values
if [ "${#varsVaried[@]}" -eq 2 ]; then
  for val1 in `seq $2 $3 $4`
  do
      for val2 in `seq $6 $7 $8`
      do
          # Evaluate the requested variables by the values
          eval ${varsVaried[0]}=${val1}
          eval ${varsVaried[1]}=${val2}

          # Check to make sure the correct variables are updated
          # echo ./go_or_grow ${tfinal} ${R} ${R0} ${dt} ${dr} ${u0} ${v0} ${eta} ${a} ${xi_u} ${ku} ${xi_v} ${kv} t.txt r.txt E.txt u.txt v.txt ${9} 0 1

          ./go_or_grow_nocutoff ${tfinal} ${R} ${R0} ${dt} ${dr} ${u0} ${v0} ${eta} ${a} ${xi_u} ${ku} ${xi_v} ${kv} td.txt rd.txt Ed.txt ud.txt v.txt ${outfile} 1 1 0.1

      done
  done
elif [ "${#varsVaried[@]}" -eq 1 ]; then
  for val1 in `seq $2 $3 $4`
  do
        # Evaluate the requested variables by the values
        eval ${varsVaried[0]}=${val1}

        # Check to make sure the correct variables are updated
        # echo ./go_or_grow_nocutoff ${tfinal} ${R} ${R0} ${dt} ${dr} ${u0} ${v0} ${eta} ${a} ${xi_u} ${ku} ${xi_v} ${kv} td.txt rd.txt Ed.txt ud.txt v.txt $5 1 1 0.1

        ./go_or_grow_nocutoff ${tfinal} ${R} ${R0} ${dt} ${dr} ${u0} ${v0} ${eta} ${a} ${xi_u} ${ku} ${xi_v} ${kv} td.txt rd.txt Ed.txt ud.txt v.txt ${outfile} 0 1 0.1
  done
fi