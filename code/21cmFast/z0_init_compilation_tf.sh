#!/bin/bash

toDoList=$1
codeLoc=$2

paramSetNum=$(($TF_TASKID + 1))
currentLine=`cat $toDoList | head -$paramSetNum | tail -1`
currentLine=( $currentLine )
currentRandSeed=${currentLine[0]}
currentSigma8=${currentLine[1]}
currenthlittle=${currentLine[2]}
currentOmegam=${currentLine[3]}
currentOmegab=${currentLine[4]}
currentns=${currentLine[5]}
currentRmfp=${currentLine[6]}
currentTvir=${currentLine[7]}
currentSquiggly=${currentLine[8]}
currentlnAs=${currentLine[9]}

currentRandSeed="`printf '%d' $currentRandSeed`"
currentSigma8="`printf '%.5f' $currentSigma8`"
currenthlittle="`printf '%.5f' $currenthlittle`"
currentOmegam="`printf '%.5f' $currentOmegam`"
currentOmegab="`printf '%.5f' $currentOmegab`"
currentns="`printf '%.5f' $currentns`"
currentRmfp="`printf '%.2f' $currentRmfp`"
currentTvir="`printf '%.1f' $currentTvir`"
currentSquiggly="`printf '%.2f' $currentSquiggly`"
currentlnAs="`printf '%.5f' $currentlnAs`"

echo "Setting up (randSeed, sigma8, h, Omm, Omb, ns, Rmfp, Tvir, Squiggly, lnAs) = ($currentRandSeed, $currentSigma8, $currenthlittle, $currentOmegam, $currentOmegab, $currentns, $currentRmfp, $currentTvir, $currentSquiggly, $currentlnAs)..."
currentDIR="RandSeed_$currentRandSeed""_Sigma8_$currentSigma8""_h_$currenthlittle""_Omm_$currentOmegam""_Omb_$currentOmegab""_ns_$currentns""_Rmfp_$currentRmfp""_Tvir_$currentTvir""_Squiggly_$currentSquiggly""_lnAs_$currentlnAs"
if [ ! -d "$currentDIR" ]; then 
    mkdir $currentDIR
fi
pushd $currentDIR > /dev/null

### Now copy over 21cmFAST program files
cp -r $codeLoc/Programs .
cp -r $codeLoc/Cosmo_c_files .
cp -r $codeLoc/Parameter_files .
cp -r $codeLoc/External_tables .
mkdir Boxes
mkdir Log_files
mkdir Output_files

pushd Parameter_files > /dev/null

cat ANAL_PARAMS.H| sed s/"ION_Tvir_MIN (double)".*/"ION_Tvir_MIN (double) $currentTvir"/ > asdf.dat
cat asdf.dat | sed s/"HII_EFF_FACTOR (float)".*/"HII_EFF_FACTOR (float) $currentSquiggly"/ > asdf2.dat
cat asdf2.dat | sed s/"R_BUBBLE_MAX (float)".*/"R_BUBBLE_MAX (float) $currentRmfp"/ > asdf3.dat
rm -f asdf.dat asdf2.dat
rm -f ANAL_PARAMS.H
mv asdf3.dat ANAL_PARAMS.H
rm -f asdf3.dat

cat COSMOLOGY.H| sed s/"SIGMA8".*/"SIGMA8 ($currentSigma8)"/ > asdf.dat
cat asdf.dat | sed s/"hlittle (".*/"hlittle ($currenthlittle)"/ > asdf2.dat
cat asdf2.dat | sed s/"OMm  (".*/"OMm ($currentOmegam)"/ > asdf3.dat
cat asdf3.dat | sed s/"OMb  (float)".*/"OMb  (float) ($currentOmegab)"/ > asdf4.dat
cat asdf4.dat | sed s/"POWER_INDEX".*/"POWER_INDEX ($currentns)"/ > asdf5.dat
rm -f asdf.dat asdf2.dat asdf3.dat asdf4.dat
rm -f COSMOLOGY.H
mv asdf5.dat COSMOLOGY.H
rm -f asdf5.dat

cat INIT_PARAMS.H| sed s/"RANDOM_SEED (long) (".*/"RANDOM_SEED (long) ($currentRandSeed)"/ > asdf.dat
rm -f INIT_PARAMS.H
mv asdf.dat INIT_PARAMS.H
rm -f asdf.dat

pushd ../Programs > /dev/null

make clean > compilation.log 2>&1
make init > compilation.log 2>&1

popd > /dev/null
popd > /dev/null
popd > /dev/null

exit
