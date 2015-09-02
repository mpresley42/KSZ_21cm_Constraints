#!/bin/bash

toDoList=$1
zStart=$2
zStep=$3
numZsteps=$4
zEnd=`echo "scale=20; $zStart + ( $numZsteps - 1 ) * $zStep" | bc -l`

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

echo "Evolving density field for (randSeed, sigma8, h, Omm, Omb, ns, Rmfp, Tvir, Squiggly, lnAs) = ($currentRandSeed, $currentSigma8, $currenthlittle, $currentOmegam, $currentOmegab, $currentns, $currentRmfp, $currentTvir, $currentSquiggly, $currentlnAs)..."
currentDIR="RandSeed_$currentRandSeed""_Sigma8_$currentSigma8""_h_$currenthlittle""_Omm_$currentOmegam""_Omb_$currentOmegab""_ns_$currentns""_Rmfp_$currentRmfp""_Tvir_$currentTvir""_Squiggly_$currentSquiggly""_lnAs_$currentlnAs"
pushd $currentDIR/Programs > /dev/null

cat drive_zscroll_noTs_evolveDensity.c | sed s/"ZSTART (".*/"ZSTART ($zStart)"/ > asdf.dat
cat asdf.dat | sed s/"ZEND (".*/"ZEND ($zEnd)"/ > asdf2.dat
cat asdf2.dat | sed s/"ZSTEP (".*/"ZSTEP ($zStep)"/ > asdf3.dat
rm -f asdf.dat asdf2.dat
rm -f drive_zscroll_noTs_evolveDensity.c
mv asdf3.dat drive_zscroll_noTs_evolveDensity.c
rm -f asdf3.dat

cat drive_zscroll_noTs_reion.c | sed s/"ZSTART (".*/"ZSTART ($zStart)"/ > asdf.dat
cat asdf.dat | sed s/"ZEND (".*/"ZEND ($zEnd)"/ > asdf2.dat
cat asdf2.dat | sed s/"ZSTEP (".*/"ZSTEP ($zStep)"/ > asdf3.dat
rm -f asdf.dat asdf2.dat
rm -f drive_zscroll_noTs_reion.c
mv asdf3.dat drive_zscroll_noTs_reion.c
rm -f asdf3.dat

cat drive_zscroll_noTs_noTb_reion.c | sed s/"ZSTART (".*/"ZSTART ($zStart)"/ > asdf.dat
cat asdf.dat | sed s/"ZEND (".*/"ZEND ($zEnd)"/ > asdf2.dat
cat asdf2.dat | sed s/"ZSTEP (".*/"ZSTEP ($zStep)"/ > asdf3.dat
rm -f asdf.dat asdf2.dat
rm -f drive_zscroll_noTs_noTb_reion.c
mv asdf3.dat drive_zscroll_noTs_noTb_reion.c
rm -f asdf3.dat

#make clean > compilation.log 2>&1
make drive_zscroll_noTs_reion >& compilation.log 2>&1
make drive_zscroll_noTs_noTb_reion >& compilation.log 2>&1
make drive_zscroll_noTs_evolveDensity >& compilation.log 2>&1
make perturb_field >& compilation.log 2>&1
make delta_T >& compilation.log 2>&1
make find_HII_bubbles >& compilation.log 2>&1
            
popd > /dev/null
