#!/bin/bash  

programDir=$1
zStart=$2
zStep=$3
numZsteps=$4
zEnd=`echo "scale=20; $zStart + ( $numZsteps - 1 ) * $zStep" | bc -l`

#pushd $workLoc/$currentDIR/Programs > /dev/null
pushd $programDir > /dev/null
echo $PWD

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
