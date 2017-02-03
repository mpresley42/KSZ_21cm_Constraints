#!/bin/bash

currentDIR=$1
codeLoc=$2
workLoc=$3
boxLength=$4
zStart=$5
zStep=$6
numZsteps=$7

pushd "$workLoc/$currentDIR/Boxes" > /dev/null

avTbFname="avTb_$currentDIR".dat
if [ ! -f $avTbFname ]; then
    touch $avTbFname
fi
nfFname="nf_$currentDIR".dat

if [ ! -f $nfFname ]; then
    touch $nfFname
fi
weightedxeFname="weightedxe_$currentDIR".dat
line_weightedxe_Fname="line_weightedxe_$currentDIR".dat
if [ ! -f $weightedxeFname ]; then
    touch $weightedxeFname
fi
if [ -f $line_weightedxe_Fname ]; then
    touch $line_weightedxe_Fname
fi

for (( i=0 ; $i<$numZsteps ; i=$i+1 )) ; do
        zValue=`echo "scale=20; $zStart + $i * $zStep" | bc -l`
	zString="`printf "%06.2f\n" $zValue`"
        echo "i=$i ; zValue is $zValue ; zString is $zString"

        ls delta_T_v3_no_halos_z*"$zString"* > tempFiles.tmp
        currentFname=`head -1 tempFiles.tmp`
        echo "currentFname is $currentFname"
        ls updated_smoothed_deltax_z"$zString"* > tempFiles.tmp
        currentDensityFname=`head -1 tempFiles.tmp`
        echo "currentDensityFname is $currentDensityFname"
        ls xH_nohalos_z"$zString"* > tempFiles.tmp
        currentNeutralFracFname=`head -1 tempFiles.tmp`
        echo "currentNeutralFracFname is $currentNeutralFracFname"
        rm tempFiles.tmp

        cp $currentFname temp_box.dat
        cp $currentDensityFname density_box.dat
        cp $currentNeutralFracFname neutralFrac_box.dat

        echo "temp_box.dat $boxLength temp_avTb.dat" > args.dat
	$codeLoc/./extractSpatialAv.x
        avTb=`cat temp_avTb.dat`
        echo "$zValue $avTb" >> $avTbFname
	rm -f temp_avTb.dat temp_box.dat
        #mv args.dat argsa_"$zString".dat
	rm args.dat

        echo "neutralFrac_box.dat $boxLength temp_nf.dat" > args.dat
        $codeLoc/./extractSpatialAv.x
        nf=`cat temp_nf.dat`
        echo "$zValue $nf" >> $nfFname
        rm -f temp_nf.dat
        rm args.dat

        echo "density_box.dat neutralFrac_box.dat $boxLength temp_weightedxe.dat" > args.dat
        $codeLoc/./compute_weightedIonizedFrac_fullBox.x
	weightedxe=`cat temp_weightedxe.dat`
        echo "$zValue $weightedxe" >> $weightedxeFname

        rm -f temp_weightedxe.dat density_box.dat neutralFrac_box.dat
	#mv args.dat argsb_"$zString".dat  
        rm args.dat
	cp $avTbFname "$workLoc/ionHist/$avTbFname"
        cp $nfFname "$workLoc/ionHist/$nfFname"
        cp $weightedxeFname "$workLoc/ionHist/$weightedxeFname"

        #rm $currentFname $currentDensityFname $currentNeutralFracFname 
	#rm updated_vx_z"$zString"* updated_vy_z"$zString"* updated_vz_z"$zString"*  
done

popd > /dev/null