#!/bin/bash

toDoList="varyTvir.dat"
scriptDIR="$PWD"
workLoc="/Users/mpresley/Data/21cmFAST/test" # location of output. **Needs to have varyTvir.dat file!!**
simCodeLoc="/Users/mpresley/Soft/21cmFAST"
codeLoc="$simCodeLoc/Programs" # location of fortran codes

echo "******************************************************************"
echo "******************************************************************"
echo "Here's the script that was run"
echo "$(cat $0)"
echo "******************************************************************"
echo "******************************************************************"

# Set box size parameters
boxLength=512   # number of cells for the low-res box               
boxLenHiRes=1024 # number of cells for the high-res box (must be int mult of boxLength)
boxLenMpc=1024   # this should be changed by same factor as HII_DIM                    

sed "33s/200/$boxLenMpc/" $simCodeLoc/Parameter_files/INIT_PARAMS_original.H > $simCodeLoc/Parameter_files/INIT_PARAMS_1.H
sed "34s/900/$boxLenHiRes/" $simCodeLoc/Parameter_files/INIT_PARAMS_1.H > $simCodeLoc/Parameter_files/INIT_PARAMS_2.H
sed "35s/300/$boxLength/" $simCodeLoc/Parameter_files/INIT_PARAMS_2.H > $simCodeLoc/Parameter_files/INIT_PARAMS.H

### Set parameters for z loop
zStart=30.0 # inclusive                        
zStep=-0.5
numCoarseSteps=5 # number of times it deletes files                                      
numFineSteps=10 # total num of redshifts = coarse * fine 

### Prep the analysis directories
if [ ! -d $workLoc ]; then
    mkdir $workLoc
fi
cd $workLoc
echo $PWD

toDoList="$workLoc/$toDoList"
echo $toDoList

### This portion also sets up the necessary directories
$scriptDIR/z0_init.sh $toDoList $simCodeLoc
echo "finished z0_init.sh"

### Set up current directory
currentLine=`cat $toDoList | head -1 | tail -1`
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

echo "Setting up initial density field for (randSeed, sigma8, h, Omm, Omb, ns, Rmfp, Tvir, Squiggly, lnAs) = ($currentRandSeed, $currentSigma8, $currenthlittle, $currentOmegam, $currentOmegab, $currentns, $currentRmfp, $currentTvir, $currentSquiggly, $currentlnAs)..."

currentDIR="RandSeed_$currentRandSeed""_Sigma8_$currentSigma8""_h_$currenthlittle""_Omm_$currentOmegam""_Omb_$currentOmegab""_ns_$currentns""_Rmfp_$currentRmfp""_Tvir_$currentTvir""_Squiggly_$currentSquiggly""_lnAs_$currentlnAs"

cd $currentDIR/Programs
echo $PWD
currentCodeLoc="$PWD"
echo $currentCodeLoc

### Compute the density fields at z = 0 once and for all
echo "...setting the actual program going..."
currentTime=`date`
echo "Starting this run at $currentTime"
./init 1>$workLoc/init.out 2>$workLoc/init.err
echo "finished ./init"

### Actually begin z loop
for (( i=0 ; $i<$numCoarseSteps ; i=$i+1 )) ; do
    echo "I'm on coarse step $i..."
    zTopOfChunk=`echo "scale=20; $zStart + ($numFineSteps * $i) * $zStep" | bc -l`
    zBottomOfChunk=`echo "scale=20; $zStart + ($numFineSteps * $i + $numFineSteps - 1) * $zStep" | bc -l`

    echo "Running redshift range $zTopOfChunk to $zBottomOfChunk inclusive"
    ### Now compile all the codes for the relevant redshifts

    $scriptDIR/evolveDensity.sh $currentDIR/Programs $zTopOfChunk $zStep $numFineSteps 1>$workLoc/evolveDensity_$i.out 2>$workLoc/evolveDensity_$i.out
    echo "finished evolveDensity.sh"

    ### Evolve the density field 
    echo "...setting the actual program going..."
    currentTime=`date`
    echo "Starting this run at $currentTime"
    $currentCodeLoc/drive_zscroll_noTs_evolveDensity 1>$workLoc/drive_evolveDensity_$i.out 2>$workLoc/drive_evolveDensity_$i.err 
    echo "finished drive_zscroll_noTs_evolveDensity"

    ### Run the reionization simulations 
    echo $PWD
    echo "...setting the actual program going..."

    currentTime=`date`
    echo "Starting this run at $currentTime"
    $currentCodeLoc/drive_zscroll_noTs_reion 1>$workLoc/drive_reion_$i.out 2>$workLoc/drive_reion_$i.err

    ### Extract ionization history stats and delete files to free up space 
    cd $workLoc
    echo $PWD
    if [ ! -d ionHist ]; then
        mkdir ionHist
    fi
    $scriptDIR/extract_FullBoxIonHistStats.sh $currentDIR $codeLoc $workLoc $boxLength $zTopOfChunk $zStep $numFineSteps 1>$workLoc/extract_IonHistStats_$i.out 2>$workLoc/extract_IonHistStats_$i.err
    echo "finished extract_allVar_FullBoxIonHistStats.sh"

done

# combine the boxes                                                                                
# echo "hi"
# echo $PWD 
# cd "$workLoc/$currentDIR/Boxes"
# echo $PWD

# echo "Combine All the Boxes!"

# ls updated*deltax* >  deltax_list.txt
# $currentCodeLoc/redshift_interpolate_boxes 0 deltax_list.txt

# ls updated_vx* > vx_list.txt
# $currentCodeLoc/redshift_interpolate_boxes 0 vx_list.txt

# ls updated_vy* > vy_list.txt
# $currentCodeLoc/redshift_interpolate_boxes 0 vy_list.txt

# ls updated_vz* > vz_list.txt
# $currentCodeLoc/redshift_interpolate_boxes 0 vz_list.txt

# ls xH_nohalos_* > xH_list.txt
# $currentCodeLoc/redshift_interpolate_boxes 0 xH_list.txt

# ls delta_T_* > dT_list.txt
# $currentCodeLoc/redshift_interpolate_boxes 0 dT_list.txt

# scrape 21cm data 
# cp $scriptDir/scrape_21cm.py ./
# python ./scrape_21cm.py		
