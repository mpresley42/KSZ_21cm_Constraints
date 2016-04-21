#!/bin/bash
#SBATCH --job-name=21cmSim
#SBATCH --account=m1871
#SBATCH --partition=debug
#SBATCH --nodes=4
#SBATCH --time=00:05:00
#SBATCH --mail-type==(BEGIN,END,FAIL)
#SBATCH --mail-user==mpresley@berkeley.edu

# This is the monster script that varies EVERYTHING!
# It varies:
# i) Cosmological parameters (and reruns density boxes)
# ii) Astrophysical parameters
# iii) Random seed for cosmo parameters if desired.

# NOTE: This version does not use TaskFarmer and thus
# can only run the first line of the varyTvir.dat file.

module swap PrgEnv-intel PrgEnv-gnu
module load python/2.7.3
module unload numpy/1.9.2
module load numpy/1.7.1
module load fftw
module load gsl

numNodes=4
numSamples=$numNodes
toDoList="varyTvir.dat"
scriptDIR=$PBS_O_WORKDIR # location of this script
workLoc="/scratch1/scratchdirs/mpresley/21cm_FAST_Sims/test6" # location of output. **Needs to have varyTvir.dat file!!!**
simCodeLoc="/global/homes/m/mpresley/soft/21cmFAST"
codeLoc="$simCodeLoc/Programs" # location of fortran codes
runscriptLoc="/global/homes/m/mpresley/KSZ_21cm_Constraints/code/21cmFast/runscripts"

echo "******************************************************************"
echo "******************************************************************"
echo "******************************************************************"
echo "Here's the script that was run"
echo "$(cat $0)"
echo "******************************************************************"
echo "******************************************************************"
echo "******************************************************************"

boxLength=400   # number of cells for the low-res box
boxLenHiRes=1200 # number of cells for the high-res box (must be int mult of boxLength)
boxLenMpc=200   # this should be changed by same factor as HII_DIM

sed "33s/200/$boxLenMpc/" $simCodeLoc/Parameter_files/INIT_PARAMS_original.H > $simCodeLoc/Parameter_files/INIT_PARAMS.H
sed "34s/900/$boxLenHiRes/" $simCodeLoc/Parameter_files/INIT_PARAMS_original.H > $simCodeLoc/Parameter_files/INIT_PARAMS.H
sed "35s/300/$boxLength/" $simCodeLoc/Parameter_files/INIT_PARAMS_original.H > $simCodeLoc/Parameter_files/INIT_PARAMS.H

zStart=30.0 # inclusive
zStep=-0.5
numCoarseSteps=5 # number of times it deletes files
numFineSteps=10 # total num of redshifts = coarse * fine

### Prep the analysis directories
if [ ! -d $workLoc ]; then
    mkdir $workLoc
fi
pushd $workLoc

totalNumTodo=`cat $toDoList | wc -l`
### This portion also sets up the necessary directories
#tf -t $totalNumTodo -n $numNodes -e serial.err -o serial.out $simCodeLoc/Programs/z0_init_compilation_tf.sh $toDoList $simCodeLoc
#echo "tf -t $totalNumTodo -n $numNodes -e serial.err -o serial.out $simCodeLoc/Programs/z0_init_compilation_tf.sh $toDoList $simCodeLoc"
echo "srun -e serial.err -o serial.out $runscriptLoc/z0_init_compilation.sh $toDoList $simCodeLoc"
#srun -e serial.err -o serial.out $runscriptLoc/z0_init_compilation.sh $toDoList $simCodeLoc
echo "finished z0_init_compilation.sh"

### Compute the density fields at z = 0 once and for all
let counter=0
for (( j=1 ; $j<=$totalNumTodo ; j=$j+1 )) ; do
    currentLine=`cat $toDoList | head -$j | tail -1`
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
    pushd $currentDIR/Programs > /dev/null

    echo "...setting the actual program going..."

    currentTime=`date`
    echo "Starting this run at $currentTime"
    export OMP_NUM_THREADS=24
    echo "srun -n 1 -N 1 -c 24 ./init >& "$currentDIR".log & "
    #srun -n 1 -N 1 -c 24 ./init >& "$currentDIR".log &
    echo "finished ./init"
    counter=$(($counter + 1))
    
    echo "Currently there are $counter runs going"
    if [ $counter = $numNodes ]; then
        wait
        let counter=0
    fi

    popd > /dev/null
done
wait

for (( i=0 ; $i<$numCoarseSteps ; i=$i+1 )) ; do
    zTopOfChunk=`echo "scale=20; $zStart + ($numFineSteps * $i) * $zStep" | bc -l`
    zBottomOfChunk=`echo "scale=20; $zStart + ($numFineSteps * $i + $numFineSteps - 1) * $zStep" | bc -l`

    echo "Running redshift range $zTopOfChunk to $zBottomOfChunk inclusive"
    ### Now compile all the codes for the relevant redshifts

    #tf -t $totalNumTodo -n $numNodes -e serial.err -o serial.out $simCodeLoc/Programs/evolveDensity_compilation_tf.sh $toDoList $zTopOfChunk $zStep $numFineSteps
    echo "srun -e serial.err -o serial.out $runscriptLoc/evolveDensity_compilation.sh $toDoList $zTopOfChunk $zStep $numFineSteps"
#    srun -e serial.err -o serial.out $runscriptLoc/evolveDensity_compilation.sh $toDoList $zTopOfChunk $zStep $numFineSteps 
    echo "finished evolveDensity_compilation.sh"

    mv serial.out $workLoc/serial_evolveDensity.out
    #mv tf.log $workLoc/tf_evolveDensity.log

    wait

    ### Evolve the density field
    let counter=0
    for (( j=1 ; $j<=$totalNumTodo ; j=$j+1 )) ; do
        currentLine=`cat $toDoList | head -$j | tail -1`
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

        echo "Evolving density field for for (randSeed, sigma8, h, Omm, Omb, ns, Rmfp, Tvir, Squiggly, lnAs) = ($currentRandSeed, $currentSigma8, $currenthlittle, $currentOmegam, $currentOmegab, $currentns, $currentRmfp, $currentTvir, $currentSquiggly, $currentlnAs)..."
        currentDIR="RandSeed_$currentRandSeed""_Sigma8_$currentSigma8""_h_$currenthlittle""_Omm_$currentOmegam""_Omb_$currentOmegab""_ns_$currentns""_Rmfp_$currentRmfp""_Tvir_$currentTvir""_Squiggly_$currentSquiggly""_lnAs_$currentlnAs"

        pushd $currentDIR/Programs > /dev/null

        echo "...setting the actual program going..."

        currentTime=`date`
        echo "Starting this run at $currentTime"
        export OMP_NUM_THREADS=24
	echo "srun -n 1 -N 1 -c 24 ./drive_zscroll_noTs_evolveDensity >& "$currentDIR".log &"
#        srun -n 1 -N 1 -c 24 ./drive_zscroll_noTs_evolveDensity >& "$currentDIR".log &
	echo "finished drive_zscroll_noTs_evolveDensity"
        counter=$(($counter + 1))
        
        echo "Currently there are $counter runs going"
        if [ $counter = $numNodes ]; then
            #echo "Hey I am waiting now!"
            wait
            let counter=0
        fi

        popd > /dev/null
    done

    wait

    ### Run the reionization simulations
    let counter=0
    for (( j=1 ; $j<=$totalNumTodo ; j=$j+1 )) ; do
        currentLine=`cat $toDoList | head -$j | tail -1`
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

        echo "Simulating reionization for (randSeed, sigma8, h, Omm, Omb, ns, Rmfp, Tvir, Squiggly, lnAs) = ($currentRandSeed, $currentSigma8, $currenthlittle, $currentOmegam, $currentOmegab, $currentns, $currentRmfp, $currentTvir, $currentSquiggly, $currentlnAs)..."
        currentDIR="RandSeed_$currentRandSeed""_Sigma8_$currentSigma8""_h_$currenthlittle""_Omm_$currentOmegam""_Omb_$currentOmegab""_ns_$currentns""_Rmfp_$currentRmfp""_Tvir_$currentTvir""_Squiggly_$currentSquiggly""_lnAs_$currentlnAs"
        pushd $currentDIR/Programs > /dev/null

        echo "...setting the actual program going..."

        currentTime=`date`
        echo "Starting this run at $currentTime"
        export OMP_NUM_THREADS=24
	echo "srun -n 1 -N 1 -c 24 ./drive_zscroll_noTs_reion >& "$currentDIR".log &"
#        srun -n 1 -N 1 -c 24 ./drive_zscroll_noTs_reion >& "$currentDIR".log &
	echo "finished drive_zscroll_noTs_reion"
        counter=$(($counter + 1))
        
        echo "Currently there are $counter runs going"
        if [ $counter = $numNodes ]; then
            wait
            let counter=0
        fi

        popd > /dev/null
    done

    wait

    ### Extract ionization history stats and delete files to free up space
    pushd $workLoc
    if [ ! -d ionHist ]; then
        mkdir ionHist
    fi
    #tf -t $totalNumTodo -n $numNodes -e serial.err -o serial.out $simCodeLoc/Programs/extract_allVar_FullBoxIonHistStats_tf.sh $toDoList $codeLoc $workLoc $boxLength $zTopOfChunk $zStep $numFineSteps
    echo "srun -e serial.err -o serial.out $runscriptLoc/extract_allVar_FullBoxIonHistStats.sh $toDoList $codeLoc $workLoc $boxLength $zTopOfChunk $zStep $numFineSteps"
#    srun -e serial.err -o serial.out $runscriptLoc/extract_allVar_FullBoxIonHistStats.sh $toDoList $codeLoc $workLoc $boxLength $zTopOfChunk $zStep $numFineSteps
    echo "finished extract_allVar_FullBoxIonHistStats.sh"

    mv serial.out $workLoc/serial_extract_IonHiststats.out
    #mv tf.log $workLoc/tf_extract_IonHiststats.log
    popd > /dev/null
done
 
#popd > /dev/null

# combine the boxes
pushd $currentDIR/Boxes > /dev/null

echo "Combine All the Boxes!"
echo $PWD

ls updated*deltax* >  deltax_list.txt 
/global/homes/m/mpresley/soft/21cmFAST/Programs/redshift_interpolate_boxes 0 deltax_list.txt 

ls updated_vx* > vx_list.txt
/global/homes/m/mpresley/soft/21cmFAST/Programs/redshift_interpolate_boxes 0 vx_list.txt

ls updated_vy* > vy_list.txt 
/global/homes/m/mpresley/soft/21cmFAST/Programs/redshift_interpolate_boxes 0 vy_list.txt

ls updated_vz* > vz_list.txt
/global/homes/m/mpresley/soft/21cmFAST/Programs/redshift_interpolate_boxes 0 vz_list.txt

ls xH_nohalos_* > xH_list.txt 
/global/homes/m/mpresley/soft/21cmFAST/Programs/redshift_interpolate_boxes 0 xH_list.txt

ls delta_T_* > dT_list.txt
/global/homes/m/mpresley/soft/21cmFAST/Programs/redshift_interpolate_boxes 0 dT_list.txt

# scrape 21cm data
cp $runscriptLoc/scrape_21cm.py ./ 
python ./scrape_21cm.py
