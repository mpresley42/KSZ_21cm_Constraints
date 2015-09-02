#!/bin/bash
#PBS -S /bin/bash
#PBS -N Planck_TT_TE_EE_lowP_lensing_ext_varyingTvir
#PBS -j eo
#PBS -l mppwidth=72,walltime=05:00:00
#PBS -q premium
#PBS -A m1871

# This is the monster script that varies EVERYTHING!
# It varies:
# i) Cosmological parameters (and reruns density boxes)
# ii) Astrophysical parameters
# iii) Random seed for cosmo parameters if desired.

module swap PrgEnv-intel PrgEnv-gnu
module load python/2.7.3
module load numpy/1.7.1
module load fftw
module load gsl
module load taskfarmer

numNodes=3
numSamples=$numNodes
toDoList="varyTvir.dat"
scriptDIR=$PBS_O_WORKDIR
workLoc="/global/scratch2/sd/acliu/tauConstraint/Planck_TT_TE_EE_lowP_lensing_ext_varyingTvir" #"$scriptDIR/.."
codeLoc="$scriptDIR/../code"
simCodeLoc="/global/u1/a/acliu/software/21cmFAST_edison" #"/Users/Adrian/Codes/21cmFAST"

echo "******************************************************************"
echo "******************************************************************"
echo "******************************************************************"
echo "Here's the script that was run"
echo "$(cat $0)"
echo "******************************************************************"
echo "******************************************************************"
echo "******************************************************************"

boxLength=300

# zStart=30.0 # inclusive
# zStep=-0.5
# numCoarseSteps=5
# numFineSteps=11
zStart=13.0 # inclusive
zStep=-0.1
numCoarseSteps=4
numFineSteps=18

#Tvir_fid=55000.0
#Squiggly_fid=35.0
#Rmfp_fid=25.0

### Prep the analysis directories
if [ ! -d $workLoc ]; then
    mkdir $workLoc
fi
pushd $workLoc

totalNumTodo=`cat $toDoList | wc -l`
### This portion also sets up the necessary directories
tf -t $totalNumTodo -n $numNodes -e serial.err -o serial.out $scriptDIR/./z0_init_compilation_tf.sh $toDoList $simCodeLoc

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
    aprun -n 1 -N 1 -d 24 ./init >& "$currentDIR".log &
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

    tf -t $totalNumTodo -n $numNodes -e serial.err -o serial.out $scriptDIR/./evolveDensity_compilation_tf.sh $toDoList $zTopOfChunk $zStep $numFineSteps
    
    mv serial.out $scriptDIR/serial_evolveDensity.out
    mv tf.log $scriptDIR/tf_evolveDensity.log

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
        aprun -n 1 -N 1 -d 24 ./drive_zscroll_noTs_evolveDensity >& "$currentDIR".log &
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
        aprun -n 1 -N 1 -d 24 ./drive_zscroll_noTs_reion >& "$currentDIR".log &
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
    tf -t $totalNumTodo -n $numNodes -e serial.err -o serial.out $scriptDIR/./extract_allVar_FullBoxIonHistStats_tf.sh $toDoList $codeLoc $workLoc $boxLength $zTopOfChunk $zStep $numFineSteps

    mv serial.out $scriptDIR/serial_extract_IonHiststats.out
    mv tf.log $scriptDIR/tf_extract_IonHiststats.log
    popd > /dev/null
done

popd > /dev/null
