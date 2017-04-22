#!/bin/bash 

# set -eu

if [ x${1} = x ]; then
    echo "please input PDB ID. stop."
    exit 1
fi

PDBID=${1}
MODELID=1

if [ ${AMBERHOME} = x ]; then
    echo "please set AMBERHOME variable. stop."
    exit 1
fi
export PATH=${AMBERHOME}/bin:${PATH}

# if you would like to use parallel sander,
# please set DO_PARALLEL environment variable, such as
# DO_PARALLEL="${HOME}/local/gnu/openmpi/bin/mpiexec -n 4 "

RELAX_DEBUG="--debug"

# model select --------------------------------------------------------
run_reduce()
{
    echo "run reduce for protonation..."
    PDBID=${1}
    if [ x${PDBID} = x ]; then
        echo "please input PDBID. stop."
        exit 1
    fi

    PDBFILE=${PDBID}.pdb
    PDBFILE_ORIG=${PDBID}.orig.pdb

    cp ${PDBFILE} ${PDBFILE_ORIG}
    echo "original PDB file is saved as ${PDBFILE_ORIG}"

    reduce ${PDBFILE_ORIG} > ${PDBFILE}
    echo "protonated PDB file is saved as ${PDBFILE}"
}

# ----------------------------------------------------------------------
run_leap()
{
    STEP=$1
    ${AMBERHOME}/bin/tleap -s -f leap_md${STEP}.in
}


run_sander()
{
    STEP=$1

    SANDER="${AMBERHOME}/bin/sander.MPI"
    if [ -z "${DO_PARALLEL}" ]; then
        SANDER="${AMBERHOME}/bin/sander"
    fi

    ${DO_PARALLEL} ${SANDER} -O \
        -i md${STEP}.in \
        -o md${STEP}.out \
        -c md${STEP}.inpcrd \
        -p md${STEP}.prmtop \
        -r md${STEP}.restrt \
        -x md${STEP}.mdcrd \
        -v md${STEP}.mdvel \
        -e md${STEP}.mden
}

make_pdb()
{
    ambpdb -p md${STEP}.prmtop -c md${STEP}.restrt -aatm  > md${STEP}_after.pdb
}

step1()
{
    echo ">>>> step1"
    ${PDF_HOME}/bin/relax_protein.py ${RELAX_DEBUG} ${PDBID}.noWAT.brd 1
    echo

    echo "do sander..."
    run_sander 1
    make_pdb 1
    ${PDF_HOME}/bin/pdb2brd.py md1_after.pdb md1_after.brd
    echo "<<<< step1 done."
    echo
}

step2()
{
    echo ">>>> step2"
    ${PDF_HOME}/bin/relax_protein.py  ${RELAX_DEBUG} md1_after.brd 2
    echo

    echo "do sander..."
    run_sander 2
    make_pdb 2
    ${PDF_HOME}/bin/pdb2brd.py md2_after.pdb md2_after.brd
    echo "<<<< step2 done."
    echo
}

step3()
{
    echo ">>>> step3"
    ${PDF_HOME}/bin/relax_protein.py  ${RELAX_DEBUG} md2_after.brd 3
    echo

    echo "do sander..."
    run_sander 3
    make_pdb 3
    ${PDF_HOME}/bin/pdb2brd.py md3_after.pdb md3_after.brd
    echo "<<<< step3 done."
    echo
}

# main
echo "protonate..."
run_reduce ${PDBID}
echo "done."
echo

echo "save brd file..."
${PDF_HOME}/bin/pdb2brd.py --model 1 ${PDBID}.pdb ${PDBID}.brd
echo "done."
echo

echo "remove WAT..."
${PDF_HOME}/bin/remove_wat.py -o ${PDBID}.noWAT.brd ${PDBID}.brd
echo "done."
echo

step1
step2
step3

