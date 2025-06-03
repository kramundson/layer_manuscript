#!/bin/bash

# Author: Katharina J. hoff
# Contact: katharina.hoff@uni-greifswald.de
# Date: Jan 12th 2023

# Copy this script into the folder where you want to execute it, e.g.:
# singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test3.sh .
# Then run "bash test3.sh".

THREADS=36

# Check whether braker3.sif is available

if [[ -z "${BRAKER_SIF}" ]]; then
    echo ""
    echo "Variable BRAKER_SIF is undefined."
    echo "First, build the sif-file with \"singularity build braker3.sif docker://teambraker/braker3:latest\""
    echo ""
    echo "After building, export the BRAKER_SIF environment variable on the host as follows:"
    echo ""
    echo "export BRAKER_SIF=\$PWD/braker3.sif"
    echo ""
    echo "You will have to modify the export statement if braker3.sif does not reside in \$PWD."
    echo ""
    exit 1
fi

# Check whether singularity exists

if ! command -v singularity &> /dev/null
then
    echo "Singularity could not be found."
    echo "On some HPC systems you can load it with \"module load singularity\"."
    echo "If that fails, please install singularity."
    echo "Possibly you misunderstood how to run this script. Before running it, please copy it to the directory where you want to execute it by e.g.:"
    echo "singularity exec -B \$PWD:\$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test3.sh ."
    echo "Then execute on the host with \"bash test3.sh\"".
    exit 1
fi

# remove output directory if it already exists
wd=RP_braker3

if [ -d $wd ]; then
    rm -r $wd
fi

singularity exec -B ${PWD}:${PWD} -B /share/comailab/kramundson/miniconda3/bin/:/share/comailab/kramundson/miniconda3/bin/ ${BRAKER_SIF} braker.pl \
        --genome=RP_tetra_homcov140_rmout/RP_tetra_homcov140_purged.p_ctg.fa.masked \
        --prot_seq=ref/DM_M82_ODB10.fa \
        --softmasking \
        --workingdir=${wd} \
        --threads=$THREADS \
        --bam=data/merged/C88_aboveground_stolons.bam  \
        --bam=data/merged/C88_in_vitro_roots.bam \
        --bam=data/merged/C88_stamens.bam \
        --bam=data/merged/C88_apical_buds.bam \
        --bam=data/merged/C88_in_vitro_shoots.bam \
        --bam=data/merged/C88_stems.bam \
        --bam=data/merged/C88_big_tubers.bam \
        --bam=data/merged/C88_mediate_tubers.bam \
        --bam=data/merged/C88_underground_stolons.bam \
        --bam=data/merged/C88_carpels.bam \
        --bam=data/merged/C88_petals.bam \
        --bam=data/merged/C88_very_small_tubers.bam \
        --bam=data/merged/C88_expanded_leaves.bam \
        --bam=data/merged/C88_petiole_from_expanded_leaves.bam \
        --bam=data/merged/C88_whole_mature_flowers.bam \
        --bam=data/merged/C88_flower_buds.bam \
        --bam=data/merged/C88_sepals.bam \
        --bam=data/merged/C88_young_leaves.bam \
        --bam=data/merged/PI310467_leaf_1.bam \
        --bam=data/merged/PI310467_leaf_2.bam \
        --bam=data/merged/PI310467_leaf_3.bam \
        --bam=data/merged/PI310467_leaf_4.bam \
        --bam=data/merged/PI310467_leaf_5.bam \
        --bam=data/merged/PI310467_leaf_6.bam \
        --bam=data/merged/MF93_leaf_1.bam \
        --bam=data/merged/MF93_leaf_2.bam \
        --bam=data/merged/MF93_leaf_3.bam \
        --bam=data/merged/MF93_leaf_4.bam \
        --bam=data/merged/MF93_leaf_5.bam \
        --bam=data/merged/MF93_leaf_6.bam \
        --AUGUSTUS_CONFIG_PATH=/scratch/kramundson/RP_anno/config/ \
        --JAVA_PATH=/share/comailab/kramundson/miniconda3/bin/ \
        --species=Sp_6 \
        --useexisting

	    # Important: the options --gm_max_intergenic 10000 --skipOptimize should never be applied to a real life run!!!                                   
            # They were only introduced to speed up the test. Please delete them from the script if you use it for real data analysis. 
