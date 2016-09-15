#!/bin/bash
#SBATCH
#SBATCH -o test_randomDefaultSeed.out
#SBATCH -p campus
#SBATCH -c 4
phydms_testdivpressure defaultSeed/ H1_HA_human_coding_IVR_wsn_prep_final_alignment1.fasta HA_renumber.txt HA_divpressure.txt HA_sitediffsel_pos_top5.txt --randomizations 3