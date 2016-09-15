#!/bin/bash
#SBATCH
#SBATCH -o test_noRandom.out
phydms_testdivpressure noRandom/ H1_HA_human_coding_IVR_wsn_prep_final_alignment1.fasta HA_renumber.txt HA_divpressure.txt HA_sitediffsel_pos_top5.txt
