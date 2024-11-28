!/bin/bash
#SBATCH --job-name=pyscenic_job       # Job name
#SBATCH --output=pyscenic_%j.out       # Output file name with job ID
#SBATCH --error=pyscenic_%j.err        # Error file name with job ID
#SBATCH --ntasks=1                     # Number of tasks (processes)
#SBATCH --cpus-per-task=20             # Number of CPU cores per task
#SBATCH --mem=320G 
#SBATCH --time=1-00:00:00          # Time limit (HH:MM:SS)
#SBATCH --mail-type=END,FAIL

echo "Running on hostname `hostname`"
echo "Starting run at: `date`"

# Load necessary modules
module load apptainer/1.2.4

INPUT_LOOM="./pySCENIC_auxiliary/20240923_cells_combined_scenic.loom"
ALL_TFS="./pySCENIC_auxiliary/allTFs_hg38.txt"
OUTPUT_REG="reg.csv"
OUTPUT_ADJ_AGG="consensus_grn_results.csv"
RANKINGS="./pySCENIC_auxiliary/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather ./pySCENIC_auxiliary/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
ANNOTATIONS="./pySCENIC_auxiliary/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
OUTPUT_LOOM="pyscenic_output.loom"
OUTPUT_DIR="./pySCENIC_auxiliary/grn_results"

for i in {1..100}
do
OUTPUT_ADJ="./pySCENIC_auxiliary/grn_results/grn_results_AG_run_$i.csv"
       apptainer exec pyscenic.sif pyscenic grn $INPUT_LOOM $ALL_TFS -o $OUTPUT_ADJ --num_workers 20
done

#run pyscenic_aggregrate_grn.py after the grn step, before the ctx step

apptainer exec pyscenic.sif pyscenic ctx $OUTPUT_ADJ_AGG $RANKINGS \
--annotations_fname $ANNOTATIONS \
--expression_mtx_fname $INPUT_LOOM \
--output $OUTPUT_REG \
--mask_dropouts \
--num_workers 20

apptainer exec pyscenic.sif pyscenic aucell \
$INPUT_LOOM \
$OUTPUT_REG \
--output $OUTPUT_LOOM \
--num_workers 20

echo " Job finished with exit code $? at: `date`"