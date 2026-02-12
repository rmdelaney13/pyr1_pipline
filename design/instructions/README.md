How to use ligand design pipeline:

#generating table
python /projects/ryde3462/software/ligand_alignment/create_table.py /projects/ryde3462/xanthurenic_acid/config_multiple.txt

# fill in atoms that you want to align
#decide num passes in config_multiple.txt

#dock_glycine_shaved
python /projects/ryde3462/software/ligand_alignment/grade_conformers_glycine_shaved_docking_multiple.py /projects/ryde3462/bile_acids/CA/config_multiple.txt > local_run.log 2>&1



#run docking in parallel using slurm script
#must change location of parent directiory in conformer_check_submit.sh
#set num passes to 1 in config_multiple.txt
#use multiple lines in ligand.csv to perform multiple passes for one job
sbatch conformer_check_submit.sh

#!/bin/bash
# This script replaces all occurrences of "WAT" with "TP3" in each .pdb file in the current directory so that ligandmpnn recognizes water .

for pdb in *.pdb; do
    echo "Processing ${pdb}..."
    sed -i 's/WAT/TP3/g' "${pdb}"
done




##

python /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/instructions/aggregate_scores.py /scratch/alpine/ryde3462/xan_design/relax_2 --output /scratch/alpine/ryde3462/xan_design/relax_2/relax_2.csv


python /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/instructions/relax_2_filter__allpolar_unsats.py \
  /scratch/alpine/ryde3462/CA_design_flipped/relax_2/relax_2.csv \
  /scratch/alpine/ryde3462/CA_design_flipped/relax_2 \
  /scratch/alpine/ryde3462/CA_design_flipped/filtered_2 \
  --target_n 2000 \
  --max_unsat 1 \
  --max_per_parent 30 \
  --ignore_charge


python /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/instructions/relax_2_filter__allpolar_unsats_xan.py \
  /scratch/alpine/ryde3462/xan_design/relax_2/relax_2.csv \
  /scratch/alpine/ryde3462/xan_design/relax_2 \
  /scratch/alpine/ryde3462/xan_design/filtered_2 \
  --target_n 1000 \
  --max_interface_unsat 4 \
  --max_lig_unsat 2 \
  --max_per_parent 20



python /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/split_and_mutate_to_fasta.py \
/scratch/alpine/ryde3462/xan_design/filtered_2 \
/scratch/alpine/ryde3462/xan_design/filtered_2/filtered.fasta

python /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/make_af3_jsons.py \
--template /projects/ryde3462/pyr_template/LCA/pyr1_xan_closed_template.json \
--fasta /scratch/alpine/ryde3462/xan_design/filtered_2/filtered.fasta \
--outdir /scratch/alpine/ryde3462/xan_design/af3/binary_input_2



python /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/aggregate_json_batches.py \
/scratch/alpine/ryde3462/xan_design/af3/ternary_input_2/ \
--limit 100

sbatch /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/alphafold3_alpine_gpu_dir_no_msas_array.sh




python /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/instructions/binary_analysis.py \
--inference_dir /scratch/alpine/ryde3462/xan_design/af3/binary_output_2 \
--ref_model /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/pyr_win_template.cif \
--output_csv /scratch/alpine/ryde3462/xan_design/af3/binary_output_2/xan_binary_full_metrics.csv 



#binary sort
python /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/instructions/sort_binary_models.py \
  --csv /scratch/alpine/ryde3462/xan_design/af3/binary_output_2/xan_binary_full_metrics.csv \
  --run_dir /scratch/alpine/ryde3462/xan_design/af3/binary_output_2 \
  --plddt_min 60.0 \
  --iptm_min 0.65 \
  --max_dist 3 \
  --out_dir /scratch/alpine/ryde3462/xan_design/af3/binary_output_2/top_binary \
  --out_csv /scratch/alpine/ryde3462/xan_design/af3/binary_output_2/top_binary/filtered_metrics.csv \
  --copy_mode symlink \
  --overwrite

python /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/extract_chain_sequence_from_cifs.py \
--in_dir /scratch/alpine/ryde3462/xan_design/af3/binary_output_2/top_binary \
--pattern "*_model.cif" \
--chain A \
--out_fasta /scratch/alpine/ryde3462/xan_design/af3/binary_output_2/top_binary/xan.fasta \
--out_csv /scratch/alpine/ryde3462/xan_design/af3/binary_output_2/top_binary/binary_chainA.csv

python /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/make_af3_jsons.py \
--template /projects/ryde3462/pyr_template/LCA/pyr1_ternary_xan_closed_template.json \
--fasta /scratch/alpine/ryde3462/xan_design/af3/binary_output_2/top_binary/xan.fasta \
--outdir /scratch/alpine/ryde3462/xan_design/af3/ternary_input_2



python /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/instructions/ternary_analysis.py \
--inference_dir /scratch/alpine/ryde3462/xan_design/af3/ternary_output_2 \
--ref_model /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/pyr_win_template.cif \
--output_csv /scratch/alpine/ryde3462/xan_design/af3/ternary_output_2/xan_ternary_full_metrics.csv 



python /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/instructions/aggregate_all_metrics.py \
--ternary_csv /scratch/alpine/ryde3462/xan_design/af3/ternary_output_2/xan_ternary_full_metrics.csv \
--binary_csv  /scratch/alpine/ryde3462/xan_design/af3/binary_output_2/xan_binary_full_metrics.csv \
--ternary_dir /scratch/alpine/ryde3462/xan_design/af3/ternary_output_2 \
--binary_dir  /scratch/alpine/ryde3462/xan_design/af3/binary_output_2 \
--output_csv  /scratch/alpine/ryde3462/xan_design/af3/master_summary_2.csv


python /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/instructions/sort_by_rmsd.py \
--input_csv /scratch/alpine/ryde3462/xan_design/af3/master_summary_2.csv \
--ternary_dir /scratch/alpine/ryde3462/xan_design/af3/ternary_output_2 \
--binary_dir /scratch/alpine/ryde3462/xan_design/af3/binary_output_2 \
--output_dir /scratch/alpine/ryde3462/xan_design/af3/low_rmsd_results_1_05 \
--threshold 0.7



##### extra commands that are helpful for analysis

python /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/logo_from_fasta.py \
/scratch/alpine/ryde3462/bile_seqs_20260114/final_fastas/ALL_dedup_pos59_VLGA_full.fasta \
--positions 59 81 83 92 94 108 110 117 120 122 141 159 160 163 164 167 \
--title "Final_Order (N=2146)" \
--out /scratch/alpine/ryde3462/bile_seqs_20260114/final_fastas/final_order.png



python /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/split_and_mutate_to_fasta_no_insert.py \
/scratch/alpine/ryde3462/CDCA_design_normal/af3/binary_outputs_2/top_binary_pdb \
/scratch/alpine/ryde3462/CDCA_design_normal/af3/binary_outputs_2/top_binary_pdb/filtered.fasta




python /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/cif_to_pdb.py \
--in_dir  /scratch/alpine/ryde3462/DCA_design_flipped/af3/binary_outputs/top_binary/  \
--out_dir /scratch/alpine/ryde3462/DCA_design_flipped/af3/binary_outputs/top_binary_pdb
--pattern *.cif


python /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/cifs_to_fasta.py \
  --in_dir /scratch/alpine/ryde3462/instance_oligo/LCA/top_binary \
  --out_file top_binary.fasta \
  --chain A



python /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/instructions/take_top_plddt_from_top_binary.py \
    --csv /scratch/alpine/ryde3462/kyna_design/af3/binary_output_1/top_binary/filtered_metrics.csv \
    --models_dir /scratch/alpine/ryde3462/kyna_design/af3/binary_output_1/top_binary \
    --out_dir /scratch/alpine/ryde3462/kyna_design/af3/binary_output_1/ternary_top20_plddt \
    --out_csv /scratch/alpine/ryde3462/kyna_design/af3/binary_output_1/ternary_top20_plddt/top20_plddt_metrics.csv \
    --n 20 \
    --id_col target \
    --plddt_col ligand_plddt \
    --copy_mode hardlink \
    --overwrite


python /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/instructions/rank_pair_ternary_binary.py \
  --ternary_csv /scratch/alpine/ryde3462/kyna_design/af3/ternary_output_2/kyna_ternary_full_metrics.csv \
  --binary_csv  /scratch/alpine/ryde3462/kyna_design/af3/binary_output_2/kyna_binary_full_metrics.csv \
  --ternary_dir /scratch/alpine/ryde3462/kyna_design/af3/ternary_output_2 \
  --binary_dir  /scratch/alpine/ryde3462/kyna_design/af3/binary_output_2 \
  --out_dir     /scratch/alpine/ryde3462/kyna_design/af3/paired_ranked_2 \
  --top_n 6 \
  --copy_mode hardlink


python /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/instructions/rank_pair_ternary_binary_diverse.py \
  --ternary_csv /scratch/alpine/ryde3462/CDCA_design_normal/af3/ternary_outputs_2/CDCA_ternary_full_metrics.csv \
  --binary_csv  /scratch/alpine/ryde3462/CDCA_design_normal/af3/binary_outputs_2/CDCA_full_metrics.csv \
  --ternary_dir /scratch/alpine/ryde3462/CDCA_design_normal/af3/ternary_outputs_2 \
  --binary_dir  /scratch/alpine/ryde3462/CDCA_design_normal/af3/binary_outputs_2 \
  --out_dir /scratch/alpine/ryde3462/CDCA_design_normal/af3/paired_ranked \
  --preselect_n_plddt 100 \
  --final_n_diverse 10 \
  --cluster_identity 0.95 \
  --diverse_dirname top_cifs_diverse \
  --copy_mode hardlink



python /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/instructions/make_logo_from_cif_dir.py \
  /scratch/alpine/ryde3462/kyna_design/af3/low_rmsd_results_1_2/models \
  --chain A \
  --positions 59 81 83 92 94 108 110 117 120 122 141 159 160 163 164 167 \
  --title "AF3_kyna" \
  --out /scratch/alpine/ryde3462/kyna_design/af3/low_rmsd_results_1_2/models/af3_binary_kyna.png



#cif to fasta
python /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/extract_chain_sequence_from_cifs.py \
--in_dir /scratch/alpine/ryde3462/kyna_design/af3/binary_output_1/top_binary \
--pattern "*_model.cif" \
--chain A \
--out_fasta /scratch/alpine/ryde3462/kyna_design/af3/binary_output_1/top_binary/kyna.fasta \
--out_csv /scratch/alpine/ryde3462/kyna_design/af3/binary_output_1/top_binary/binary_chainA.csv



#######################################################################################################################################################################################33


#####Aggregating AF# seqs from mpnn design steP
python /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/merge_chainA_fasta.py \
--root_dir /scratch/alpine/ryde3462/CA_design/mpnn_4_af3 \
--output_fasta /scratch/alpine/ryde3462/CA_design/mpnn_4_af3/all_chainA_merged.fa


python /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/merge_chainA_fasta_new_seqs.py \
  --root_dir /scratch/alpine/ryde3462/CA_design_final/mpnn_1_normal \
  --output_fasta /scratch/alpine/ryde3462/CA_design_final/mpnn_1_normal_merged.fa



python /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/make_af3_jsons.py \
--template /projects/ryde3462/pyr_template/LCA/pyr1_ternary_kyna_closed_template.json \
--fasta /scratch/alpine/ryde3462/kyna_design/af3/binary_output_2/top_binary/kyna.fasta \
--outdir /scratch/alpine/ryde3462/kyna_design/af3/ternary_intput_2/

python /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/ternary_analysis.py \
  --dir_ternary_outputs /scratch/alpine/ryde3462/instance_oligo/CDCA_testing/ternary_outputs \
  --output_csv /scratch/alpine/ryde3462/instance_oligo/CDCA_testing/ternary_outputs/output.csv

python /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/select_for_md.py \
    --results_csv /scratch/alpine/ryde3462/LCA_1_no_hab/af3_inputs_1/output/final_comprehensive_results.csv \
    --source_cif_dir /scratch/alpine/ryde3462/LCA_1_no_hab/af3_inputs_1/output/filtered_models_60_3.4_5.5_hab/msas/output/collected_models \
    --dest_dir /scratch/alpine/ryde3462/LCA_1_no_hab/top_20_for_md






##