#!/bin/bash
#
#SBATCH -J def_dedup
#SBATCH -o /proj/omics4tb2/Global_Search/Pilot_Fail_New/scripts/logs/%x_log.out
#SBATCH -e /proj/omics4tb2/Global_Search/Pilot_Fail_New/scripts/logs/%x_log.out

genome_dir="/proj/omics4tb2/Global_Search/reference_genomes/past_smic"
data_root="/proj/omics4tb2/Global_Search/Pilot_Fail_New/raw_data"
data_folder="R1"
out_folder="/proj/omics4tb2/Global_Search/Pilot_Fail_New/results"
outFilterMismatchNmax=10
outFilterMismatchNoverLmax=0.01
outFilterScoreMinOverLread=0.05
outFilterMatchNmin=10
dedup_prefix="dedup"
if [ $dedup_prefix == "dedup" ]
then
    dedup_cmd="--dedup"
else
    dedup_cmd=""
fi    

star_prefix="star_${outFilterMismatchNmax}_${outFilterMismatchNoverLmax}_${outFilterScoreMinOverLread}_${outFilterMatchNmin}_${dedup_prefix}"
salmon_prefix="salmon_${outFilterMismatchNmax}_${outFilterMismatchNoverLmax}_${outFilterScoreMinOverLread}_${outFilterMatchNmin}_${dedup_prefix}"

echo $star_prefix
echo $salmon_prefix
#python run_STAR_Salmon.py $genome_dir $data_root $data_folder $out_folder
python run_STAR_Salmon.py $genome_dir $data_root $data_folder $out_folder --outFilterMismatchNmax $outFilterMismatchNmax --outFilterMismatchNoverLmax $outFilterMismatchNoverLmax --outFilterScoreMinOverLread $outFilterScoreMinOverLread --outFilterMatchNmin $outFilterMatchNmin $dedup_cmd --starPrefix $star_prefix --salmonPrefix $salmon_prefix


