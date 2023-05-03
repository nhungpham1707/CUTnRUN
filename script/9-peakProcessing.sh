#!/bin/bash
# this script extract peaks from each condition after peakcalling. Either by intersect, merge or multiinter
# Nhung 11 04 2023

########## find overlap peaks among replicates in the same condition #################
## LUCIFERASE
# get sample paths
declare -a luc_IDs=()
total_sample=${#luciferase[@]}
n=0
for sample_ID in ${luciferase[@]}; do
    n=$((n+1))
    echo "-----------running $n out of $total_sample luc samples---------------------- "

    luc_IDs+=( "${peak_dir}/${sample_ID}/narrow/${sample_ID}_paired_control_peaks.narrowPeak" )

done
# generate overlap report
bedtools intersect -a ${luc_IDs[0]} -b ${luc_IDs[@]:1} -names ${luc_IDs[@]:1} -sorted > ${peak_analysis_dir}/luciferase_intersect.txt

bedtools multiinter -i ${luc_IDs[@]} -header | awk '{if ($4 == 2) print}' > luciferase_multiinter.txt

## TFE3
declare -a tfe3_IDs=()
total_sample=${#tfe3[@]}
n=0
for sample_ID in ${tfe3[@]}; do
    n=$((n+1))
    echo "-----------running $n out of $total_sample tfe3 samples---------------------- "

    luc_IDs+=( "${peak_dir}/${sample_ID}/narrow/${sample_ID}_paired_control_peaks.narrowPeak" )

done
# generate overlap report
bedtools intersect -a ${tfe3_IDs[0]} -b ${tfe3_IDs[@]:1} -names ${tfe3_IDs[@]:1} -sorted > ${peak_analysis_dir}/tfe3_intersect.txt

bedtools multiinter -i ${tfe3_IDs[@]} -header | awk '{if ($4 == 2) print}' > tfe3_multiinter.txt

## FUSION
declare -a fusion_IDs=()
total_sample=${#fusion[@]}
n=0
for sample_ID in ${fusion[@]}; do
    n=$((n+1))
    echo "-----------running $n out of $total_sample fusion samples---------------------- "

    fusion_IDs+=( "${peak_dir}/${sample_ID}/narrow/${sample_ID}_paired_control_peaks.narrowPeak" )

done

bedtools intersect -a ${fusion_IDs[0]} -b ${fusion_IDs[@]:1} -names ${fusion_IDs[@]:1} -sorted > ${peak_analysis_dir}/fusion_intersect.txt

bedtools multiinter -i ${fusion_IDs[@]} -header | awk '{if ($4 == 2) print}' > fusion_multiinter.txt

## HISTONE SAMPLES
declare -a histone_IDs=()
total_sample=${#h3k4me3_all[@]}
n=0
for sample_ID in ${h3k4me3_all[@]}; do
    n=$((n+1))
    echo "-----------running $n out of $total_sample histone samples---------------------- "

    histone_IDs+=( "${peak_no_control_dir}/${sample_ID}/narrow/${sample_ID}_paired_peaks.narrowPeak" )

done

 bedtools intersect -a ${histone_IDs[0]} -b ${histone_IDs[@]:1} -names ${histone_IDs[@]:1} -sorted > ${peak_analysis_dir}/h3k4me3_intersect.txt

bedtools multiinter -i ${histone_IDs[@]} -header > ${peak_analysis_dir}/h3k4me3_multiinter.txt
