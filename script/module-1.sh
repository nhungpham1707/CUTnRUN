# Module 1: Data processing. quality checking, trimming, alignment, replicates correlation
echo "start running module 1- data processing for cut and run analysis at $(date)"
# step 1. sequencing quality check 
echo "------------------step1. running quality check--------------"
. ./1-qualityCheck.sh

# step 2. adapter and bad reads trimming 
echo "-------------------step 2. running trimming-----------------"
. ./2-trimming.sh 

# step 3. alignment- map to hg38 genome 
echo "-------------------step 3. running alignment----------------"
. ./3-alignment.sh

# step 4. filtering: remove duplciates and reads < 20bp
echo "-------------------step 4. running filtering----------------"
. ./4-filtering.sh

# step 5. check sample correlation. To decide if any replicate should be removed
echo "-------------------step 5. running sample correlation--"
. ./5-samplesCorrelation.sh

# step 6. merge and transform bam file to bigwig
echo "-------------------step 6. running transform bam to bigwig--"
. ./6-bam2bigwig.sh

# step 7a. peak calling raw data with control
echo "---------step 7a. running peak calling for raw data tfe3---------"
. ./peakCalling_raw_data.sh "$tfe3C" "${tfe3[@]}" 

echo "---------step 7a. running peak calling for raw data luc---------"
. ./peakCalling_raw_data.sh "$luciferaseC" "${luciferase[@]}" 

echo "---------step 7a. running peak calling for raw data fusion---------"
. ./peakCalling_raw_data.sh "$fusionC" "${fusion[@]}" 

# step 7b. peak calling raw data without control
echo "---------step 7b. running peak calling for raw data without control ---------"
. ./peakCalling_no_control.sh

# step 8. calculate FRiP
echo "---------step 8. running FRiP calculation---------"
. ./Calculate_FRiP.sh  

# step 9. merge replicates # fix this step 
. ./merge_replicate.sh "tfe3_test2" "${tfe3[@]}"

# . ./merge_replicate.sh "luc_test" "${luciferase[@]}"

# . ./merge_replicate.sh "fusion_test" "${fusion[@]}"