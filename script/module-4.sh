# motif analysis 
echo "-------------------step 14. running motif finding----------------------"

# define where to save motif result
motif_sub_dir=${motif_dir}/luc_top_rpkm20_antitfe3Normalization
mkdir -p $motif_sub_dir

# define bed file of peaks regions to find motif (in this case, the bed file is created from peak annotation script)
bed_dir=${res_dir}/diffBind_analysis/antiTFE3_Samples/luc_top_sites_RPKM_thredshold_20.bed 
# find motif. change size and len if needed
findMotifsGenome.pl $bed_dir $fasta_genome_dir ${motif_sub_dir} -size 100 -len 8,10,12

# find peak location from a motif (e.g. motif1 in homerResults)
findMotifsGenome.pl $bed_dir $fasta_genome_dir ${motif_sub_dir} -size 200 -len 8,10,12 -find $motif_sub_dir/homerResults/motif1.motif > $motif_sub_dir/motif1_location

