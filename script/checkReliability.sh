# Check the reliability of the cut and run experiment by checking peaks in histone samples. If peaks in these samples reflect what has been known ( for example, if Integrated Genome Viewer on known gene regions is similar with literature) then the experiment is reliable and can be used for further analysis.  
# Nhung 11 04 2023
# get hg38_gene from here http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1605079819_JqZEx4ChvvwMuyzegwcGS5Q6AxPr&clade=mammal&org=0&db=0&hgta_group=genes&hgta_track=refSeqComposite&hgta_table=ncbiRefSeq&hgta_regionType=genome&position=&hgta_outputType=primaryTable&hgta_outFileName= and save as bed 

# plot heat map 
h3k4me3=( "SCC-ChIC-PMC-DRO-FH", "SCC-ChIC-PMC-DRO-LH", "SCC-ChIC-PMC-DRO-TH" )
h3k4me3_2=( "bulkChIC-PMC-DRO-011", "bulkChIC-PMC-DRO-012", "bulkChIC-PMC-DRO-013") 

computeMatrix reference-point --referencePoint TSS \
  -b 1000 -a 1000 \
  -R hg38_gene \
-S ${merged_bigwig_dir}/fusion_merged.bw\
 ${merged_bigwig_dir}/tfe3_merged.bw\
 ${merged_bigwig_dir}/luciferase_merged.bw\
 --skipZeros -o ${merged_bigwig_dir}/matrix_gene_hg38.mat.gz -p $cores

plotHeatmap -m ${merged_bigwig_dir}/matrix_gene_hg38.mat.gz -out ${merged_bigwig_dir}/transcript_hg38.png --sortUsing sum