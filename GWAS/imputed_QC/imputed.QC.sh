#***********************
# 1. prepare index file
#***********************
imputed_data="ABCD/ABCD_imputed/imputed"
ls $imputed_data/*.vcf.gz | while read file; do fname=$(basename $file); chr=$(echo $fname | awk '{split($1, a, ".dose"); print a[1]}'); echo -e "$chr\tdata/$fname"; done > index.tsv 

#***********************
# 2. run nextflow pipeline for filtering imputed variants
# dependencies needed are: tabix & bcftools
#***********************
nextflow run -bg filter.imputed.variants.nf --index index.tsv --outFolder results -with-trace -resume -c nextflow.config > pipeline.log

