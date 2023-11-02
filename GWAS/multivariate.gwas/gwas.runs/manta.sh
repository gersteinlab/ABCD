# install  mvgwas-nf (nf implementation of MANTA for GWAS)

git clone https://github.com/dgarrimar/mvgwas-nf
cd mvgwas-nf/


# run mvgwas 
nextflow run mvgwas.nf --l 1000 --geno all.chr.vcf.gz --pheno phenotypes.tsv --cov covariates.tsv --out mvgwas.tsv -resume -with-singularity -with-trace -bg -with-mpi
