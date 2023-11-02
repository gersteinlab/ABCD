
# Digital phenotyping from wearables using AI characterizes psychiatric disorders and identifies genetic associations

## Requirements
```
tsai==0.3.1
sktime==0.14.1
```


## Installation
Create a conda environment with python 3.7:
```
conda create -n ABCD python=3.7
```

Activate the environment:
```
conda activate ABCD
```

Install ABCD and requirements:
```
pip install -e . --user
```

## Example Running
cd into the ABCD folder:
```
cd ABCD
```

Run model with preprocessing:
```
python run.py --preprocess True --add_genomics False \
--raw_data_path [PATH_OF_RAW_LABEL] --group [GROUP_NAME] --label_path [PATH_OF_LABEL] \
--cov_path [PATH_OF_COVARIATES] --out_path [OUTPUT_PATH] --save_path [SAVE_PATH]
```

Run model with preprocessing and genomics data:
```
python run.py --preprocess True --add_genomics True --raw_data_path [PATH_OF_RAW_LABEL] \
--group [GROUP_NAME] --label_path [PATH_OF_LABEL] \
--cov_path [PATH_OF_COVARIATES] --out_path [OUTPUT_PATH] --save_path [SAVE_PATH]
```

Load preprocessed files and run model:
```
python run.py --X_path [PATH_OF_X] --Y_path [PATH_OF_Y]
```

## Try Different Models
Supported models include 'RNN', 'TST', 'InceptionTimePlus', 'XceptionTimePlus', 'MiniRocket'.
Run different models:

Load preprocessed files and run model:
```
python run.py --X_path [PATH_OF_X] --Y_path [PATH_OF_Y] --model [MODEL_NAME]
```

## Model Interpretability
For step importance \& feature importance, make sure you are running tsai 0.3.1.


### Step Importance
```
python run.py --X_path [PATH_OF_X] --Y_path [PATH_OF_Y] --step_importance True
```

### Feature Importance
```
python run.py --X_path [PATH_OF_X] --Y_path [PATH_OF_Y] --feature_importance True
```

### GradCAM
```
python run.py --X_path [PATH_OF_X] --Y_path [PATH_OF_Y] --grad_cam True
```

## GWAS 

### Univariate GWAS
For univariate GWAS we employed [plink2](https://www.cog-genomics.org/plink/2.0/).

#### Binary trait 
```
plink2 \
  --threads 6 \
  --out out/gwas_sumstats \
  --keep indivs2keep.txt \
  --pfile imputed/pgen.files/genotype \
  --pheno phenotypes.binary.tsv \
  --chr 1-22,X \
  --covar covariates.tsv \
  --covar-variance-standardize \
  --glm firth-fallback hide-covar omit-ref no-x-sex
```
#### Continuous trait
```plink2 \
  --threads 6 \
  --out out/gwas_sumstats \
  --keep indivs2keep.txt \
  --pfile imputed/pgen.files/genotype \
  --pheno phenotypes.quantitative.tsv \
  --pheno-name cbcl_externalizing,cbcl_internalizing,liability_Xception_without_cbcl,liability_Xception_with_CBCL,liability_Xception_without_cbcl_v2,liability_Xception_with_CBCL_v2,XGB_without_cbcl,XGB_with_cbcl,XGB_without_cbcl_v2,XGB_with_cbcl_v2 \
  --chr 1-22,X \
  --covar covariates.tsv \
  --covar-variance-standardize \
  --glm hide-covar omit-ref no-x-sex
```

### Multivariate GWAS
For multivariate GWAS we employed [mvgwas-nf](https://github.com/dgarrimar/mvgwas-nf).
```
nextflow run mvgwas.nf --l 1000 --geno all.chr.vcf.gz --pheno phenotypes.tsv --cov covariates.tsv --out mvgwas.tsv -resume -with-singularity -with-trace -bg -with-mpi
```
