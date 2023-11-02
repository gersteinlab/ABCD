
# this script is to simulate the power of a 
# continuous vs. binary trait

.libPaths("/gpfs/gibbs/pi/gerstein/bb926/R/libraries")


#************
# LIBRARIES *
#************

library(ggplot2)
library(MASS)

#************
# OPTION PARSING *
#************


suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  
  make_option( c ( "-m", "--maf" ), default=NULL,
               help = "The specified maf."),
  
  make_option( c( "-o", "--output_folder" ), default = NULL,
               help = "Output folder where to save the residual results from each simulation[default = %default]." )
  
)

parser <- OptionParser(
  usage = "%prog [options]", 
  option_list=option_list,
  description = "\n Power simulation for linear and logistic models for gwas."
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options


#********
# BEGIN *
#********

#--------------------------------------------------------------------------
# 1. simulate a cohort of 1,500 individuals                               |
# with a genotype at a particular SNP following a specified MAF           |
# we follow: https://github.com/dgarrimar/manta-sim/blob/sim/bin/binGT.R  |
#--------------------------------------------------------------------------

## 1.1. specify MAF
set.seed(123)
maf = as.numeric(opt$maf)


## 1.2. specify number of individuals
n = 1500 

## 1.3. we model this with a binomial distribution
## setting the number of trials equal to 2 (because we have two alleles for each individual)
## the probability of success (i.e. probability of having an alternative allele)
## is given by MAF
## x is the vector of genotypes of the 1,500 individuals
## 0 = 0|0 (e.g. AA)
## 1 = 1|0 or 0|1 (e.g. AT or TA)
## 2 = 1|1 (e.g. TT)
x = rbinom(n, 2, maf)

table(x)


# 2. instead of simulating different degrees of correlation between genotype and trait
# we simulate 10 different betas ranging between 0 and 1
set.seed(123)
betas <- sort(runif(50, 0, 1))


# 3. dataframe to store results & define alpha (false positive rate)
res <- c()
alpha <- 0.05


# 4. run simulations for each beta
for ( b in betas ) {
  
  # 4.1. vector to store p-values from lm
  p.lm <- c()
  
  # 4.2. vector to store p-values from logistic regression
  p.log <- c()
  
  # 4.3. run 10,000 simulations, each time with a different set of residuals
  set.seed(123)
  for (i in 1:10000){
    
    # 4.3.1. simulate 1,500 normally distributed residuals (one residual for each individual)
    e <- rnorm(n)
    
    # 4.3.2. construct y for lm and run lm
    y.lm <- x*b + e
    p.lm <- c(p.lm, summary(lm(y.lm ~ x))$coefficients[2,4]) 
    
    # 4.3.3. binarize y and run logistic regression
    y.log <- ifelse(y.lm > median(y.lm), 1, 0)
    p.log <- c(p.log, summary(glm(y.log ~ x, family = "binomial"))$coefficients[8])
    
  }
  
  # 4.4. perform multiple-testing correction on p-values
  p.lm.fdr <- p.adjust(p.lm, method = "BH")
  p.log.fdr <- p.adjust(p.log, method = "BH")
  
  res <- rbind(res, data.frame(b, 
                               power_lm = length(p.lm.fdr[p.lm.fdr < alpha]),
                               power_log = length(p.log.fdr[p.log.fdr < alpha])))
  
}

#View(res)

#a = as.numeric(opt$mafs)
#save.image(res, file = sprintf("/simulations.%s.RData", a))

#save.RDS(res, file = sprintf("/residual.%s.tsv.", a))

write.table(res, file = paste0(opt$output_folder, "/res_", maf, ".tsv"), 
            row.names =F, col.names = T, sep="\t", quote = F)

