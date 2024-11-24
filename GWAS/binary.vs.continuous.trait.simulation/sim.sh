# 1. run simulations of binary and continuous traits
# for ten different values of MAF

for i in {0.01,0.025,0.05,0.07,0.09,0.1,0.25,0.3,0.4,0.5}; do 

        Rscript sim.R -m $i -o results;

done 

