#!/bin/bash
#Performs all tests for MatrixeQTLGLM
#load snps and genes into Rdata objects
queue="mini"
time="5"
memory=1000

bsub -q $queue -W $time -K -R "rusage[mem=$memory]" -R "select[model=XeonL5640]" 'Rscript  ~/glm_eqtl/MatrixeQTLGLM/load_static.R SNPEXP F testSNP.txt testexp.txt test.Rdata'&
bsub -q $queue -W $time -K -R "rusage[mem=$memory]" -R "select[model=XeonL5640]" 'Rscript ~/glm_eqtl/MatrixeQTLGLM/load_static.R ANNO testSNPanno.txt testgeneanno.txt testanno.Rdata'&
wait
bsub -q $queue -W $time -K -R "rusage[mem=$memory]" -R "select[model=XeonL5640]" 'Rscript ~/glm_eqtl/MatrixeQTLGLM/maineQTL.R test_eqtls . eqtlfolds/ testanno.Rdata test.Rdata 100 10 1:00 CISTRA'&
wait
bsub -q $queue -W $time -R -K "rusage[mem=$memory]" -R "select[model=XeonL5640]" 'Rscript ../load_eqtls.R eqtlfolds testdb.db 10 100 testSNP.txt T' &
wait
head -1 testSNP.txt > neqtlsnps.txt
bsub -q $queue -W $time -R "rusage[mem=$memory]" -R "select[model=XeonL5640]" 'grep -Fw -f eqtlsnps.txt testSNP.txt >> neqtlsnps.txt'&
head -1 testexp.txt > neqtlgenes.txt
bsub -q $queue -W $time -R "rusage[mem=$memory]" -R "select[model=XeonL5640]" 'grep -Fw -f eqtlgenes.txt testexp.txt >> neqtlgenes.txt' &
wait

bsub -q $queue -W $time -K -R "rusage[mem=$memory]" -R "select[model=XeonL5640]" 'Rscript  ~/glm_eqtl/MatrixeQTLGLM/snp_exp_melt.R neqtlsnps.txt SNPS 2 761 test_snp_melt.txt'&
echo "SNPs melted!"
bsub -q $queue -W $time -K -R "rusage[mem=$memory]" -R "select[model=XeonL5640]" 'Rscript ~/glm_eqtl/MatrixeQTLGLM/snp_exp_melt.R neqtlgenes.txt GENE test_gene_melt.txt' &
wait
bsub -q $queue -W $time -K -R "rusage[mem=$memory]" -R "select[model=XeonL5640]" '~/glm_eqtl/MatrixeQTLGLM/import-index.sh -i -s -f test_snp_melt.txt -d testdb.db' &
echo "imported SNPS started!"
bsub -q $queue -W $time -K -R "rusage[mem=$memory]" -R "select[model=XeonL5640]" '~/glm_eqtl/MatrixeQTLGLM/import-index.sh -i -g -f test_gene_melt.txt -d testdb.db' &
echo "imported genes started!"
wait
echo "import finished!"
bsub -q $queue -W $time -K -R "rusage[mem=$memory]" -R "select[model=XeonL5640]" '~/glm_eqtl/MatrixeQTLGLM/import-index.sh -x -s -d testdb.db'&
bsub -q $queue -W $time -K -R "rusage[mem=$memory]" -R "select[model=XeonL5640]" '~/glm_eqtl/MatrixeQTLGLM/import-index.sh -x -g -d testdb.db'&
wait
echo "index finished too!"

