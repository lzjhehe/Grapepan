#PCA 
plink --vcf 2.snp.vcf --make-bed --out allsnp
plink --bfile allsnp --pca 20 --out pca

#IBD
plink --bfile allsnp --genome --out IBD


#pi
for i in Wild VL Vvwild Table1 Wine Table2
do 
echo "vcftools --vcf \${path}/2.snp.vcf  --keep \${path}/${i}.txt --window-pi 10000 --window-pi-step 5000 --out ${i}.pi.txt" | cat header - > ${i}.pi.sh
sh ${i}.pi.sh
done

#Fst
for i in Wild VL Vvwild Table1 Wine Table2
do for a in Wild VL Vvwild Table1 Wine Table2
do 
if [ $i != $a ];then 
echo "vcftools --vcf \${path}/2.snp.vcf --weir-fst-pop \${path}/${i}.txt --weir-fst-pop \${path}/${a}.txt --fst-window-size 10000 --fst-window-step 5000 --out ${i}_vs_${a}_fst.txt" | cat header - > ${i}_vs_${a}_fst.sh
sh ${i}_vs_${a}_fst.sh
fi
done
done

#Tajima'D
for i in Wild VL Vvwild Table1 Wine Table2
do echo "vcftools --vcf \${path}/2.snp.vcf --keep \${path}/${i}.txt --TajimaD 5000 --out ${i}_TajimaD" | cat header - > ${i}_TajimaD.sh
sh ${i}_TajimaD.sh
done

#tree

#admixture


#deleterrious burden


#GWAS
bcftools view -S gwassample -O v 2.snp.vcf.gz > 3.gwas.snp.vcf
bcftools view -S gwassample -O v 2.sv.vcf.gz > 3.gwas.sv.vcf
bcftools concat -a --threads 20 -O v 3.gwas.snp.vcf.gz 3.gwas.sv.vcf.gz > 3.gwas.snp.sv.vcf
plink --vcf 3.gwas.snp.sv.vcf --keep-allele-order --make-bed --const-fid --maf 0.05 --hwe 1e-5 --out allsnpsv.maf.hwe
plink --bfile allsnpsv.maf.hwe --pca 20 --out pca.maf.hwe
python3 /permGWAS/create_h5_file.py -x allsnpsv.maf.hwe.bed
for i in {1..29};do
for j in 2016 2017;do
python /permGWAS/permGWAS.py -x allsnpsv.maf.hwe..h5 -y /pheno${j}/pheno${i}.txt --cov pca.maf.hwe.csv --plot --qqplot --perm 1000 --out_file GWAS${i} --out_dir ${j}
done
done

#her
ldak5.2.linux --bfile allsnpsv.maf.hwe --window-prune 0.98 --window-kb 100 --thin first
awk < first.in '{print $1, 1}' > weights.thin
awk '{print $2}' allsnpsv.maf.hwe.bim| grep 'snp' > snplist
awk '{print $2}' allsnpsv.maf.hwe.bim| grep 'sv' > svlist
plink --bfile allsnpsv.maf.hwe --extract snplist --make-bed --keep-allele-order --out snp
plink --bfile allsnpsv.maf.hwe --extract svlist --make-bed --keep-allele-order --out sv

for j in {-20..10}; do
alpha=`echo $j | awk '{print $1/20}'`; echo $alpha
ldak5.2.linux --calc-kins-direct LDAK-Thin-snp-$alpha --bfile ../snp --weights ../weights.thin --power $alpha
ldak5.2.linux --calc-kins-direct LDAK-Thin-sv-$alpha --bfile ../sv --weights ../weights.thin --power $alpha
done

for i in {1..29};
do seq -20 10| awk '{print $1/20}'| parallel -j 10 "ldak5.2.linux --reml ./reml/alpha{}_T${i}_snp --pheno trait_1 --mpheno ${i} --grm LDAK-Thin-snp-{} --constrain YES --covar pca5"
for j in {-20..10}; do alpha=`echo $j | awk '{print $1/20}'`; grep Alt_Likelihood ./reml/alpha${alpha}_T${i}_snp.reml |awk -v alpha=${alpha} '{print alpha, $2}' >> ${i}_snp_Alt_Likelihood.txt ; done
ldak5.2.linux --find-gaussian ${i}_snp_alpha --likelihoods ${i}_snp_Alt_Likelihood.txt
done
for i in {1..29};
do seq -20 10| awk '{print $1/20}'| parallel -j 10 "ldak5.2.linux --reml ./reml/alpha{}_T${i}_sv --pheno trait_1 --mpheno ${i} --grm LDAK-Thin-sv-{} --constrain YES --covar pca5"
for j in {-20..10}; do alpha=`echo $j | awk '{print $1/20}'`; grep Alt_Likelihood ./reml/alpha${alpha}_T${i}_sv.reml |awk -v alpha=${alpha} '{print alpha, $2}' >> ${i}_sv_Alt_Likelihood.txt ; done
ldak5.2.linux --find-gaussian ${i}_sv_alpha --likelihoods ${i}_sv_Alt_Likelihood.txt
done

grep Power *.power > powerall.txt

echo LDAK-Thin-snp >list
echo LDAK-Thin-sv >>list
seq 1 29 | parallel -j 10 'ldak5.2.linux --reml reml/trait1_{}_snpsv --pheno trait_1 --mpheno {} --mgrm list --constrain YES --covar pca5'
for i in {1..29};do echo $i >> trait1_allresult;tail -n 5 trait1_${i}_snpsv.reml >> trait1_allresult;done

for i in 1 2;do
for j in {1..29};do
awk -F ',' -v OFS='\t' '{print \$7,\$3}' p_values_P${i}GWAS${j}.csv|sed '1d' |paste left - | cat header - > P${i}_${j}.summary
Rscript LDAKsum.R P${i}_${j}.summary 332 P${i}_${j} 
done
done

for i in {1..29};do 
    for j in {1..29};do 
    if [[ "${i}" -lt "${j}" ]];then
        echo "ldak5.2.linux --sum-cors ./cor/gencor_P1_${i}_${j} --summary P1_${i}.sumstats --summary2 P1_${j}.sumstats --tagfile LDAK-Thin.tagging" >>c.sh
        #echo $i
    fi
    done
done
sh c.sh

echo "Component Value SD">>allresult
for i in {1..29};do 
    for j in {1..29};do 
    if [[ "${i}" -lt "${j}" ]];then
        echo -e "$i ${j}asd" >> Coher_All.txt
        grep "Coher_All" ./cor/gencor_P1_${i}_${j}.cors >>Coher_All.txt
        echo -e "$i ${j}asd" >> Cor_All.txt
        grep "Cor_All" ./cor/gencor_P1_${i}_${j}.cors >>Cor_All.txt
    fi
    done
done
