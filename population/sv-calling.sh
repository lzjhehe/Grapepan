PanGenie-index -v grapeMC.vcf -r MC.PN.fa -t 10 -o /ref/PNref
for i in $(cat fqlist);do
zcat /fq/${i}_1.clean.fastq.gz > /fq/${i}.fq
zcat /fq/${i}_2.clean.fastq.gz >> /fq/${i}.fq
PanGenie -f /ref/PNref -i /fq/${i}.fq -s ${i} -t 19 -o /res/${i} && rm /fq/${i}.fq
bgzip -c -@ 5 /res/${i}_genotyping.vcf > /res/${i}_genotyping.vcf.gz && tabix -p /res/${i}_genotyping.vcf.gz