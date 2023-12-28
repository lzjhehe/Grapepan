for i in $(cat fqlist);do
vg_64 giraffe -p -t 32 -Z grapeMC.d2.gbz -d grapeMC.d2.dist -m grapeMC.d2.min \
-f /fq/${i}_1.clean.fastq.gz -f /fq/${i}_2.clean.fastq.gz -o BAM -N ${i} > ${i}_MC.bam
samtools sort -O BAM -@ 6 -o /sortbam/${i}.sort.bam ${i}_MC.bam
samtools index /sortbam/${i}.sort.bam
gtx vc -r PN.MC.fa -g -i /sortbam/${i}_MC.sort.bam -o ${i}_MC.g.vcf.gz
done

gtx joint \
  -r PN.MC.fa \
  --sample-name-map sample_gvcf_map.tsv \
  -o allsnp.vcf.gz

bcftools view -i 'F_MISSING <= 0.2 & MAF > 0.05 & N_ALT = 1 & TYPE="snp" & F_PASS(FORMAT/DP >=5 & GT!="mis") >0.8' --threads 30 allsnp.vcf.gz | bcftools annotate -x INFO,FORMAT -O v> 2.snp.vcf
bgzip -@ 30 2.snp.vcf && tabix 2.snp.vcf.gz


bcftools view -i 'F_MISSING <= 0.2 & MAF > 0.05 & N_ALT = 1 & TYPE="indel" & F_PASS(FORMAT/DP >=5 & GT!="mis") >0.8' --threads 30 allsnp.vcf.gz | bcftools annotate -x INFO,FORMAT -O v> 2.indel.vcf
bgzip -@ 96  2.indel.vcf && tabix 2.indel.vcf.gz
