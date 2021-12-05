#find the fastq.gz files
#head gz file to get the platform IDs: 

	gzip -cd Sample1_Lane2_Read1.fastq.gz | head


	cd /mnt/floyd_server/Daniel_Luo/Raw_data/gbm_049_merged/fastq

	sudo gzip -dk 1.fastq.gz
	sudo gzip -dk 2.fastq.gz
	sudo gzip -dk 3.fastq.gz
	sudo gzip -dk 4.fastq.gz
	sudo gzip -dk 5.fastq.gz
	sudo gzip -dk 6.fastq.gz
	sudo gzip -dk 7.fastq.gz
	sudo gzip -dk 8.fastq.gz
	sudo gzip -dk 9.fastq.gz
	sudo gzip -dk 10.fastq.gz



	cd /home/daniel/ubuntu/workspace/all_049/gbm_049_merged/alignments   #cd into allignments even though fastq files are not here. this is for samtools. 


#keep in mind hisat2 can only write to a local directory, not a networked directory. 
#keep in mind that samtools can be very resource intensive if you set the number of temp files it creates to 8 or more via -@ n.


	hisat2 -p 8 --rg-id=H53NGBBXY.1 --rg SM:1 --rg LB:1_NAGTGCTT+NCACCTCA --rg PL:ILLUMINA --rg PU:H53NGBBXY.1.NAGTGCTT+NCACCTCA -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness RF /mnt/floyd_server/Daniel_Luo/Raw_data/gbm_049_merged/fastq/1.fastq -S /home/daniel/ubuntu/workspace/all_049/gbm_049_merged/alignments/1.sam

	samtools sort -@ 4 -o 1.bam 1.sam

#make sure that sorting was successful
	rm 1.sam


	hisat2 -p 8 --rg-id=H53NGBBXY.2 --rg SM:2 --rg LB:2_CAGTGCTT+ACACCTCA --rg PL:ILLUMINA --rg PU:H53NGBBXY.2.CAGTGCTT+ACACCTCA -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness RF /mnt/floyd_server/Daniel_Luo/Raw_data/gbm_049_merged/fastq/2.fastq -S /home/daniel/ubuntu/workspace/all_049/gbm_049_merged/alignments/2.sam
	samtools sort -@ 4 -o 2.bam 2.sam
	rm 2.sam


	hisat2 -p 8 --rg-id=H53NGBBXY.1 --rg SM:3 --rg LB:3_NTGATCCA+NCTCAACG --rg PL:ILLUMINA --rg PU:H53NGBBXY.1.NTGATCCA+NCTCAACG -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness RF /mnt/floyd_server/Daniel_Luo/Raw_data/gbm_049_merged/fastq/3.fastq -S /home/daniel/ubuntu/workspace/all_049/gbm_049_merged/alignments/3.sam
	samtools sort -@ 4 -o 3.bam 3.sam
	rm 3.sam


	hisat2 -p 8 --rg-id=H53NGBBXY.2 --rg SM:4 --rg LB:4_GTGATCCA+ACTCAACG --rg PL:ILLUMINA --rg PU:H53NGBBXY.2.GTGATCCA+ACTCAACG -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness RF /mnt/floyd_server/Daniel_Luo/Raw_data/gbm_049_merged/fastq/4.fastq -S /home/daniel/ubuntu/workspace/all_049/gbm_049_merged/alignments/4.sam
	samtools sort -@ 4 -o 4.bam 4.sam

	rm 4.sam


	hisat2 -p 8 --rg-id=H53NGBBXY.1 --rg SM:5 --rg LB:5_NCAGCAAG+NTCTGCAA --rg PL:ILLUMINA --rg PU:H53NGBBXY.1.NCAGCAAG+NTCTGCAA -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness RF /mnt/floyd_server/Daniel_Luo/Raw_data/gbm_049_merged/fastq/5.fastq -S /home/daniel/ubuntu/workspace/all_049/gbm_049_merged/alignments/5.sam
	samtools sort -@ 4 -o 5.bam 5.sam
	rm 5.sam

	hisat2 -p 8 --rg-id=H53NGBBXY.2 --rg SM:6 --rg LB:6_ACAGCAAG+GTCTGCAA --rg PL:ILLUMINA --rg PU:H53NGBBXY.2.ACAGCAAG+GTCTGCAA -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness RF /mnt/floyd_server/Daniel_Luo/Raw_data/gbm_049_merged/fastq/6.fastq -S /home/daniel/ubuntu/workspace/all_049/gbm_049_merged/alignments/6.sam
	samtools sort -@ 4 -o 6.bam 6.sam
	rm 6.sam

	hisat2 -p 8 --rg-id=H53NGBBXY.1 --rg SM:7 --rg LB:7_NATACGGA+NACACCAC --rg PL:ILLUMINA --rg PU:H53NGBBXY.1.NATACGGA+NACACCAC -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness RF /mnt/floyd_server/Daniel_Luo/Raw_data/gbm_049_merged/fastq/7.fastq -S /home/daniel/ubuntu/workspace/all_049/gbm_049_merged/alignments/7.sam
	samtools sort -@ 4 -o 7.bam 7.sam
	rm 7.sam

	hisat2 -p 8 --rg-id=H53NGBBXY.2 --rg SM:8 --rg LB:8_CATACGGA+AACACCAC --rg PL:ILLUMINA --rg PU:H53NGBBXY.2.CATACGGA+AACACCAC -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness RF /mnt/floyd_server/Daniel_Luo/Raw_data/gbm_049_merged/fastq/8.fastq -S /home/daniel/ubuntu/workspace/all_049/gbm_049_merged/alignments/8.sam
	samtools sort -@ 4 -o 8.bam 8.sam
	rm 8.sam


	hisat2 -p 8 --rg-id=H53NGBBXY.1 --rg SM:9 --rg LB:9_NCGTGCAT+NTGTACCA --rg PL:ILLUMINA --rg PU:H53NGBBXY.1.NCGTGCAT+NTGTACCA -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness RF /mnt/floyd_server/Daniel_Luo/Raw_data/gbm_049_merged/fastq/9.fastq -S /home/daniel/ubuntu/workspace/all_049/gbm_049_merged/alignments/9.sam
	samtools sort -@ 4 -o 9.bam 9.sam
	rm 9.sam

	hisat2 -p 8 --rg-id=H53NGBBXY.2 --rg SM:10 --rg LB:10_TCGTGCAT+CTGTACCA --rg PL:ILLUMINA --rg PU:H53NGBBXY.2.TCGTGCAT+CTGTACCA -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness RF /mnt/floyd_server/Daniel_Luo/Raw_data/gbm_049_merged/fastq/10.fastq -S /home/daniel/ubuntu/workspace/all_049/gbm_049_merged/alignments/10.sam
	samtools sort -@ 4 -o 10.bam 10.sam
	rm 10.sam



#stay in folder with the bam files
#merge the files

	java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=049_t.bam INPUT=1.bam INPUT=2.bam
	java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=049_s.bam INPUT=3.bam INPUT=4.bam
	java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=049_n.bam INPUT=5.bam INPUT=6.bam
	java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=049_o.bam INPUT=7.bam INPUT=8.bam
	java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=049_v.bam INPUT=9.bam INPUT=10.bam


	samtools index 049_t.bam
	samtools index 049_s.bam
	samtools index 049_n.bam
	samtools index 049_o.bam
	samtools index 049_v.bam


	mkdir -p /home/daniel/ubuntu/workspace/all_049/gbm_049_merged/expression/stringtie/ref_only/
	cd /home/daniel/ubuntu/workspace/all_049/gbm_049_merged/expression/stringtie/ref_only/

	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 049_t/transcripts.gtf -A 049_t/gene_abundances.tsv /home/daniel/ubuntu/workspace/all_049/gbm_049_merged/alignments/049_t.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 049_s/transcripts.gtf -A 049_s/gene_abundances.tsv /home/daniel/ubuntu/workspace/all_049/gbm_049_merged/alignments/049_s.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 049_n/transcripts.gtf -A 049_n/gene_abundances.tsv /home/daniel/ubuntu/workspace/all_049/gbm_049_merged/alignments/049_n.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 049_o/transcripts.gtf -A 049_o/gene_abundances.tsv /home/daniel/ubuntu/workspace/all_049/gbm_049_merged/alignments/049_o.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 049_v/transcripts.gtf -A 049_v/gene_abundances.tsv /home/daniel/ubuntu/workspace/all_049/gbm_049_merged/alignments/049_v.bam



	mkdir -p /home/daniel/ubuntu/workspace/all_049/gbm_049_merged/expression/htseq_counts
	cd /home/daniel/ubuntu/workspace/all_049/gbm_049_merged/expression/htseq_counts


	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id /home/daniel/ubuntu/workspace/all_049/gbm_049_merged/alignments/049_t.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 049_t.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id /home/daniel/ubuntu/workspace/all_049/gbm_049_merged/alignments/049_s.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 049_s.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id /home/daniel/ubuntu/workspace/all_049/gbm_049_merged/alignments/049_n.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 049_n.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id /home/daniel/ubuntu/workspace/all_049/gbm_049_merged/alignments/049_o.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 049_o.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id /home/daniel/ubuntu/workspace/all_049/gbm_049_merged/alignments/049_v.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 049_v.tsv



#then continue with the R pipeline in next document

