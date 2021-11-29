

cd $gbm_049/alignments


sudo gzip -dk 1.fastq.gz
sudo gzip -dk 2.fastq.gz
sudo gzip -dk 3.fastq.gz
sudo gzip -dk 4.fastq.gz
sudo gzip -dk 5.fastq.gz
sudo gzip -dk 6.fastq.gz
sudo gzip -dk 7.fastq.gz
sudo gzip -dk 8.fastq.gz
sudo gzip -dk 9_rep1.fastq.gz
sudo gzip -dk 9_rep2.fastq.gz
sudo gzip -dk 10_rep1.fastq.gz
sudo gzip -dk 10_rep2.fastq.gz


hisat2 -p 8 --rg-id=H53NGBBXY.1 --rg SM:1 --rg LB:1_NAGTGCTT+NCACCTCA --rg PL:ILLUMINA --rg PU:H53NGBBXY.1.NAGTGCTT+NCACCTCA -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness RF $gbm_049_server/fastq/1.fastq -S $gbm_049/alignments/1.sam
samtools sort -@ 8 -o 1.bam 1.sam
rm 1.sam

hisat2 -p 8 --rg-id=H53NGBBXY.2 --rg SM:2 --rg LB:2_CAGTGCTT+ACACCTCA --rg PL:ILLUMINA --rg PU:H53NGBBXY.2.CAGTGCTT+ACACCTCA -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness RF $gbm_049_server/fastq/2.fastq -S $gbm_049/alignments/2.sam
samtools sort -@ 8 -o 2.bam 2.sam
rm 2.sam

hisat2 -p 8 --rg-id=H53NGBBXY.1 --rg SM:3 --rg LB:3_NTGATCCA+NCTCAACG --rg PL:ILLUMINA --rg PU:H53NGBBXY.1.NTGATCCA+NCTCAACG -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness RF $gbm_049_server/fastq/3.fastq -S $gbm_049/alignments/3.sam
samtools sort -@ 8 -o 3.bam 3.sam
rm 3.sam

hisat2 -p 8 --rg-id=H53NGBBXY.2 --rg SM:4 --rg LB:4_GTGATCCA+ACTCAACG --rg PL:ILLUMINA --rg PU:H53NGBBXY.2.GTGATCCA+ACTCAACG -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness RF $gbm_049_server/fastq/4.fastq -S $gbm_049/alignments/4.sam
samtools sort -@ 8 -o 4.bam 4.sam
rm 4.sam

hisat2 -p 8 --rg-id=H53NGBBXY.1 --rg SM:5 --rg LB:5_NCAGCAAG+NTCTGCAA --rg PL:ILLUMINA --rg PU:H53NGBBXY.1.NCAGCAAG+NTCTGCAA -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness RF $gbm_049_server/fastq/5.fastq -S $gbm_049/alignments/5.sam
samtools sort -@ 8 -o 5.bam 5.sam
rm 5.sam

hisat2 -p 8 --rg-id=H53NGBBXY.2 --rg SM:6 --rg LB:6_ACAGCAAG+GTCTGCAA --rg PL:ILLUMINA --rg PU:H53NGBBXY.2.ACAGCAAG+GTCTGCAA -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness RF $gbm_049_server/fastq/6.fastq -S $gbm_049/alignments/6.sam
samtools sort -@ 8 -o 6.bam 6.sam
rm 6.sam

hisat2 -p 8 --rg-id=H53NGBBXY.1 --rg SM:7 --rg LB:7_NATACGGA+NACACCAC --rg PL:ILLUMINA --rg PU:H53NGBBXY.1.NATACGGA+NACACCAC -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness RF $gbm_049_server/fastq/7.fastq -S $gbm_049/alignments/7.sam
samtools sort -@ 8 -o 7.bam 7.sam
rm 7.sam

hisat2 -p 8 --rg-id=H53NGBBXY.2 --rg SM:8 --rg LB:8_CATACGGA+AACACCAC --rg PL:ILLUMINA --rg PU:H53NGBBXY.2.CATACGGA+AACACCAC -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness RF $gbm_049_server/fastq/8.fastq -S $gbm_049/alignments/8.sam
samtools sort -@ 8 -o 8.bam 8.sam
rm 8.sam

hisat2 -p 8 --rg-id=H53NGBBXY.1 --rg SM:9_rep1 --rg LB:9_NCGTGCAT+NTGTACCA --rg PL:ILLUMINA --rg PU:H53NGBBXY.1.NCGTGCAT+NTGTACCA -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness RF $gbm_049_server/fastq/9_rep1.fastq -S $gbm_049/alignments/9_rep1.sam
samtools sort -@ 8 -o 9_rep1.bam 9_rep1.sam
rm 9_rep1.sam

hisat2 -p 8 --rg-id=H53NGBBXY.2 --rg SM:9_rep2 --rg LB:9_TCGTGCAT+CTGTACCA --rg PL:ILLUMINA --rg PU:H53NGBBXY.2.TCGTGCAT+CTGTACCA -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness RF $gbm_049_server/fastq/9_rep2.fastq -S $gbm_049/alignments/9_rep2.sam
samtools sort -@ 8 -o 9_rep2.bam 9_rep2.sam
rm 9_rep2.sam


hisat2 -p 8 --rg-id=H53NGBBXY.1 --rg SM:10_rep1 --rg LB:10_NTACATCC+NCACCTAG --rg PL:ILLUMINA --rg PU:H53NGBBXY.1.NTACATCC+NCACCTAG -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness RF $gbm_049_server/fastq/10_rep1.fastq -S $gbm_049/alignments/10_rep1.sam
samtools sort -@ 8 -o 10_rep1.bam 10_rep1.sam
rm 10_rep1.sam


hisat2 -p 8 --rg-id=H53NGBBXY.2 --rg SM:10_rep2 --rg LB:10_CTATATCC+TCACCTAG --rg PL:ILLUMINA --rg PU:H53NGBBXY.2.CTATATCC+TCACCTAG -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness RF $gbm_049_server/fastq/10_rep2.fastq -S $gbm_049/alignments/10_rep2.sam
samtools sort -@ 8 -o 10_rep2.bam 10_rep2.sam
rm 10_rep2.sam




















#navigate to folder where the bam files are
cd $gbm_049/alignments

java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=9.bam INPUT=9_rep1.bam INPUT=9_rep2.bam
java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=10.bam INPUT=10_rep1.bam INPUT=10_rep2.bam



samtools index 1.bam
samtools index 2.bam
samtools index 3.bam
samtools index 4.bam
samtools index 5.bam
samtools index 6.bam
samtools index 7.bam
samtools index 8.bam
samtools index 9.bam
samtools index 10.bam



mkdir -p $gbm_049/expression/stringtie/ref_only/
cd $gbm_049/expression/stringtie/ref_only/

stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 1/transcripts.gtf -A 1/gene_abundances.tsv $gbm/alignments/1.bam
stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 2/transcripts.gtf -A 2/gene_abundances.tsv $gbm/alignments/2.bam

stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 3/transcripts.gtf -A 3/gene_abundances.tsv $gbm/alignments/3.bam
stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 4/transcripts.gtf -A 4/gene_abundances.tsv $gbm/alignments/4.bam

stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 5/transcripts.gtf -A 5/gene_abundances.tsv $gbm/alignments/5.bam
stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 6/transcripts.gtf -A 6/gene_abundances.tsv $gbm/alignments/6.bam

stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 7/transcripts.gtf -A 7/gene_abundances.tsv $gbm/alignments/7.bam
stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 8/transcripts.gtf -A 8/gene_abundances.tsv $gbm/alignments/8.bam

stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 9/transcripts.gtf -A 9/gene_abundances.tsv $gbm/alignments/9.bam
stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 10/transcripts.gtf -A 10/gene_abundances.tsv $gbm/alignments/10.bam





mkdir -p $gbm_049/expression/htseq_counts
cd $gbm_049/expression/htseq_counts


htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_049/alignments/1.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 1.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_049/alignments/2.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 2.tsv

htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_049/alignments/3.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 3.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_049/alignments/4.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 4.tsv

htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_049/alignments/5.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 5.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_049/alignments/6.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 6.tsv

htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_049/alignments/7.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 7.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_049/alignments/8.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 8.tsv

htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_049/alignments/9.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 9.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_049/alignments/10.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 10.tsv




mkdir -p $gbm_049/de/ballgown/ref_only
cd $gbm_049/de/ballgown/ref_only/


ids type path-to-file-011_invitro_1 011 $gbm/expression/stringtie/1 011_invitro_2 011 $gbm/expression/stringtie/2 ... ...


printf "\"ids\",\"type\",\"path
\"\n\"1\",\"049_tissue\",\"$gbm_049/expression/stringtie/ref_only/1
\"\n\"2\",\"049_tissue\",\"$gbm_049/expression/stringtie/ref_only/2
\"\n\"3\",\"049_slice\",\"$gbm_049/expression/stringtie/ref_only/3
\"\n\"4\",\"049_slice\",\"$gbm_049/expression/stringtie/ref_only/4
\"\n\"5\",\"049_discells\",\"$gbm_049/expression/stringtie/ref_only/5
\"\n\"6\",\"049_discells\",\"$gbm_049/expression/stringtie/ref_only/6
\"\n\"7\",\"049_organoid\",\"$gbm_049/expression/stringtie/ref_only/7
\"\n\"8\",\"049_organoid\",\"$gbm_049/expression/stringtie/ref_only/8
\"\n\"9\",\"049_invitro\",\"$gbm_049/expression/stringtie/ref_only/9
\"\n\"10\",\"049_invitro\",\"$gbm_049/expression/stringtie/ref_only/10
\"\n" > GBM_049_all.csv


printf "\"ids\",\"type\",\"path\"\n\"1\",\"049_tissue\",\"$gbm_049/expression/stringtie/ref_only/1\"\n\"2\",\"049_tissue\",\"$gbm_049/expression/stringtie/ref_only/2\"\n\"3\",\"049_slice\",\"$gbm_049/expression/stringtie/ref_only/3\"\n\"4\",\"049_slice\",\"$gbm_049/expression/stringtie/ref_only/4\"\n\"5\",\"049_discells\",\"$gbm_049/expression/stringtie/ref_only/5\"\n\"6\",\"049_discells\",\"$gbm_049/expression/stringtie/ref_only/6\"\n\"7\",\"049_organoid\",\"$gbm_049/expression/stringtie/ref_only/7\"\n\"8\",\"049_organoid\",\"$gbm_049/expression/stringtie/ref_only/8\"\n\"9\",\"049_invitro\",\"$gbm_049/expression/stringtie/ref_only/9\"\n\"10\",\"049_invitro\",\"$gbm_049/expression/stringtie/ref_only/10\"\n" > GBM_049_all.csv



R --no-restore
	library(ballgown)
	library(genefilter)
	library(dplyr)
	library(devtools)
	library(ggplot2)
	library(gplots)
	library(GenomicRanges)

	pheno_data = read.csv("GBM_049_all.csv")  


	bg = ballgown(samples=as.vector(pheno_data$path), pData=pheno_data)
	bg


		bg_table = texpr(bg, 'all')


		bg_gene_names = unique(bg_table[, 9:10])
		head(bg_gene_names)


		save(bg, file='bg.rda')

		bg



		pdf(file="GBM_049_R_output.pdf")

		working_dir = "~/workspace/gbm_049/de/ballgown/ref_only"
		setwd(working_dir)
		dir()



		
		load('bg.rda')



		
		bg_table = texpr(bg, 'all')
		bg_gene_names = unique(bg_table[, 9:10])




		gene_expression = as.data.frame(gexpr(bg))
		head(gene_expression)



		colnames(gene_expression)


		
		row.names(gene_expression)



		dim(gene_expression)




		i = row.names(gene_expression) == "BRD4"
		gene_expression[i,]

# output:
#       FPKM.1   FPKM.2   FPKM.3   FPKM.4   FPKM.5   FPKM.6   FPKM.7   FPKM.8
#BRD4 21.61065 27.15311 7.702215 8.867953 4.289441 9.266627 8.182026 9.918504
#       FPKM.9  FPKM.10
#BRD4 15.50656 15.95751









		transcript_gene_table = indexes(bg)$t2g
		head(transcript_gene_table)
		
		length(row.names(transcript_gene_table)) #Transcript count
		length(unique(transcript_gene_table[,"g_id"])) #Unique Gene count








	#Plot the number of transcripts per gene. Many genes will have only 1 transcript, some genes will have several transcripts. Use the 'table()' command #to count the number of times each gene symbol occurs (i.e. the # of transcripts that have each gene symbol). Then use the 'hist' command to create a #histogram of these counts

		counts=table(transcript_gene_table[,"g_id"])
		c_one = length(which(counts == 1))
		c_more_than_one = length(which(counts > 1))
		c_max = max(counts)
		hist(counts, breaks=50, col="bisque4", xlab="Transcripts per gene", main="Distribution of transcript count per gene")
		legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
		legend("topright", legend_text, lty=NULL)





	#Plot the distribution of transcript sizes as a histogram. lengths will be those of known transcripts. Good QC step: we had a low coverage library, or #other problems, we might get short 'transcripts' that are actually only pieces of real transcripts.

		full_table <- texpr(bg , 'all')
		hist(full_table$length, breaks=500, xlab="Transcript length (bp)", main="Distribution of transcript lengths", col="steelblue")


		
		min(gene_expression[,"FPKM.1"])
		max(gene_expression[,"FPKM.2"])


		
		min_nonzero=1



		data_columns=c(1:10)
		short_names=c("tissue_1","tissue_2","slice_1","slice2","discells_1", "discells_2","organoid_1","organoid_2","invitro_1","invitro_2")





		colors()
		data_colors=c("tomato1","tomato2","royalblue1","royalblue2","grey1","grey2","seagreen1","seagreen2","grey3","grey4")

		boxplot(log2(gene_expression[,data_columns]+min_nonzero), col=data_colors, names=short_names, las=2, ylab="log2(FPKM)", main="Distribution of FPKMs for all 10 sample libraries")









---------------------------
optional

#plot a pair of replicates to assess reproducibility of technical replicates. 
#Tranform the data by converting to log2 scale after adding an arbitrary small value to avoid log2(0). Also add a straight line of slope 1, and #intercept 0. Also calculate the correlation coefficient and display in a legend.

	x = gene_expression[,"FPKM.1"]
	y = gene_expression[,"FPKM.2"]
	plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="FPKM (049_slice, Replicate 1)", ylab="FPKM (049_slice, Replicate 2)", main="Comparison of expression values for replicates of slice samples")
	abline(a=0,b=1)
	rs=cor(x,y)^2
	legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")


			check organoid samples

				x = gene_expression[,"FPKM.3"]
				y = gene_expression[,"FPKM.4"]
				plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="FPKM (011_organoid, Replicate 1)", ylab="FPKM (011_organoid, Replicate 2)", main="Comparison of expression values for replicates of organoid samples")
				abline(a=0,b=1)
				rs=cor(x,y)^2
				legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")

			check tissue samples

				x = gene_expression[,"FPKM.5"]
				y = gene_expression[,"FPKM.6"]
				plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="FPKM (011_tissue, Replicate 1)", ylab="FPKM (011_tissue, Replicate 2)", main="Comparison of expression values for replicates of tissue samples")
				abline(a=0,b=1)
				rs=cor(x,y)^2
				legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")


			check in vitro samples

				x = gene_expression[,"FPKM.7"]
				y = gene_expression[,"FPKM.8"]
				plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="FPKM (011_invitro, Replicate 1)", ylab="FPKM (011_invitro, Replicate 2)", main="Comparison of expression values for replicates of in vitro samples")
				abline(a=0,b=1)
				rs=cor(x,y)^2
				legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")

---------------------------------









## Compare the correlation distance between all replicates

#Calculate the FPKM sum for all 10 libraries

	gene_expression[,"sum"]=apply(gene_expression[,data_columns], 1, sum)



#Filter out genes with a grand sum FPKM of less than 5 to filter out low expressing genes

	i = which(gene_expression[,"sum"] > 5)

# visually check raw FPKM counts 
	gene_expression[order(gene_expression$sum, decreasing = TRUE),]  %>% head



#Calculate the correlation between all pairs of data

	r=cor(gene_expression[i,data_columns], use="pairwise.complete.obs", method="pearson")
	r
	
	
## Plot MDS.
#Convert correlation to distance, and use 'multi-dimensional scaling' to plot the relative differences between libraries, by calculating 2-dimensional #coordinates to plot points for each library using eigenvectors (eig=TRUE). d, k=2 means 2 dimensions
	
	d=1-r
	mds=cmdscale(d, k=2, eig=TRUE)
	par(mfrow=c(1,1))
	plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all non-zero genes)", xlim=c(-0.4,0.4), ylim=c(-0.4,0.4))
	points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
	text(mds$points[,1], mds$points[,2], short_names, col=data_colors)




#Calculate the differential expression results including significance

	results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
	results_genes = merge(results_genes,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))



#close out the PDF
dev.off()


