# 049_re-analysis

## below scripts are for data in /home/daniel/ubuntu/workspace/all_049/gbm_049_unmerged

	nano .bashrc
	export gbm_049=/home/daniel/ubuntu/workspace/all_049/gbm_049_unmerged


#navigate to folder where the bam files are

	java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=9.bam INPUT=9_rep1.bam INPUT=9_rep2.bam
	java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=10.bam INPUT=10_rep1.bam INPUT=10_rep2.bam


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


format: ids type path-to-file-011_invitro_1 011 $gbm/expression/stringtie/path


printf "\"ids\",\"type\",\"path
\"\n\"1\",\"049_tissue\",\"$gbm_049/expression/stringtie/ref_only/1
\"\n\"2\",\"049_tissue\",\"$gbm_049/expression/stringtie/ref_only/2
\"\n\"3\",\"049_slice\",\"$gbm_049/expression/stringtie/ref_only/3
\"\n\"4\",\"049_slice\",\"$gbm_049/expression/stringtie/ref_only/4
\"\n\"5\",\"049_neurosphere\",\"$gbm_049/expression/stringtie/ref_only/5
\"\n\"6\",\"049_neurosphere\",\"$gbm_049/expression/stringtie/ref_only/6
\"\n\"7\",\"049_organoid\",\"$gbm_049/expression/stringtie/ref_only/7
\"\n\"8\",\"049_organoid\",\"$gbm_049/expression/stringtie/ref_only/8
\"\n\"9\",\"049_invitro\",\"$gbm_049/expression/stringtie/ref_only/9
\"\n\"10\",\"049_invitro\",\"$gbm_049/expression/stringtie/ref_only/10
\"\n" > GBM049_all.csv


script:

	printf "\"ids\",\"type\",\"path\"\n\"1\",\"049_tissue\",\"$gbm_049/expression/stringtie/ref_only/1\"\n\"2\",\"049_tissue\",\"$gbm_049/expression/stringtie/ref_only/2\"\n\"3\",\"049_slice\",\"$gbm_049/expression/stringtie/ref_only/3\"\n\"4\",\"049_slice\",\"$gbm_049/expression/stringtie/ref_only/4\"\n\"5\",\"049_neurosphere\",\"$gbm_049/expression/stringtie/ref_only/5\"\n\"6\",\"049_neurosphere\",\"$gbm_049/expression/stringtie/ref_only/6\"\n\"7\",\"049_organoid\",\"$gbm_049/expression/stringtie/ref_only/7\"\n\"8\",\"049_organoid\",\"$gbm_049/expression/stringtie/ref_only/8\"\n\"9\",\"049_invitro\",\"$gbm_049/expression/stringtie/ref_only/9\"\n\"10\",\"049_invitro\",\"$gbm_049/expression/stringtie/ref_only/10\"\n" > GBM049_all.csv




	R --no-restore
	library(ballgown)
	library(genefilter)
	library(dplyr)
	library(devtools)
	library(ggplot2)
	library(gplots)
	library(GenomicRanges)

	pheno_data = read.csv("GBM049_all.csv")  


	bg = ballgown(samples=as.vector(pheno_data$path), pData=pheno_data)
	bg


	bg_table = texpr(bg, 'all')


	bg_gene_names = unique(bg_table[, 9:10])
	head(bg_gene_names)


	save(bg, file='bg.rda')

	bg


	

	working_dir = "/home/daniel/ubuntu/workspace/all_049/gbm_049_unmerged/de/ballgown/ref_only"
	setwd(working_dir)
	dir()



		
	load('bg.rda')
		
	bg_table = texpr(bg, 'all')
	bg_gene_names = unique(bg_table[, 9:10])




	gene_expression = as.data.frame(gexpr(bg))
	head(gene_expression)



	colnames(gene_expression)


	dim(gene_expression)
	
#The data frame is always the format table[rows 1 through whatever, c(columns 1 through whatever,any other column you need added)]. To view the first 3 rows of gene #expression data for the first 4 samples plus the sample in the 6th column (1:3 means 1 through 3).
	
	gene_expression[1:3, c(1:4,6)]

#To view gene expression for a single gene by pulling out row names that matches (or ==) BRD4, across all columns (that weird [i,] which means i want every column, saves #typing c(1:10) by just not typing anything)

	i = row.names(gene_expression) == "BRD4"
	gene_expression[i,]



	transcript_gene_table = indexes(bg)$t2g
	head(transcript_gene_table)
		
	length(row.names(transcript_gene_table)) #Transcript count
	length(unique(transcript_gene_table[,"g_id"])) #Unique Gene count




## Plot the number of transcripts per gene. 


	pdf(file="GBM049_R_output.pdf") #open up a pdf writer

#Many genes will have only 1 transcript, some genes will have several transcripts. Use the 'table()' command #to count the number of times each gene symbol occurs (i.e. #the # of transcripts that have each gene symbol). Then use the 'hist' command to create a #histogram of these counts

	
	counts=table(transcript_gene_table[,"g_id"])
	c_one = length(which(counts == 1))
	c_more_than_one = length(which(counts > 1))
	c_max = max(counts)
	hist(counts, breaks=50, col="bisque4", xlab="Transcripts per gene", main="Distribution of transcript count per gene")
	legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
	legend("topright", legend_text, lty=NULL)





## Plot the distribution of transcript sizes as a histogram. lengths will be those of known transcripts. 
#Good QC step: we had a low coverage library, or other problems, we might get short 'transcripts' that are actually only pieces of real transcripts.

	full_table <- texpr(bg , 'all')
	hist(full_table$length, breaks=500, xlab="Transcript length (bp)", main="Distribution of transcript lengths", col="steelblue")


		
	min(gene_expression[,"FPKM.1"])
	max(gene_expression[,"FPKM.2"])


		
	min_nonzero=1



	data_columns=c(1:10)
	short_names=c("tissue_1","tissue_2","slice_1","slice2","discells_1", "discells_2","organoid_1","organoid_2","invitro_1","invitro_2")




	#colors()
	data_colors=c("tomato1","tomato2","royalblue1","royalblue2","seagreen1","seagreen2","grey1","grey2","brown1","brown2")

	boxplot(log2(gene_expression[,data_columns]+min_nonzero), col=data_colors, names=short_names, las=2, ylab="log2(FPKM)", main="Distribution of FPKMs for libraries of all 10 samples")




---------------------------
optional

#plot a pair of replicates to assess how far apart are my technical replicates and the reproducibility of these technical replicates. 
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



	dev.off() #close off the pdf and start a new one for MDS
	

	
	pdf(file="GBM049_R_output_MDS.pdf")


## Compare the correlation distance between all replicates


#Scree plot: Determine the amount of variance coming from each principal component in a table:

	pc <- princomp(gene_expression[,data_columns],cor=TRUE,scores=TRUE)
	summary(pc)
	plot(pc,type='lines')
	

#Calculate the FPKM sum for all 10 libraries

	gene_expression[,"sum"]=apply(gene_expression[,data_columns], 1, sum)

#Filter out genes with a grand sum FPKM of less than 10

	i = which(gene_expression[,"sum"] > 10)


#Calculate the correlation between all pairs of data

	r=cor(gene_expression[i,data_columns], use="pairwise.complete.obs", method="pearson")
	r
	


## Plot MDS.
#Convert correlation to distance, and use 'multi-dimensional scaling' to plot the relative differences between libraries, by calculating 2-dimensional #coordinates to #plot points for each library using eigenvectors (eig=TRUE). d, k=2 means 2 dimensions
	
	d=1-r
	mds=cmdscale(d, k=2, eig=TRUE)
	par(mfrow=c(1,1))
	plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all non-zero genes)", xlim=c(-0.4,0.4), ylim=c(-0.2,0.2))
	points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
	text(mds$points[,1], mds$points[,2], short_names, col=data_colors)



	dev.off()






## comparing invitro to non-invitro using stattest function


Need to make a new header file and recode the "type" header,since this will determine what gets compared in a stattest. the stattest function will use type as a covariate and use fpkm as a meansurement. since this function can't compare multiple things, need to make another file called GBM049_all_stattest.csv and make the "type" tissue vs. non-tissue. Then heatmaps can be done comparing invitro to everything else. 

	cd /home/daniel/ubuntu/workspace/all_049/gbm_049_unmerged/de/ballgown/ref_only
	mkdir stattest
	cd stattest
	


printf "\"ids\",\"type\",\"path
\"\n\"1\",\"tissue\",\"$gbm_049/expression/stringtie/ref_only/1
\"\n\"2\",\"tissue\",\"$gbm_049/expression/stringtie/ref_only/2
\"\n\"3\",\"non-tissue\",\"$gbm_049/expression/stringtie/ref_only/3
\"\n\"4\",\"non-tissue\",\"$gbm_049/expression/stringtie/ref_only/4
\"\n\"5\",\"non-tissue\",\"$gbm_049/expression/stringtie/ref_only/5
\"\n\"6\",\"non-tissue\",\"$gbm_049/expression/stringtie/ref_only/6
\"\n\"7\",\"non-tissue\",\"$gbm_049/expression/stringtie/ref_only/7
\"\n\"8\",\"non-tissue\",\"$gbm_049/expression/stringtie/ref_only/8
\"\n\"9\",\"non-tissue\",\"$gbm_049/expression/stringtie/ref_only/9
\"\n\"10\",\"non-tissue\",\"$gbm_049/expression/stringtie/ref_only/10
\"\n" > GBM049_all_stattest.csv


script:


	printf "\"ids\",\"type\",\"path\"\n\"1\",\"tissue\",\"$gbm_049/expression/stringtie/ref_only/1\"\n\"2\",\"tissue\",\"$gbm_049/expression/stringtie/ref_only/2\"\n\"3\",\"non-tissue\",\"$gbm_049/expression/stringtie/ref_only/3\"\n\"4\",\"non-tissue\",\"$gbm_049/expression/stringtie/ref_only/4\"\n\"5\",\"non-tissue\",\"$gbm_049/expression/stringtie/ref_only/5\"\n\"6\",\"non-tissue\",\"$gbm_049/expression/stringtie/ref_only/6\"\n\"7\",\"non-tissue\",\"$gbm_049/expression/stringtie/ref_only/7\"\n\"8\",\"non-tissue\",\"$gbm_049/expression/stringtie/ref_only/8\"\n\"9\",\"non-tissue\",\"$gbm_049/expression/stringtie/ref_only/9\"\n\"10\",\"non-tissue\",\"$gbm_049/expression/stringtie/ref_only/10\"\n" > GBM049_all.csv

	cat GBM049_all_stattest.csv
	
now rerun all the R scripts to see if the resulting MDS looks weird (it should) and see if stattest and heat map now works. 


	R --no-restore
	library(ballgown)
	library(genefilter)
	library(dplyr)
	library(devtools)
	library(ggplot2)
	library(gplots)
	library(GenomicRanges)

	pheno_data = read.csv("GBM049_all_stattest.csv")  

	bg = ballgown(samples=as.vector(pheno_data$path), pData=pheno_data)
	bg

	bg_table = texpr(bg, 'all')

	bg_gene_names = unique(bg_table[, 9:10])
	head(bg_gene_names)

	save(bg, file='bg.rda')

	bg

	pdf(file="GBM049_R_output_stattest.pdf")

	working_dir = "~/workspace/gbm_049_unmerged/de/ballgown/ref_only/stattest"
	setwd(working_dir)
	dir()
		
	load('bg.rda')
		
	bg_table = texpr(bg, 'all')
	bg_gene_names = unique(bg_table[, 9:10])

	gene_expression = as.data.frame(gexpr(bg))
	head(gene_expression)

	colnames(gene_expression)
	dim(gene_expression)

	i = row.names(gene_expression) == "BRD4"
	gene_expression[i,]

	transcript_gene_table = indexes(bg)$t2g
	head(transcript_gene_table)
		
	length(row.names(transcript_gene_table)) #Transcript count
	length(unique(transcript_gene_table[,"g_id"])) #Unique Gene count

	counts=table(transcript_gene_table[,"g_id"])
	c_one = length(which(counts == 1))
	c_more_than_one = length(which(counts > 1))
	c_max = max(counts)
	hist(counts, breaks=50, col="bisque4", xlab="Transcripts per gene", main="Distribution of transcript count per gene")
	legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
	legend("topright", legend_text, lty=NULL)

	full_table <- texpr(bg , 'all')
	hist(full_table$length, breaks=500, xlab="Transcript length (bp)", main="Distribution of transcript lengths", col="steelblue")
		
	min(gene_expression[,"FPKM.1"])
	max(gene_expression[,"FPKM.2"])
		
	min_nonzero=1

	data_columns=c(1:10)
	short_names=c("tissue_1","tissue_2","slice_1","slice_2","neurospheres_1", "neurospheres_2","organoid_1","organoid_2","invitro_1","invitro_2")


	#colors()
	
	data_colors=c("tomato1","tomato2","royalblue1","royalblue2","seagreen1","seagreen2","grey1","grey2","brown1","brown2")

	boxplot(log2(gene_expression[,data_columns]+min_nonzero), col=data_colors, names=short_names, las=2, ylab="log2(FPKM)", main="Distribution of FPKMs for libraries of all 10 samples")

	



	gene_expression[,"sum"]=apply(gene_expression[,data_columns], 1, sum)
	
	i = which(gene_expression[,"sum"] > 10)

	r=cor(gene_expression[i,data_columns], use="pairwise.complete.obs", method="pearson")
	r



	d=1-r
	mds=cmdscale(d, k=2, eig=TRUE)
	par(mfrow=c(1,1))
	plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all non-zero genes)", xlim=c(-0.4,0.4), ylim=c(-0.4,0.4))
	points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
	text(mds$points[,1], mds$points[,2], short_names, col=data_colors)

	plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all non-zero genes)", xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
	points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
	text(mds$points[,1], mds$points[,2], short_names, col=data_colors)



## Calculate the differential expression results including significance

	results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
	results_genes = merge(results_genes,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))


## View the distribution of differential expression values as a histogram

#Display only those that are significant according to Ballgown

	sig=which(results_genes$pval<0.05)
	results_genes[,"de"] = log2(results_genes[,"fc"])
	hist(results_genes[sig,"de"], breaks=50, col="seagreen", xlab="log2(Fold change) tissue vs non-tissue", main="Distribution of differential expression values")
	abline(v=-2, col="black", lwd=2, lty=2)
	abline(v=2, col="black", lwd=2, lty=2)
	legend("topleft", "Fold-change > 4", lwd=2, lty=2)



#Display the grand expression values from tissue and non-tissue and mark those that are significantly differentially expressed

	gene_expression[,"tissue"]=apply(gene_expression[,c(1:2)], 1, mean)
	gene_expression[,"non-tissue"]=apply(gene_expression[,c(3:10)], 1, mean)

	x=log2(gene_expression[,"invitro"]+min_nonzero)
	y=log2(gene_expression[,"non-invitro"]+min_nonzero)
	plot(x=x, y=y, pch=16, cex=0.25, xlab="tissue FPKM (log2)", ylab="non-tissue FPKM (log2)", main="tissue vs non-tissue FPKMs")
	abline(a=0, b=1)
	xsig=x[sig]
	ysig=y[sig]
	points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
	legend("topleft", "Significant", col="magenta", pch=16)



Get the gene symbols for the top N (according to corrected p-value) and display them on the plot
	
	topn = order(abs(results_genes[sig,"fc"]), decreasing=TRUE)[1:25]
	topn = order(results_genes[sig,"qval"])[1:25]
	text(x[topn], y[topn], results_genes[topn,"gene_name"], col="black", cex=0.75, srt=45)


Write a simple table of differentially expressed transcripts to an output file

Each should be significant with a log2 fold-change >= 2

	sigpi = which(results_genes[,"pval"]<0.05)
	sigp = results_genes[sigpi,]
	sigde = which(abs(sigp[,"de"]) >= 2)
	sig_tn_de = sigp[sigde,]


	o = order(sig_tn_de[,"qval"], -abs(sig_tn_de[,"de"]), decreasing=FALSE) #Order the output by or p-value and then break ties using fold-change

	output = sig_tn_de[o,c("gene_name","id","fc","pval","qval","de")]
	write.table(output, file="SigDE_R_ballgown.txt", sep="\t", row.names=FALSE, quote=FALSE)

#View selected columns of the first 25 lines of output
output[1:25,c(1,4,5)]



#### Plot #11 - Create a heatmap to vizualize expression differences between the eight samples

#Define custom dist and hclust functions for use with heatmaps

	mydist=function(c) {dist(c,method="euclidian")}
	myclust=function(c) {hclust(c,method="average")}

	main_title="sig DE Transcripts"
	par(cex.main=0.8)
	sig_genes_de=sig_tn_de[,"id"]
	sig_gene_names_de=sig_tn_de[,"gene_name"]

	data=log2(as.matrix(gene_expression[as.vector(sig_genes_de),data_columns])+1)
	heatmap.2(data, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both", margins=c(10,4), Rowv=TRUE, Colv=TRUE, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, cexRow=0.3, cexCol=1, labRow=sig_gene_names_de,col=rev(heat.colors(75)))


dev.off()



