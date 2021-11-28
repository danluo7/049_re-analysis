#find the fastq.gz files
#head gz file to get the platform IDs: 

	gzip -cd Sample1_Lane2_Read1.fastq.gz | head


hisat2 -p 8 --rg-id=H53NGBBXY.1 --rg SM:1 --rg LB:1_NAGTGCTT+NCACCTCA --rg PL:ILLUMINA --rg PU:H53NGBBXY.1.NAGTGCTT+NCACCTCA -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR /home/drake/workspace/049_RNAseq/FASTQ/Sample6_Lane2_Read1.fastq.gz -S /home/drake/workspace/049_reanalysis/alignments/1.sam







NAGTGCTT+NCACCTCA
CAGTGCTT+ACACCTCA
NTGATCCA+NCTCAACG
GTGATCCA+ACTCAACG
NCAGCAAG+NTCTGCAA
ACAGCAAG+GTCTGCAA
NATACGGA+NACACCAC
CATACGGA+AACACCAC
NCGTGCAT+NTGTACCA
TCGTGCAT+CTGTACCA
NTACATCC+NCACCTAG
CTATATCC+TCACCTAG
