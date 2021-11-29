#find the fastq.gz files
#head gz file to get the platform IDs: 

	gzip -cd Sample1_Lane2_Read1.fastq.gz | head

hisat2 -p 8 --rg-id=H53NGBBXY.1 --rg SM:1 --rg LB:1_NAGTGCTT+NCACCTCA --rg PL:ILLUMINA --rg PU:H53NGBBXY.1.NAGTGCTT+NCACCTCA -x $gbm/RNA_REF_FA/hg38/genome --dta -U /home/floyd/ubuntu/workspace/gbm/049/FASTQ/Sample6_Lane2_Read1.fastq.gz -S /home/floyd/ubuntu/workspace/gbm/049/alignments/1.sam

hisat2 -p 8 --rg-id=H53NGBBXY.1 --rg SM:1 --rg LB:1_NAGTGCTT+NCACCTCA --rg PL:ILLUMINA --rg PU:H53NGBBXY.1.NAGTGCTT+NCACCTCA -x $gbm/RNA_REF_FA/hg38/genome --dta -U /home/daniel/ubuntu/workspace/all_049/gbm_049_unmerged/FASTQ/Sample6_Lane2_Read1.fastq.gz -S /home/daniel/ubuntu/workspace/all_049/gbm_049_unmerged/alignments/1.sam




hisat2 -p 8 --rg-id=H53NGBBXY.1 --rg SM:1 --rg LB:1_NAGTGCTT+NCACCTCA --rg PL:ILLUMINA --rg PU:H53NGBBXY.1.NAGTGCTT+NCACCTCA -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR /home/floyd/ubuntu/workspace/gbm/049/FASTQ/Sample6_Lane2_Read1.fastq.gz -S /home/floyd/ubuntu/workspace/gbm/049/alignments/1.sam

hisat2 -p 8 --rg-id=H53NGBBXY.1 --rg SM:1 --rg LB:1_CAGTGCTT+ACACCTCA --rg PL:ILLUMINA --rg PU:H53NGBBXY.1.CAGTGCTT+ACACCTCA -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR /home/floyd/ubuntu/workspace/gbm/049/FASTQ/Sample6_Lane3_Read1.fastq.gz -S /home/floyd/ubuntu/workspace/gbm/049/alignments/2.sam




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
