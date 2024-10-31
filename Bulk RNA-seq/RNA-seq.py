1. Data processing
# Trim
trim_galore --gzip --clip_R1 10 --clip_R2 10 --illumina --paired ../Rawdata/${sample}/${sample}_R1.fq.gz ../Rawdata/${sample}/${sample}_R2.fq.gz

#mapping
STAR --genomeDir ~/genome/hg38/hg38.index --runThreadN 20 --readFilesIn ../Cleandata/${sample}_R1_val_1.fq.gz ../Cleandata/${sample}_R2_val_2.fq.gz --readFilesCommand zcat --outFileNamePrefix ${sample} --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 10

#Count
htseq-count -f bam -r pos --max-reads-in-buffer 1000000 --stranded no --minaqual 10 --type exon --idattr gene_id --mode union --nonunique none --secondary-alignments ignore --supplementary-alignments ignore --counts_output ${sample}_count.tsv ../mapping/${sample}Aligned.sortedByCoord.out.bam ~/data/gene/hg38.ncbiRefSeq.gtf

2. Data analysis
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import zscore, ranksums, pearsonr, spearmanr
import seaborn as sns
from seaborn import clustermap
from scipy.stats import ttest_ind
from statsmodels.sandbox.stats.multicomp import multipletests as FDR


# TPM matrix
data1 = data*1e6/np.sum(data, axis=0)

# merged data
stage1 = np.array(['U87-DMSO', 'U87-O', 'U87-24'])
meta1 = np.array(['U87-DMSO', 'U87-DMSO', 'U87-O', 'U87-O','U87-24', 'U87-24'])

data_merge = pd.DataFrame(index = data1.index, columns = stage1)
for x in stage1:
    data_merge[x] = np.mean(data1.iloc[:,meta1==x], axis=1)


#Whole transcriptome heatmap
cg = clustermap(zscore(data1, axis=1), cmap='seismic', xticklabels=data1.columns, yticklabels=[], metric='euclidean', method='ward', cbar_kws={'label': 'Normalized Gene Expression'},figsize=(3,4))

#DEG analysis and violin plot
pval = np.array([ttest_ind(data1[['U87-DMSO-1', 'U87-DMSO-2']].iloc[i], data1[['U87-O-1', 'U87-O-2']].iloc[i]) for i in range(len(data1))])
fdr = FDR1(pval[:,1])
fc = np.log2(data_merge['U87-O']+1)-np.log2(data_merge['U87-DMSO']+1)


fig, ax = plt.subplots(figsize = (4,3))
DMfilter = np.logical_and(pval[:,1]<0.05, abs(fc)>np.log2(1.5))
ax.scatter(fc[~DMfilter], -np.log2(pval[:,1])[~DMfilter], c='grey', s=8, edgecolors='none', alpha=0.8)
DMfilter = np.logical_and(pval[:,1]<0.05, fc>np.log2(1.5))
ax.scatter(fc[DMfilter], -np.log2(pval[:,1])[DMfilter], c='r', s=8, edgecolors='none', alpha=0.8)
ax.text(-4.5,2,'Up: ' + str(np.sum(DMfilter)), c='r', fontsize=12)
DMfilter = np.logical_and(pval[:,1]<0.05, fc<-np.log2(1.5))
ax.scatter(fc[DMfilter], -np.log2(pval[:,1])[DMfilter], c='blue', s=8, edgecolors='none', alpha=0.8)
ax.text(-4.5,0.5,'Down: ' + str(np.sum(DMfilter)), c='blue', fontsize=12)
ax.plot([np.log2(1.5),np.log2(1.5)],[0,20], c='red', linestyle='dashed', linewidth=0.5)
ax.plot([-np.log2(1.5),-np.log2(1.5)],[0,20], c='red', linestyle='dashed', linewidth=0.5)
ax.plot([-20,20],[-np.log2(0.05),-np.log2(0.05)], c='red', linestyle='dashed', linewidth=0.5)
ax.set_xlabel('log2(FC)')
ax.set_ylabel('-log2(p-value)')
ax.set_xlim([-5,5])
ax.set_title('Osimertinib v.s. DMSO')
plt.tight_layout()
plt.savefig('volcano.pdf', transparent=True)
plt.close()

#pathway gene expression heatmap
genes1 = np.array(['HSPA5', 'XBP1', 'ERN1', 'EIF2AK3', 'CANX', 'DDIT3', 'ATF4', 'ATF6', 'HERPUD1', 'DNAJC3', 'PPP1R15A'])
genes2 = np.array(['BAX', 'TP53', 'CASP1','CASP4','CASP8', 'CASP9','CASP10', 'FAS',  'TNFRSF10B', 'CYCS', 'PMAIP1', 'HRK', 'PTEN', 'GADD45A'])
genes3 = np.array(['MKI67', 'CCND1','CCND2','CCND3',  'CCNE2','CDK2', 'CDK6','CDK10', 'CDK15','MTOR', 'JAK2', 'E2F6','E2F8', 'ORC1'])
genes4 = np.array(['ATR', 'RAD51','DDB1','DDB2', 'CHEK1', 'MSH2', 'PCNA', 'RPA1', 'SMC3', 'DCLRE1B'])
genes5 = np.array(['IL1A', 'IL1B', 'IL6', 'IL33', 'LIF', 'TGFB1', 'TGFB2', 'TGFB3','IL11',  'CCL2'])

cg = clustermap(zscore(data1.loc[genes1], axis=1), cmap='seismic', xticklabels=data1.columns,row_cluster=False, col_cluster=False, yticklabels=genes1, metric='euclidean', method='ward', cbar_kws={'label': 'Normalized\nGene Expression'},figsize=(3,4), vmin = -1.5, vmax =1.5)















