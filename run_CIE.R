# Please cite the following two papers if using CIE:
# 1) Farahman, S., O'Connor, C., Macoska, J., Zarringhalam, K. 'Causal Inference Engine: A platform for directional gene set enrichment analysis and inference of active transcriptional regulators', NAR doi: https://doi.org/10.1093/nar/gkz1046
# 2) Fakhry CT, Choudhary P, Gutteridge A, Sidders B, Chen P, Ziemek D, Zarringhalam, K. 'Interpreting transcriptional changes using causal graphs: new methods and their practical utility on public networks'. BMC Bioinformatics. 2016;17(1):318

library(CIE)
library(org.Hs.eg.db)
source('CIE/bin/mycie.R')

#database tissueCorrectedChIP, can be downloaded from https://umbibio.math.umb.edu/cie/app
tissChIP  <- filterChIPAtlas(NA, NA, NA, cellLineType="all", tissueCorrect=TRUE, databaseDir="/CIE/data/human/tissueCorrectedChIP/")

# path = the folder where you store the DESeq2 analysis results
path <- 'deseq_re'
for(f in list.files(path, '.csv')){
  dge <- read.csv(file.path(path, f))
  dge <- subset(dge, pvalue <= 0.05) 
  names(dge)[1] <- 'name'
  tmp <- data.frame(name=dge$name, pval=dge$pvalue, fc=dge$log2FoldChange)

  #entrez_human.csv, can be downloded from Ensembl BioMart. It has two columns: entrez and name, where the entrez column contains Entrez IDs and the name column contains gene symbols
  genemap <- read.csv('CIE/data/entrez_human.csv')
  DGEs <- merge(genemap, tmp, by='name')
  print(paste('working on', f, '--', nrow(DGEs)))

  re <- runCIE(DGEs=DGEs,
               methods='Ternary',
               ents=tissChIP$ChIPall.ents,
               rels=tissChIP$ChIPall.rels,
               useFile=F,
               hypTabs="1")

  write.csv(re, paste0('re/cie_ChIP_Ternary_', f), row.names=F)
}
