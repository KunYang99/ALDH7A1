library("DESeq2")

pheno <- read.csv('data/pheno.csv')
edata <- read.csv('data/count_tb.csv', row.names=1)
edata <- edata[, as.vector(pheno$id)]
table(names(edata) == pheno$id)

# run DESeq
dds <- DESeqDataSetFromMatrix(countData = edata, colData = pheno, design= ~ group)
dds <- DESeq(dds)
re <- results(dds, contrast=c("group","KO","WT"), alpha=0.05)
re <- re[order(re$pvalue),]
sig <- subset(re, padj < 0.05)
print(paste('all:', nrow(sig), 'sig genes found'))
write.csv(re, 're/DGE_KOvsWT.csv')
write.csv(sig, 're/sig_KOvsWT.csv')

