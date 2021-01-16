

#########################################
#########################################
#anno <- read.delim("../Desktop/stud/7sem/kursinis2/ewas/EWAS_Atlas_probe_annotations.tsv", sep = "\t", header = TRUE)
#anno <- remove.factors(anno)

#length(anno$Probe.id)
#length(unique(anno$Probe.id))
data <- readRDS("../Desktop/stud/7sem/kursinis2/TF_cg_study.rds")
anno <- read.delim("../Desktop/stud/7sem/kursinis2/HM450.hg19.manifest.original.tsv", as.is=TRUE)
length(anno$Name)  # CpG id
#anno$Name[1:5]

anno2 <- anno 
length(anno$Name[grepl("cg\\d{8}", anno$Name)])
anno <- anno[grepl("cg\\d{8}", anno$Name),]

#str(anno2$Probe.id)

#rasta 3k kurie nera cg
length(anno$Name[!grepl("cg\\d{8}", anno$Name)])
#anno$Name[!grepl("cg\\d{8}", anno$Name)][1:10]

col <- colnames(data)# [1:5]
#odds  <- matrix(nrow = length(unique(col)), ncol = length(unique(col)),
#                       dimnames = list(unique(col),  unique(col)))
#pvals <- matrix(nrow = length(unique(col)), ncol = length(unique(col)),
#                dimnames = list(unique(col),  unique(col)))

odds <- readRDS("../Desktop/stud/7sem/kursinis2/unadjusted_odds_5.rds")
pvals <- readRDS("../Desktop/stud/7sem/kursinis2/unadjusted_pvals_5.rds")


#length(data[data[,1] == TRUE])

#unique(anno2$Probe.id %in% rownames(data)[data[,1]])
#col1 <- anno2$Probe.id %in% rownames(data)[data[,1]]
#col2 <- anno2$Probe.id %in% rownames(data)[data[,3]]
#tab <- table(col1, col2)
#tab
#sum(data[,3])
#colnames(data)[1]
#rownames(data)[data[,2] == TRUE]
#fish <- fisher.test(tab)
#fish$p.value
#fish$estimate

Sys.time()
n1 <- 612
n2 <- 232
#per 6 valandas
for(i in n1:length(col)){
  col1 <- anno$Name %in% rownames(data)[data[,i]]
  col1 <- factor(col1, levels = c("TRUE", "FALSE"))
  for(j in i:length(col)){
    col2 <- anno$Name %in% rownames(data)[data[,j]]
    col2 <- factor(col2, levels = c("TRUE", "FALSE"))
    tab  <- table(col1, col2)
    fish <- fisher.test(tab)
    
    odds[i,j]  <- fish$estimate
    pvals[i,j] <- fish$p.value
  }
}
Sys.time()

#tab
#i
#j
#col1<- factor(col1, levels = c("TRUE", "FALSE"))
#table(col1,factor(col2, levels = c("TRUE", "FALSE")))
#odds
saveRDS(odds, "../Desktop/stud/7sem/kursinis2/unadjusted_odds_6.rds")
#pvals
saveRDS(pvals, "../Desktop/stud/7sem/kursinis2/unadjusted_pvals_6.rds")
pvals2 <- pvals
pvals3 <- pvals

Sys.time()
pvals2[] <- p.adjust(pvals, method = "fdr") < 0.05 & odds > 1 
#pvals2
Sys.time()
saveRDS(pvals2, "../Desktop/stud/7sem/kursinis2/adjusted_pval_more_1_odds_6.rds")
Sys.time()

pvals3[] <- p.adjust(pvals, method = "fdr") < 0.05 & odds < 1 
#pvals3
Sys.time()
saveRDS(pvals3, "../Desktop/stud/7sem/kursinis2/adjusted_pval_less_1_odds_6.rds")


cat("finish", '\n')
Sys.time()

