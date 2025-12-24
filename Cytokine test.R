library(dplyr)
library(heatmaply)
library(agricolae)
library(dplyr)

## Import cytokine test data
Cytokine_human <- read_excel("~/Desktop/20230715 Cytokine test.xlsx")
Cytokine_human_nolabel <- select(Cytokine_human, - Label, - Sample, - Group) 

## Heatmap visualization
# Calculating Z scores for each cytokine
z_scores <- scale(Cytokine_human_nolabel)
Cytokine_Z <- cbind(Cytokine_human$Group, z_scores)
Cytokine_Z <- as.data.frame(Cytokine_Z)
colnames(Cytokine_Z)[1] <- c("Group")
z_scores[is.na(z_scores)] <- 0

# Heatmap
row.names(z_scores) <-  Cytokine_human$Sample
heatmaply(z_scores, na.rm = TRUE, Colv=FALSE, Rowv = FALSE)

## P value by group
Cytokine_human <- as.data.frame(Cytokine_human)
Cytokine_human$Group <- as.factor(Cytokine_human$Group)
pvalue <- c()
Star <- 5
Over <- 100
for ( i in c(Star:Over)){
  M.aov<-aov(Cytokine_human[[i]]~Group,data=Cytokine_human) 
  aov.sum<-summary(M.aov)
  pvalue[i-Star+1]<-aov.sum[[1]]$`Pr(>F)`[1]
}
pvalue <- as.data.frame(pvalue)
Cytokine_P <- t(pvalue)
colnames(Cytokine_P) <- colnames(Cytokine_human_nolabel)[2:97]
Cytokine_P <- t(Cytokine_P)


## Correlation between cFos and cytokine
# Delete sytokines with all 0
Cytokine_cFos <- na.omit(Cytokine_human)
Cytokine_cFos <- as.data.frame(Cytokine_cFos)
Cytokine_cFos <- Cytokine_cFos[,apply(Cytokine_cFos,2,function(x) !all(x==0))]

# Calculate R2 and P of simple linear regression
Cytokine_cFos_Stats <- c()
Cytokine_cFos_Stats <- as.data.frame(Cytokine_cFos_Stats)
for (i in 5:93){
  lm_fit <- lm( Cytokine_cFos[, 4] ~ Cytokine_cFos[, i], data = Cytokine_cFos)
  R2 <- summary(lm_fit)$r.squared
  Cytokine_cFos_Stats[1, i-4] <- R2
  p <- summary(lm_fit)$coefficients[2,4]
  Cytokine_cFos_Stats[2, i-4] <- p
}
colnames(Cytokine_cFos_Stats) <- colnames(Cytokine_cFos)[5:93]
row.names(Cytokine_cFos_Stats) <-  c("R2", "p")
Cytokine_cFos_Stats <- t(Cytokine_cFos_Stats)

## Combining stats
Cytokine_results <- merge(Cytokine_cFos_Stats, Cytokine_P, by = "row.names", all = T)




