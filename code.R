require('KnowSeq')
require(caret)
require(dplyr)
require(ROCR)
require(limma)
require(readxl)
require(GEOquery)
require(class)
set.seed(50)

# PREPROCESs.

# GSE156063
GSE156063 <- getGEO("GSE156063", destdir = 'GSE156063')
GSE156063_age <- round(as.numeric(GSE156063$GSE156063_series_matrix.txt.gz$`age:ch1`))
GSE156063_gender <- GSE156063$GSE156063_series_matrix.txt.gz$`gender:ch1`
GSE156063_rpm <- GSE156063$GSE156063_series_matrix.txt.gz$`sars-cov-2 pcr:ch1`
GSE156063_pcr <- GSE156063$GSE156063_series_matrix.txt.gz$`sars-cov-2 rpm:ch1`
GSE156063_labels <- GSE156063$GSE156063_series_matrix.txt.gz$`disease state:ch1`  # 
GSE156063_counts <- read.csv('GSE156063_swab_gene_counts.csv')
rownames <- GSE156063_counts[,1]
GSE156063_counts[,1] <- NULL
rownames(GSE156063_counts) <- rownames
Annotation_gene_GSE156063 <- getGenesAnnotation(rownames(GSE156063_counts))
Annotation_gene_GSE156063 <- Annotation_gene_GSE156063[order(Annotation_gene_GSE156063$ensembl_gene_id),]
GSE156063_expressionMatrix <- calculateGeneExpressionValues(as.matrix(GSE156063_counts), Annotation_gene_GSE156063, genesNames = TRUE) # 
GSE156063_severity <- read_excel('Patient_class_inpatient_outpatient_MickKammetal.xlsx')
GSE156063_outliers <- RNAseqQA(GSE156063_expressionMatrix,toPNG = FALSE, toPDF = FALSE,toRemoval = TRUE)
GSE156063_age <- GSE156063_age[-which(colnames(GSE156063_expressionMatrix)%in%GSE156063_outliers$outliers)]
GSE156063_gender <- GSE156063_gender[-which(colnames(GSE156063_expressionMatrix)%in%GSE156063_outliers$outliers)]
GSE156063_rpm <- GSE156063_rpm[-which(colnames(GSE156063_expressionMatrix)%in%GSE156063_outliers$outliers)]
GSE156063_pcr <- GSE156063_pcr[-which(colnames(GSE156063_expressionMatrix)%in%GSE156063_outliers$outliers)]
GSE156063_labels <- GSE156063_labels[-which(colnames(GSE156063_expressionMatrix)%in%GSE156063_outliers$outliers)] 
GSE156063_severity <-  GSE156063_severity[-which(GSE156063_severity$CZB_ID%in%GSE156063_outliers$outliers),]
GSE156063_expressionMatrix <- GSE156063_outliers$matrix
GSE156063_severity <- GSE156063_severity[which(GSE156063_severity$Viral_status=='SARS-CoV-2'),]
GSE156063_severity <- GSE156063_severity[which(GSE156063_severity$Patient_class=='Outpatient'|GSE156063_severity$Patient_class=='Inpatient'|GSE156063_severity$Patient_class=='Emergency'),]
GSE156063_age <- GSE156063_age[which(colnames(GSE156063_expressionMatrix)%in%GSE156063_severity$CZB_ID)]
GSE156063_gender <- GSE156063_gender[which(colnames(GSE156063_expressionMatrix)%in%GSE156063_severity$CZB_ID)]
GSE156063_rpm <- GSE156063_rpm[which(colnames(GSE156063_expressionMatrix)%in%GSE156063_severity$CZB_ID)]
GSE156063_pcr <- GSE156063_pcr[which(colnames(GSE156063_expressionMatrix)%in%GSE156063_severity$CZB_ID)]
GSE156063_labels <- GSE156063_labels[which(colnames(GSE156063_expressionMatrix)%in%GSE156063_severity$CZB_ID)] 
GSE156063_expressionMatrix <- GSE156063_expressionMatrix[,which(colnames(GSE156063_expressionMatrix)%in%GSE156063_severity$CZB_ID)]
GSE156063_age <- GSE156063_age[order(colnames(GSE156063_expressionMatrix))]
GSE156063_gender <- GSE156063_gender[order(colnames(GSE156063_expressionMatrix))]
GSE156063_rpm <- GSE156063_rpm[order(colnames(GSE156063_expressionMatrix))]
GSE156063_pcr <- GSE156063_pcr[order(colnames(GSE156063_expressionMatrix))]
GSE156063_expressionMatrix <- GSE156063_expressionMatrix[,order(colnames(GSE156063_expressionMatrix))]
GSE156063_severity <- GSE156063_severity[order(GSE156063_severity$CZB_ID),]
GSE156063_labels_severity <- GSE156063_severity$Patient_class
GSE156063_labels_severity[which(GSE156063_severity$ICU=='ICU')]<- 'ICU'
GSE156063_expressionMatrix <- GSE156063_expressionMatrix[,-which(GSE156063_labels_severity=='Emergency')]
GSE156063_age <- GSE156063_age[-which(GSE156063_labels_severity=='Emergency')]
GSE156063_gender <- GSE156063_gender[-which(GSE156063_labels_severity=='Emergency')]
GSE156063_rpm <- GSE156063_rpm[-which(GSE156063_labels_severity=='Emergency')]
GSE156063_pcr <- GSE156063_pcr[-which(GSE156063_labels_severity=='Emergency')]
GSE156063_labels_severity <- GSE156063_labels_severity[-which(GSE156063_labels_severity=='Emergency')]

# GSE162835 
GSE162835 <- getGEO("GSE162835", destdir = 'GSE162835')
GSE162835_sup <- as.data.frame(read_excel('sup_info.xlsx'))
GSE162835_age <- GSE162835_sup$...2[3:52]
GSE162835_gender <- GSE162835_sup$...3[3:52]
GSE162835_labels <- GSE162835$GSE162835_series_matrix.txt.gz$`disease:ch1`
for (i in 1:50){
  GSE162835_labels[i] <- 'SC2' 
}
GSE162835_labels_severity <- GSE162835$GSE162835_series_matrix.txt.gz$`disease severity:ch1`
for (i in 1:length(GSE162835_labels_severity)){
  if (GSE162835_labels_severity[i]=='Asymptomatic/Mild'){
    GSE162835_labels_severity[i] <- 'Outpatient'
  } else if (GSE162835_labels_severity[i]=='Moderate'){
    GSE162835_labels_severity[i] <- 'Inpatient'
  } else if (GSE162835_labels_severity[i]=='Severe/Critical'){
    GSE162835_labels_severity[i] <- 'ICU'
  }
}
GSE162835_expressionMatrix <- as.matrix(read_excel('GSE162835_COVID_GEO_processed.xlsx'))
rownames <- GSE162835_expressionMatrix[,1]
rownames(GSE162835_expressionMatrix) <- rownames
GSE162835_expressionMatrix <- GSE162835_expressionMatrix[,2:51]
GSE162835_expressionMatrix <- apply(GSE162835_expressionMatrix,2,as.numeric)
rownames(GSE162835_expressionMatrix) <- rownames
GSE162835_outliers <- RNAseqQA(GSE162835_expressionMatrix,toPNG = FALSE, toPDF = FALSE,toRemoval = TRUE) 
GSE162835_labels <- GSE162835_labels[-which(colnames(GSE162835_expressionMatrix) %in% GSE162835_outliers$outliers)]
GSE162835_labels_severity <- GSE162835_labels_severity[-which(colnames(GSE162835_expressionMatrix) %in% GSE162835_outliers$outliers)]
GSE162835_expressionMatrix <- GSE162835_outliers$matrix
GSE162835_age <- GSE162835_age[-c(48,50)]
GSE162835_gender <- GSE162835_gender[-c(48,50)]

# GSE152075
GSE152075 <- getGEO("GSE152075", destdir = 'GSE152075')
GSE152075_age <- GSE152075$GSE152075_series_matrix.txt.gz$`age:ch1`
GSE152075_gender <- GSE152075$GSE152075_series_matrix.txt.gz$`gender:ch1`
GSE152075_labels <- GSE152075$GSE152075_series_matrix.txt.gz$`sars-cov-2 positivity:ch1`
for (i in 1:length(GSE152075_labels)){
  if (GSE152075_labels[i] == 'pos'){
    GSE152075_labels[i] <- 'SC2'
  } else {
    GSE152075_labels[i] <- 'Control'
  }
}
GSE152075_counts <- as.matrix(read.table('GSE152075_raw_counts_GEO.txt', header =TRUE))
Annotation_gene_GSE152075 <- getGenesAnnotation(rownames(GSE152075_counts), filter = 'external_gene_name')
Annotation_gene_GSE152075 <- Annotation_gene_GSE152075[order(Annotation_gene_GSE152075$external_gene_name),]
GSE152075_counts<- GSE152075_counts[order(rownames(GSE152075_counts)),]
GSE152075_counts_1 <- GSE152075_counts[which(rownames(GSE152075_counts) %in% Annotation_gene_GSE152075[,2]),]
for (i in 1:length(rownames(GSE152075_counts_1))){
  rownames(GSE152075_counts_1)[i] <- Annotation_gene_GSE152075[which(Annotation_gene_GSE152075[,2] == rownames(GSE152075_counts_1)[i])[1],1] 
}
Annotation_gene_GSE152075_1 <- getGenesAnnotation(rownames(GSE152075_counts_1))
Annotation_gene_GSE152075_1 <- Annotation_gene_GSE152075_1[order(Annotation_gene_GSE152075_1$ensembl_gene_id),]
GSE152075_counts_1<- GSE152075_counts_1[order(rownames(GSE152075_counts_1)),]
GSE152075_expressionMatrix <- calculateGeneExpressionValues(GSE152075_counts_1, Annotation_gene_GSE152075_1, genesNames = TRUE) #
GSE152075_severity <- read.csv('2021-03-19_Rojas.csv')
GSE152075_outliers <- RNAseqQA(GSE152075_expressionMatrix,toPNG = FALSE, toPDF = FALSE,toRemoval = TRUE) # 
GSE152075_age <- GSE152075_age[-which(colnames(GSE152075_expressionMatrix)%in%GSE152075_outliers$outliers)]
GSE152075_gender <- GSE152075_gender[-which(colnames(GSE152075_expressionMatrix)%in%GSE152075_outliers$outliers)]
GSE152075_labels <- GSE152075_labels[-which(colnames(GSE152075_expressionMatrix)%in%GSE152075_outliers$outliers)] 
GSE152075_severity <-  GSE152075_severity[-which(GSE152075_severity$ï..alt_name%in%GSE152075_outliers$outliers),]
GSE152075_expressionMatrix <- GSE152075_outliers$matrix
GSE152075_severity <- GSE152075_severity[which(GSE152075_severity$covid_status=='pos'),]
GSE152075_severity <- GSE152075_severity[which(GSE152075_severity$admitted.to.hospital..not.just.ED..at.time.of.initial.test.=='yes'|GSE152075_severity$admitted.to.hospital..not.just.ED..at.time.of.initial.test.=='no'),]
GSE152075_expressionMatrix_severity <- GSE152075_expressionMatrix[,which(colnames(GSE152075_expressionMatrix)%in%GSE152075_severity$ï..alt_name)]
GSE152075_expressionMatrix_severity <- cbind(GSE152075_expressionMatrix_severity,GSE152075_expressionMatrix[,which(GSE152075_labels=='Control')])
GSE152075_labels_severity <- c(GSE152075_severity$admitted.to.hospital..not.just.ED..at.time.of.initial.test., rep('Control',52))
for (i in 1:length(GSE152075_labels_severity)){
  if (GSE152075_labels_severity[i]=='no'){
    GSE152075_labels_severity[i]<-'Outpatient'
  } else if (GSE152075_labels_severity[i]=='yes'){
    GSE152075_labels_severity[i]<-'Inpatient'
  }
}
GSE152075_labels_severity[which(GSE152075_severity$ICU=='yes')]<- 'ICU'


# GSE152075 52 CONTROL / 4 UCI / 5 INPATIENT / 43 OUTPATIENT
# GSE162835 3 UCI / 10 INPATIENT / 35 OUTPATIENT
# GSE156063 4 UCI / 4 INPATIENT / 52 OUTPATIENT 

I1 <- intersect(rownames(GSE156063_expressionMatrix), rownames(GSE152075_expressionMatrix_severity))
I2 <- intersect(I1, rownames(GSE162835_expressionMatrix))
GSE156063_I <- GSE156063_expressionMatrix[which(rownames(GSE156063_expressionMatrix) %in%  I2),]
GSE152075_I <- GSE152075_expressionMatrix_severity[which(rownames(GSE152075_expressionMatrix_severity) %in%  I2),]
GSE162835_I <- GSE162835_expressionMatrix[which(rownames(GSE162835_expressionMatrix) %in%  I2),]
GSE156063_I <- GSE156063_I[order(rownames(GSE156063_I)),] 
GSE152075_I <- GSE152075_I[order(rownames(GSE152075_I)),] 
labels <- c(GSE156063_labels_severity,GSE152075_labels_severity,GSE162835_labels_severity)
expression_matrix <- cbind(GSE156063_I,GSE152075_I,GSE162835_I)

#NORMALIZATION BETWEEN ARRAYS
expression_matrix_norm_scale <- normalizeBetweenArrays(expression_matrix, method = 'scale')

# BACTH EFFECT
expression_matrix_norm_scale_fix <- batchEffectRemoval(expression_matrix_norm_scale,labels, method = 'sva')

#OUTLIERS REMOVAL
outliers_scale <- RNAseqQA(expression_matrix_norm_scale_fix,toRemoval = TRUE, toPNG = FALSE, toPDF = FALSE) #3
expression_matrix_norm_scale_fix_out <- expression_matrix_norm_scale_fix[,-which(colnames(expression_matrix_norm_scale_fix) %in% outliers_scale$outliers)]
labels_scale <- labels[-which(colnames(expression_matrix_norm_scale_fix) %in% outliers_scale$outliers)]


#train-test
Index_train_test_scale <- createDataPartition(labels_scale, p = .80, list = FALSE, times = 1)
train_labels_scale <- labels_scale[Index_train_test_scale]
test_labels_scale <- labels_scale[-Index_train_test_scale]
train_matrix_scale <- expression_matrix_norm_scale_fix_out[,Index_train_test_scale]
test_matrix_scale <- expression_matrix_norm_scale_fix_out[,-Index_train_test_scale]

#DEGS
folds <-5
cvIndex_scale <- createDataPartition(train_labels_scale, p = .80, list = FALSE, times = folds)
cvResults_scale <- list()
cvDEGs_scale <- list ()
for (i in seq_len(folds)){
  cvResults_scale[[i]] <- DEGsExtraction(train_matrix_scale[,cvIndex_scale[,i]], as.factor(train_labels_scale[cvIndex_scale[,i]]), lfc=1, cov=2, pvalue = 0.05, number = Inf)
  cvDEGs_scale[[i]] <- rownames(cvResults_scale[[i]]$DEGsMatrix)
}
DEGs_scale <- Reduce(f='intersect', cvDEGs_scale) # lfc 1 cov 2 pvalue 0.05 #136 genes

#feature selection algorithms
gene_mrmr_scale <- featureSelection(t(train_matrix_scale),train_labels_scale,DEGs_scale, mode ='mrmr')


#k-nn
set.seed(200)
knn_train_mrmr_scale <- knn_trn(t(train_matrix_scale), as.factor(train_labels_scale), names(gene_mrmr_scale), 5) 
knn_test_mrmr_scale <- knn_test(t(train_matrix_scale), as.factor(train_labels_scale), t(test_matrix_scale), as.factor(test_labels_scale), names(gene_mrmr_scale), bestK =  knn_train_mrmr_scale$bestK)

#validation plot
plot(knn_train_mrmr_scale$accuracyInfo$meanAccuracy[1:15], type = 'l', col= 'black', ylab='Metric Performance', xlab='Genes', lwd=2, ylim = c(0.89,1), panel.first = grid(col='gray45'),cex.axis=1.2,cex.lab=1.2)
lines(knn_train_mrmr_scale$sensitivityInfo$meanSensitivity[1:15], col='blue', lwd=2, lty=2)
lines(knn_train_mrmr_scale$specificityInfo$meanSpecificity[1:15], col='#FF8B00', lwd=2, lty=4)
lines(knn_train_mrmr_scale$F1Info$meanF1[1:15], col='red', lwd=2, lty=4)
legend(x=11.9 ,y =0.9137, c('Accuracy', 'Sensitivity','Specificity','F1-Score'), lty = c(1,2,4,5), col = c('black','blue','#FF8B00','red'), cex=1.2)


#LOOCV
LOOCV <- knn.cv(t(expression_matrix_norm_scale_fix_out[which(rownames(expression_matrix_norm_scale_fix_out) %in% names(gene_mrmr_scale)[1:4]),]), cl = labels_scale, k=7, use.all = TRUE)
confusionMatrix(data = LOOCV, reference = as.factor(labels_scale))
dataPlot(confusionMatrix(data = LOOCV, reference = as.factor(labels_scale))$table, labels = labels_scale  ,mode = "confusionMatrix",toPNG = FALSE, toPDF = FALSE)



# AUC 
response <- as.factor(test_labels_scale)
aucs <- rep(NA, length(levels(response))) # store AUCs
legendLabels <- as.character()
colours <- c('red','blue','green','black')

par(oma = c(5, 1, 0, 1))
plot(x=NA, y=NA, xlim=c(0,1), ylim=c(0,1),
     ylab="Sensitivity",
     xlab="1 - Specificity",
     bty='n',
     cex.axis=1.3,
     cex.lab=1.3)

for (i in seq_along(levels(response))) {
  cur.class <- levels(response)[i]
  binaryTraining.labels <- as.factor(train_labels_scale == cur.class)
  binaryTest.labels <- as.factor(test_labels_scale == cur.class)
  
  binary_knn_cv_mrmr_results <- knn_trn(t(train_matrix_scale), binaryTraining.labels, names(gene_mrmr_scale)[1:4])
  
  binary_knn_test_mrmr_results <- knn_test(t(train_matrix_scale), binaryTraining.labels, t(test_matrix_scale), binaryTest.labels, names(gene_mrmr_scale)[1:4], bestK = binary_knn_cv_mrmr_results$bestK)
  
  score <- binary_knn_test_mrmr_results$predictions[[4]]
  score <- as.vector(score)
  score[score=='FALSE'] <- 0
  score[score=='TRUE'] <- 1
  binaryTest.labels <- as.vector(binaryTest.labels)
  binaryTest.labels[binaryTest.labels=='FALSE'] <- 0
  binaryTest.labels[binaryTest.labels=='TRUE'] <- 1
  pred <- prediction(as.numeric(score), as.numeric(binaryTest.labels))
  perf <- performance(pred, "tpr", "fpr")
  roc.x <- unlist(perf@x.values)
  roc.y <- unlist(perf@y.values)
  lines(roc.y ~ roc.x, col = colours[i], lwd = 2)
  # store AUC
  auc <- performance(pred, "auc")
  auc <- unlist(slot(auc, "y.values"))
  aucs[i] <- auc
  legendLabels[i] <- paste(levels(response)[i], " AUC: ",format(round(aucs[i], 4), nsmall = 3),sep = "")
}

print(paste0("Mean AUC under the precision-recall curve is: ", round(mean(aucs), 2)))

lines(x=c(0,1), c(0,1))
legend(x=0.3951 ,y =0.305, legendLabels, lty=1, ncol= 1,inset = c(0,0),  col = colours, cex = 1.3,lwd=3)

par(fig = c(0, 1, 0, 1), oma = c(0.6, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legendLabels, lty=1, ncol= 1,inset = c(0,0),  col = colours, cex = 1.1,lwd=3)



#CLUSTERING K-MEANS
set.seed(50)
# los 7 primeros genes del ranking
fviz_nbclust(t(expression_matrix_norm_scale_fix_out[which(rownames(expression_matrix_norm_scale_fix_out) %in% names(gene_mrmr_scale)[1:7]),]), kmeans, method = "wss")
#todos los DEGS
fviz_nbclust(t(expression_matrix_norm_scale_fix_out[which(rownames(expression_matrix_norm_scale_fix_out) %in% names(gene_mrmr_scale)),]), kmeans, method = "wss")
# los 4 primeros
fviz_nbclust(t(expression_matrix_norm_scale_fix_out[which(rownames(expression_matrix_norm_scale_fix_out) %in% names(gene_mrmr_scale)[1:4]),]), kmeans, method = "wss")



# K-means clustering con K = 4 y 20 asignaciones aleatorias de clústeres iniciales 
k.means <- kmeans(x = t(expression_matrix_norm_scale_fix_out[which(rownames(expression_matrix_norm_scale_fix_out) %in% names(gene_mrmr_scale)[1:7]),]), centers = 4, nstart = 25)
fviz_cluster(k.means, data = t(expression_matrix_norm_scale_fix_out[which(rownames(expression_matrix_norm_scale_fix_out) %in% names(gene_mrmr_scale)[1:7]),]),
             palette = c("#2E9FDF", "#00AFBB", "#E7B800",'red'), 
             geom = "point",
             ellipse.type = "euclid", 
             ggtheme = theme_bw()
)

# Dimension reduction using PCA
pca <- prcomp(t(expression_matrix_norm_scale_fix_out[which(rownames(expression_matrix_norm_scale_fix_out) %in% names(gene_mrmr_scale)[1:7]),]),  scale =TRUE)
# Coordinates of individuals
ind.coord <- as.data.frame(get_pca_ind(pca)$coord)
# Add clusters obtained using the K-means algorithm
ind.coord$cluster <- factor(k.means$cluster)
# Add Species groups from the original data sett
ind.coord$label <- labels_scale
# Data inspection
head(ind.coord)
# Percentage of variance explained by dimensions
eigenvalue <- round(get_eigenvalue(pca), 1)
variance.percent <- eigenvalue$variance.percent
head(eigenvalue)
ggscatter(
  ind.coord, x = "Dim.1", y = "Dim.2", 
  color = "cluster", palette = c("#2E9FDF", "#00AFBB", "#E7B800",'red'), ellipse = TRUE, ellipse.type = "euclid", shape = 'label',
  size = 2,  legend = "right", ggtheme = theme_light(),
  xlab = paste0("Dim 1 (", variance.percent[1], "% )" ),
  ylab = paste0("Dim 2 (", variance.percent[2], "% )" )
) + stat_mean(aes(color = cluster), size = 4) + coord_fixed(ratio = 1)



#heatmap
library("gplots")
library("heatmap.plus")
library("RColorBrewer")


condition_colors <- unlist(lapply(labels_scale,function(x){
  if(grepl("ICU",x)) 'blue' #pink
  else if (grepl('Control',x)) 'red' #grey
  else if (grepl('Inpatient',x)) 'green'
  else if (grepl('Outpatient',x)) 'black'
}))

# I like to have a line just to assign input in one step
input <- as.matrix(expression_matrix_norm_scale_fix_out[which(rownames(expression_matrix_norm_scale_fix_out) %in% names(gene_mrmr_scale)),])
par(oma = c(5, 1, 0, 1))
heatmap.2(input, trace="none", density="none", col=bluered(20), cexRow=1, cexCol=0.2,
          ColSideColors=condition_colors, scale="row",
          hclust=function(x) hclust(x))
par(fig = c(0, 1, 0, 1), oma = c(27.9, 35.5, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
labelegend <- c('ICU','Control','Inpatient','Outpatient')
legend("bottom", labelegend, ncol= 1,inset = .02, fill = c('blue','red','green','black'), cex = 1.1)


