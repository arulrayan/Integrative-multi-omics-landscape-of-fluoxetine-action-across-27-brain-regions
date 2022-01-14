
 library("DESeq2")
library(stringr ) ;

dirs = c("PLC" , "ILC" , "CGC" , "InsularCortex" , "NaCshell" , "NaCcore" , "LS" , "MS" , "BNST" , "MPOA" , "PVN" , "MDT" , "LH" , "VMH" , "BLA" , "CMA" , "dorDG" , "dorCA1" , "HB" , "venDG" , "venCA1" , "Arc" , "VTA" , "PAG" , "Raphe" , "LC" , "Snt");

chdir = getwd() ;


for (tis in 1:27) 
{
setwd(dirs[tis])
scondition = read.table("sampletype", sep="\t") 
numsample = dim(scondition)[1] ;

chdata <- matrix(scan("correctedScore", 0), ncol= numsample, byrow=TRUE)


chipname = as.matrix(read.table("chipnames")) ;
pos = which(str_detect(chipname , "Rep3" ) == TRUE  ) ;
chdata = chdata[, -pos] ;
scondition = scondition[-pos,1] ;
numsample= dim(chdata)[2] ;


maxch = apply(chdata , 1, 'max') ;
maxch90 = quantile(maxch, probs = 0.90) ;
rpos = which(maxch > maxch90) ;


colnames(chdata)<-c(as.character(1: numsample))
mat1 = as.character(1:numsample) ;
mat2 = as.matrix(scondition)

libType<-rep("single-end", dim(chdata)[2])

 Design = data.frame(  row.names = mat1,  condition =  mat2 ,libType = libType )
colnames(Design) = c("condition" , "type") ;

chdata = ceiling(chdata) ;

dds <- DESeqDataSetFromMatrix(countData = chdata,  colData = Design, design = ~ condition)

colData(dds)$condition <- factor( colData(dds)$condition,levels=c("FT","Sham") ) ;

dds$condition <- relevel(dds$condition, ref="Sham") ;
#dds <- estimateSizeFactors(dds, controlGenes=rpos ) ;

dds <- DESeq(dds)
res <- results(dds)


write.table(res , file="DEseq-res.txt", col.names = T, row.names = T ) ;
setwd(chdir) ;
}





