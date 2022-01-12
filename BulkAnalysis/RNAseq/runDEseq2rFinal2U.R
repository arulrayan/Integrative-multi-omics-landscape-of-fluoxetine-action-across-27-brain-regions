library(DESeq2)


dirs = c("PLC" , "ILC" , "CGC" , "InsularCortex" , "NaCshell" , "NaCcore" , "LS" , "MS" , "BNST" , "MPOA" , "PVN" , "MDT" , "LH" , "VMH" , "BLA" , "CMA" , "dorDG" , "dorCA1" , "HB" , "venDG" , "venCA1" , "Arc" , "VTA" , "PAG" , "Raphe" , "LC" , "SNT");


gnames = read.table( "Rn5-ensemble2GeneName.txt") ;
gnames = as.matrix(gnames) ;
rownames(gnames) = gnames[, 1] ;

sampleDataRNA = read.delim('LibraryList_RNA_Seq.txt')

pos = which(sampleDataRNA$Replicate != 'Rep2') ;
sampleDataRNA  = sampleDataRNA[pos,] ; 

fro = nrow( sampleDataRNA) ;
files = as.matrix(sampleDataRNA$LibID) ;


expr = matrix(0, 29188, fro) ;

for (i in 1:fro)
{
texpr = read.table(files[i] , header = FALSE ) ;
expr[, i] = as.matrix(texpr[1:29188, 2]) ;
}

colnames(expr) = files ;
rownames(expr) = texpr[1:29188,1] ;
NMid = as.matrix(texpr[1:29188,1]) ;

proms = read.table("Rn5-ensembleGenes-promoters.txt", header=FALSE) ;
NMp = as.matrix(proms[, 1]) ;
rownames(proms) =  as.matrix(proms[, 1]) ;


pos = as.matrix(which( NMid %in% NMp )) ;
expr = expr[NMid[pos] , ] ;
nprom = proms[NMid[pos],2:5 ] ;


expr =  expr[ , order(colnames(expr) ) ] ;
oexpr = expr ;

sampleDataRNA=sampleDataRNA[order((sampleDataRNA$LibID)),] ;

rt = read.table('remRep2Alt-2020.txt', fill = T) ;

dispfile = "inputoutput/dispersion.pdf" ;
pdf(dispfile) ;

allFDR = matrix(0, 29093, 27) ; 
allFC = matrix( 0, 29093, 27) ;

for (tis in 1:27 ) 
{
rt1 = as.matrix(rt[tis,]) ;
pos = which( sampleDataRNA$Brain.region == dirs[tis] );
rpos = which(rt1 == TRUE) ;
pos = pos[rpos] ;

tbams = as.matrix(sampleDataRNA$LibID[pos]) ;
chdata  =  expr[ ,tbams] ;
chdatao = expr[, tbams] ;

chro = dim(chdata)[1] ;
chco = dim(chdata)[2] ;
rt2 = rt1[1:chco] ;

maxch = apply(chdatao, 1, 'max') ;
maxch60 = quantile(maxch, probs = 0.60) ;
rpos = which(maxch > maxch60) ;

nprom1 = nprom ;

batchs = sampleDataRNA$Replicate[pos] ;
rbatch = unique(batchs) ;
scondition = sampleDataRNA$Treatment[pos] ;


mat1 = as.matrix(rownames(expr)) ;
mat2 = as.matrix(scondition) ;


libType<-rep("paired-end", dim(chdata)[2] ) ;
Design = data.frame(  row.names = ,  condition =  mat2 , batch= batchs ,libType = libType ) ;
colnames(Design) = c("condition" , "batch" , "type") ;
chdata = ceiling(chdata) ;
dds <- DESeqDataSetFromMatrix(countData = chdata,  colData = Design, design = ~ batch + condition ) ;
dds$condition <- factor(dds$condition , levels=c("Fluoxetine", "Sham"))
dds$batch <- factor( dds$batch , levels= rbatch) ;
dds$condition <- relevel(dds$condition, ref="Sham");
 dds <- DESeq(dds,  full=design(dds), reduced = ~ batch) ;

plotDispEsts( dds);


res <- results(dds, contrast=c("condition" , "Fluoxetine" , "Sham")) ;

roname = as.matrix(rownames(res)) ;
roname1 = roname ;
pos = which( roname %in% gnames[,1] ) ;
roname1[pos] = gnames[roname[pos],2] ;
chdata = chdata[ , order(colnames(chdata) ) ] ;

d2save = data.frame(roname1, nprom1 ,res ) ;

filestring = paste(  "inputoutput/" , dirs[tis] , ".txt",  sep = "") ;
write.table( d2save , file=filestring , col.names = F , row.names = T) ;

allFDR[, tis ] = res[,6] ; 
allFC[, tis] = res[,2] ;

}

write.table(allFDR , file="inputoutput/allFDR.txt" ) ;
write.table(allFC , file="inputoutput/allFDR.txt") ;


dev.off() ;
