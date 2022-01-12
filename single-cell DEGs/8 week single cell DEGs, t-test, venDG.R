log2addsmall <- function(x){ return(log2(x+0.01)) }
addsmall <- function(x){ return((x+0.01)) }

datafile <- "/../Single-cell collated tp10k_venDG.txt"

data <- read.table(datafile, stringsAsFactors=F)

underscores <- unlist(gregexpr("_", names(data)))
underscore1 <- underscores[seq(from=1, to=length(underscores), by=2)]
celltype <- substring(names(data), 1, underscore1-1)

threshold <- 0.5
for(cell in unique(celltype)){
	reps <- lapply(names(data), function(ch) grep(cell, ch))
	dataF <- data[, which(reps==1)]
	avgs <- unlist(apply(dataF, 1, mean))
	dataF2 <- dataF[avgs>=threshold, ]
	genes <- rownames(dataF2)
	
	Sham <- lapply(names(dataF2), function(ch) grep('Sham', ch))
	ShamF <- dataF2[, which(Sham==1)]
	FTF <- dataF2[, which(lapply(Sham, length)==0)]
	
	### transform to log2 with pseudocount
	ShamF2 <- as.data.frame(apply(ShamF, 1:2, log2addsmall))
	FTF2 <- as.data.frame(apply(FTF, 1:2, log2addsmall))
		
	masterttest <- NULL
	for(i in 1:nrow(ShamF2)){
		ttest <- t.test(ShamF2[i,], FTF2[i,])
		if(i==1){
			masterttest <- unlist(ttest)
		} else {
			masterttest <- rbind(masterttest, unlist(ttest))
		}
	}
	masterttest <- as.data.frame(masterttest, row.names=genes)
	for(j in 1:9){ ### 1st nine columns of ttest output need to be converted to numeric
		masterttest[[j]] <- as.numeric(levels(masterttest[[j]]))[masterttest[[j]]]
	}
	
	### calculate log2FC based on non-log2 values (but with a pseudocount)
	FTFaddsmall <- apply(FTF, 1:2, addsmall)
	ShamFaddsmall <- apply(ShamF, 1:2, addsmall)
	log2FC <- log2(apply(FTFaddsmall, 1, mean)/apply(ShamFaddsmall, 1, mean))
	adjpval <- p.adjust(masterttest[[3]], 'fdr')
	masterttest <- cbind(masterttest, log2FC, adjpval)
	
	filename <- paste0(cell, '_tp10k>=0.5_ttest_stats_fdr.csv')
	write.csv(masterttest, filename)
}	
