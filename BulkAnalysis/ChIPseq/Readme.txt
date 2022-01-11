The folders contain read-counts and other information for every brain Region
Notice here older version of DEseq2 was used. 
Results may slightly vary based on version of DEseq2 however overall it will remain same to some extent. 


Every folder has following files
1. CorrectedScore  : it contains read-counts after GC correction for all peaks in the union list for all regions (total 48006 peaks) 

3. sampletype: contain information about if they ar sham or fluoxetine treated   

4. Chipnames : contains information about the replicate number corresponding to every column in read-count matrix (notice Replicate 3 was dropped because of bad libraries).   

DEseq-res.txt is the output of the code runDEseq.R for peaks from unionlist for all region 

DEseq-res2r-2rep.txt is the output of the code runDEseq.R for peaks chosen for ever region


Notice here that Rep3 was dropped during quality control. Only Rep1 and Rep2 were used. Their original fastq files can downloaded from GEO database. 

Bigger files are zipped, so they must be unzipped to be readable by code. 
The Vstep files are not here as they are very big.
but they can be made from bam files using the code  arrange2binsUCSC(need compilation)
 and command

 arrange2binsUCSC input.bed/bam/bedGraph output.vstep binlimitsFile bed/bam/sam/bedGraph paired_end










 