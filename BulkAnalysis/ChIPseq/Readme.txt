The folders contain read-counts and other information for every brain Region
Notice here about correct version of DEseq2 to replicate the exact result. 
 
Results may slightly vary based on version of DEseq2 however overall it will remain similar to some extent. 


Every folder corresponding to name of a Brain region has following files
1. CorrectedScore  : it contains read-counts after GC correction for all peaks in the union list for all regions  

3. sampletype: contain information about if they ar sham or fluoxetine treated   

4. chipnames : contains information about the replicate number corresponding to every column in read-count matrix (notice Replicate 3 was dropped because of quality issue).   


DEseq-res.txt is the output of the code runDEseq.R for peaks from unionlist in corresponding region folders


$Note$ Bigger files are zipped, so they must be unzipped to be readable by code. 

*********In order to caculate the CorrectedScore (read-count data) on union peak list  (in file allpeaksPval5FC4)for every brain region *********************

Their original fastq files can downloaded from GEO database and aligned to rn5 genome to get a bam file. 
The bam files are converted to vstep files which contain read-counts in 200 bp bins genome wide.

The vstep files are not here as they are very big.
but they can be made from bam files using the code  arrange2binsUCSC(need compilation)
 and command

 arrange2binsUCSC input.bam output.vstep rn5/bins-200-property-log bam 0

the arrange2binsUCSC.c program can be compiled using following command

cc -o arrange2binsUCSC arrange2binsUCSC.c -I PATH/samtools-0.1.17  PATH/samtools-0.1.17/libbam.a  -lz -lpthread

where PATH is the path of folder where samtools-0.1.17 is compiled in the user compter/server.

the vstepfiles should be put in their corresponding folder, the information of which is provided in file called as "chipnames"  and "inputnames" in folder for every brain region. 


after puting the correponding vstepfiles in respective folders, the matlab code makeMatrixtracks.m can be used. 

makeMatrixTracks.m does GC-bias correction of binned read-counts in 200 bins in vstep files using the input library.

Notice the code makeMatrixTracks.m was also used to save the GC-bias-corrected tracks of binned read-counts (in 200 bins). 

Following trackhub of the ChIP-seq tracks can be uploaded to UCSC genome browser
https://s3-ap-southeast-1.amazonaws.com/ratbraindata.store.genome.sg/ChipWig/hub.txt
or 
might be accessible as a session in UCSC genome browser
https://genome-asia.ucsc.edu/s/Nirmala/rn5_ChIPseq_FT



**************************  


                                                                                                                   