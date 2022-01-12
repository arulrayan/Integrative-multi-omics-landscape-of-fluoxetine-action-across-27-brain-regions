#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sam.h"



 void readBAMFile(char* filename, char *outfilename, char *limitsf, int PE ) ;

void readBedFile(char* filename, char *outfilename, char *limitsfi) ;
void readBedGraph(char* filename, char *outfilename, char *limitsf) ;

main(int argc, char **argv  )
{
int PE ;
char command[300] , fname[100];
if(argc < 4) { printf("\nUsage : arrange2binsUCSC input.bed output.vstep binlimitsFile bed/bam/sam/bedGraph paired_end \n" ) ; exit(1); } ;

if(argc > 5) {  PE = atoi(argv[5]) ;}

if(argc > 4)
{
if(strcmp(argv[4], "bed") == 0) { readBedFile( argv[1], argv[2], argv[3]) ; }

if((strcmp(argv[4], "bam") == 0) || (strcmp(argv[4], "sam") == 0)) { 
 readBAMFile( argv[1], argv[2], argv[3], PE ) ;}
 }

if(strcmp(argv[4], "bedGraph") == 0) { readBedGraph( argv[1], argv[2], argv[3]) ; }
}

void readBedFile(char* filename, char *outfilename, char *limitsf)
{
FILE *infile, *ofile, *ofile1, *logfile ;
int i, j, line ,count, chrloc, chrloc1, idx, idx1, totalbins, lineidx[40], chrinfo[40][2] , chrnum;
long maxloc , minloc;
char instr[200], tempstr[200], prevchr[200], currchr[200] , sign[4], chromosomes[40][20] ;
char trackname[300] ;
int *locarray , *chrarray, tagcounts, wsize, taglength , idxj, chrcount, columns ;
float *data, mean, ratio ;


maxloc = 0 ;
minloc = 10000000000 ;
wsize = 200 ;
taglength = 200 ;


logfile = fopen( limitsf, "r") ;

chrcount = 1 ;
fscanf(logfile, "%s", instr) ; sprintf(chromosomes[chrcount],"%s" ,instr) ;
while(!feof(logfile))
{
fscanf(logfile, "%s", instr) ; chrinfo[chrcount][0] = atoi(instr) ; 
fscanf(logfile, "%s", instr) ; chrinfo[chrcount][1] = atoi(instr) ;
chrcount = chrcount + 1;
fscanf(logfile, "%s", instr) ; sprintf(chromosomes[chrcount],"%s" ,instr) ;
}
fclose(logfile) ;

totalbins = 0 ;
lineidx[0] = 0 ;

for(i = 1 ; i < chrcount ; i++)
{
lineidx[i] = lineidx[i-1] + (chrinfo[i][1] - chrinfo[i][0]) / wsize + 1;
printf("%d  ", lineidx[i]) ;
totalbins = lineidx[i] ;
}


line = 0 ;
count = 0 ;



//printf("loaded file") ; scanf("%s" , tempstr) ;

data = (float*)malloc(totalbins * sizeof(float)) ;
locarray = (int*)malloc(totalbins * sizeof(int)) ;

//printf("allocated memory %d " ,totalbins) ; scanf("%s" , tempstr) ;

for(i = 1 ; i < chrcount ; i++)
{
printf("%s %d %d %d %d\n", chromosomes[i], lineidx[i-1], chrinfo[i][0], lineidx[i] -1, chrinfo[i][1]) ;
for(j = lineidx[i-1] ; j < lineidx[i]; j++)
{
locarray[j] = chrinfo[i][0] + (j - lineidx[i-1]) * wsize ;  
data[j] = 0 ;
}
}

//printf("going to load file") ; scanf("%s" , tempstr) ;

infile = fopen(filename , "r") ;
sprintf( prevchr, "chr0") ; 

while(!feof(infile))
{
line = line + 1 ;
fscanf(infile, "%s", currchr) ;
fscanf(infile, "%s", instr) ; chrloc = atoi(instr) ;
fscanf(infile, "%s", instr) ; chrloc1 = atoi(instr) ; 
fscanf(infile, "%s", instr) ; 
fscanf(infile, "%s", instr) ; 
fscanf(infile, "%s" , sign) ;


/*
printf("%s %d %d %d\n", currchr, line, minloc, maxloc ) ;
*/

if(strcmp(currchr , prevchr) != 0)
{
// printf("%s %d %d %d\n", prevchr, line, minloc, maxloc ) ;
count = 50 ;
for(j = 1 ; j < chrcount ; j++)
{
if(strcmp(chromosomes[j], currchr) == 0) { count = j ; }
}
strcpy(prevchr,  currchr) ;
}

if(count < chrcount)
{

if(strcmp(sign, "+") == 0)
{
chrloc1 = chrloc + taglength ;
}
else
{
chrloc = chrloc1 - taglength ;  
}

idx = lineidx[count-1] + (chrloc - chrinfo[count][0])/ wsize; 
idx1 = lineidx[count-1] + (chrloc1 - chrinfo[count][0])/ wsize; 

if((chrloc > chrinfo[count][0]) && (chrloc < chrinfo[count][1]) && (idx >= lineidx[count-1]) &&( idx < lineidx[count])) 
{
ratio = (locarray[idx] + wsize - chrloc)/ (taglength + 0.1) ;
data[idx] = data[idx] + ratio ;
}

if((chrloc1 > chrinfo[count][0]) && (chrloc1 < chrinfo[count][1])&& (idx1 >= lineidx[count-1]) &&( idx1 < lineidx[count])) 
{

if(idx != idx1) {
ratio = (chrloc1 - locarray[idx1])/ (taglength + 0.1) ;  
data[idx1] = data[idx1] + ratio ;
}

if((idx1 - idx) > 1)
{
for(idxj = idx + 1 ; idxj < idx1 ; idxj++)
{ data[idxj] = data[idxj] + wsize / (taglength + 0.1) ; }
}

}

}

}
fclose(infile) ;

/*
mean = 0 ;
for(j = 0 ; j < totalbins ;j++)
{
mean = mean + data[j] / totalbins ; 
}
*/


ofile = fopen( outfilename, "w") ;
//sprintf(trackname, "track-%s", argv[2]) ;
//ofile1 = fopen(trackname, "w") ;

//fprintf( ofile1, "track type=wiggle_0 name=\"%s\" description=\"%s\"\n", trackname, trackname ) ;

for(i = 1 ; i < chrcount ; i++)
{

//fprintf( ofile1, "variableStep chrom=%s span=200\n", chromosomes[i] ) ;
for(j = lineidx[i-1] ; j < lineidx[i]; j++)
{
tagcounts = (int)(data[j] + 0.3) ;
fprintf( ofile, "%d %d %d\n", chrnum, locarray[j], tagcounts ) ;
//fprintf( ofile, "%d\n", tagcounts ) ;
//fprintf( ofile1, "%d %d\n", locarray[j], tagcounts ) ;
}
}

fclose(ofile) ;
//fclose(ofile1) ;
}



void readBAMFile(char* filename, char *outfilename, char *limitsf, int PE )
{
    
    FILE *infile, *ofile, *ofile1, *logfile ;
int i, j, line ,count, chrloc, chrloc1, idx, idx1, totalbins, lineidx[40], chrinfo[40][2] , chrnum;
long maxloc , minloc;
char instr[200], tempstr[200], prevchr[100], currchr[100] , sign[4], chromosomes[40][20];
char trackname[300] ;
int *locarray , *chrarray, tagcounts, wsize, taglength , idxj, chrcount, columns ;
float *data, mean, ratio ;
int *idmap,  signInt,  signInt1;


maxloc = 0 ;
minloc = 10000000000 ;
wsize = 200 ;
taglength = 200 ;


logfile = fopen( limitsf, "r") ;

chrcount = 1 ;
fscanf(logfile, "%s", instr) ; sprintf(chromosomes[chrcount],"%s" ,instr) ;
while(!feof(logfile))
{
fscanf(logfile, "%s", instr) ; chrinfo[chrcount][0] = atoi(instr) ; 
fscanf(logfile, "%s", instr) ; chrinfo[chrcount][1] = atoi(instr) ;
chrcount = chrcount + 1;
fscanf(logfile, "%s", instr) ; sprintf(chromosomes[chrcount],"%s" ,instr) ;
}
fclose(logfile) ;

totalbins = 0 ;
lineidx[0] = 0 ;

for(i = 1 ; i < chrcount ; i++)
{
lineidx[i] = lineidx[i-1] + (chrinfo[i][1] - chrinfo[i][0]) / wsize + 1;
printf("%d  ", lineidx[i]) ;
totalbins = lineidx[i] ;
}


line = 0 ;
count = 0 ;


//printf("loaded file") ; scanf("%s" , tempstr) ;

data = (float*)malloc(totalbins * sizeof(float)) ;
locarray = (int*)malloc(totalbins * sizeof(int)) ;

//printf("allocated memory %d " ,totalbins) ; scanf("%s" , tempstr) ;

for(i = 1 ; i < chrcount ; i++)
{
printf("%s %d %d %d %d\n", chromosomes[i], lineidx[i-1], chrinfo[i][0], lineidx[i] -1, chrinfo[i][1]) ;
for(j = lineidx[i-1] ; j < lineidx[i]; j++)
{
locarray[j] = chrinfo[i][0] + (j - lineidx[i-1]) * wsize ;  
data[j] = 0 ;
}
}
    
        
  samfile_t *fp;
  if ((fp = samopen(filename, "rb", 0)) == 0) 
  {
 printf("\n could not open file %s ", filename );
    exit(1);
  }

 
 idmap = (int*) malloc(fp->header->n_targets * sizeof(int)) ;
  printf("\n  reading %s \n" ,filename );
 
  for ( i = 0; i <  fp->header->n_targets ;  i++)
  {
  char *currchr = fp->header->target_name[i]; 
idmap[i] = 0 ;
  for (j = 0 ; j < chrcount ; j++) { 
 
if( strcmp ( chromosomes[j], currchr)  == 0 )
  {
  idmap[i] = j ;
  break ;
  } 
  }
}

 // printf("\n  read header %s \n" ,filename ); scanf("%s", tempstr) ;

if (PE == 1)
{
bam1_t *b = bam_init1();
bam1_t *b1 = bam_init1();

while (samread(fp, b) >= 0)
 {
if(BAM_FPROPER_PAIR)
{
samread(fp, b1) ;
const bam1_core_t *c = &b->core;
const bam1_core_t *c1 = &b1->core;

    if ( (c->tid != -1) & (c1->tid != -1))
    {
      if (c->flag & BAM_FREVERSE) { chrloc = bam_calend(c, bam1_cigar(b)) ;  signInt = -1 ;}
            else {    chrloc  =   c->pos;  signInt = 1 ;  }

   
      if (c1->flag & BAM_FREVERSE) { chrloc1 = bam_calend(c1 , bam1_cigar(b1)) ;  signInt1 = -1 ;}
            else {    chrloc1  =   c1->pos;  signInt1 = 1 ;  }

if( ((signInt * signInt1)  == -1) & ( c->tid == c1->tid) )
{
count = 0 ;
    count = idmap[c->tid] ;
if(count > 0)
{
  
idx = lineidx[count-1] + (chrloc - chrinfo[count][0])/ wsize; 
idx1 = lineidx[count-1] + (chrloc1 - chrinfo[count][0])/ wsize; 

if((chrloc > chrinfo[count][0]) && (chrloc < chrinfo[count][1]) && (idx >= lineidx[count-1]) &&( idx < lineidx[count])) 
{
ratio = (locarray[idx] + wsize - chrloc)/ (taglength + 0.1) ;
data[idx] = data[idx] + ratio ;
}

if((chrloc1 > chrinfo[count][0]) && (chrloc1 < chrinfo[count][1])&& (idx1 >= lineidx[count-1]) &&( idx1 < lineidx[count])) 
{

if(idx != idx1) {
ratio = (chrloc1 - locarray[idx1])/ (taglength + 0.1) ;  
data[idx1] = data[idx1] + ratio ;
}

if((idx1 - idx) > 1)
{
for(idxj = idx + 1 ; idxj < idx1 ; idxj++)
{ data[idxj] = data[idxj] + wsize / (taglength + 0.1) ; }
}

}
 }
}} }
  bam_destroy1(b);
    bam_destroy1(b1);
  
}
}

else
{
bam1_t *b = bam_init1();
 
  while (samread(fp, b) >= 0)
  {
    const bam1_core_t *c = &b->core;
    if (c->tid != -1)
    {
      if (c->flag & BAM_FREVERSE) { chrloc = bam_calend(c, bam1_cigar(b)) ;  signInt = -1 ;}
            else {    chrloc  =   c->pos;  signInt = 1 ;  }
  
count = 0 ;
    count = idmap[c->tid] ;

if(count > 0)
{
  
idx = lineidx[count-1] + (chrloc - chrinfo[count][0])/ wsize; 

if((chrloc > chrinfo[count][0]) && (idx >= lineidx[count-1]) &&( idx < lineidx[count])) 
{
data[idx] = data[idx] + 1 ;
}

}
 
    }     
  }  
   bam_destroy1(b);
}
  
  samclose(fp);
  
  
ofile = fopen( outfilename, "w") ;

for(i = 1 ; i < chrcount ; i++)
{

for(j = lineidx[i-1] ; j < lineidx[i]; j++)
{
tagcounts = (int)(data[j] + 0.3) ;
fprintf( ofile, "%d\n", tagcounts ) ;
}
}

fclose(ofile) ;

}



void readBedGraph(char* filename, char *outfilename, char *limitsf)
{
FILE *infile, *ofile, *ofile1, *logfile ;
int i, j, line ,count, chrloc, chrloc1, idx, idx1, totalbins, lineidx[40], chrinfo[40][2] , chrnum;
long maxloc , minloc;
char instr[200], tempstr[200], prevchr[200], currchr[200] , sign[4], chromosomes[40][20] ;
char trackname[300] ;
int *locarray , *chrarray, tagcounts, wsize, taglength , idxj, chrcount, columns ;
float *data, mean, ratio , ltags;


maxloc = 0 ;
minloc = 10000000000 ;
wsize = 200 ;
taglength = 200 ;


logfile = fopen( limitsf, "r") ;

chrcount = 1 ;
fscanf(logfile, "%s", instr) ; sprintf(chromosomes[chrcount],"%s" ,instr) ;
while(!feof(logfile))
{
fscanf(logfile, "%s", instr) ; chrinfo[chrcount][0] = atoi(instr) ; 
fscanf(logfile, "%s", instr) ; chrinfo[chrcount][1] = atoi(instr) ;
chrcount = chrcount + 1;
fscanf(logfile, "%s", instr) ; sprintf(chromosomes[chrcount],"%s" ,instr) ;
}
fclose(logfile) ;

totalbins = 0 ;
lineidx[0] = 0 ;

for(i = 1 ; i < chrcount ; i++)
{
lineidx[i] = lineidx[i-1] + (chrinfo[i][1] - chrinfo[i][0]) / wsize + 1;
printf("%d  ", lineidx[i]) ;
totalbins = lineidx[i] ;
}


line = 0 ;
count = 0 ;



//printf("loaded file") ; scanf("%s" , tempstr) ;

data = (float*)malloc(totalbins * sizeof(float)) ;
locarray = (int*)malloc(totalbins * sizeof(int)) ;

//printf("allocated memory %d " ,totalbins) ; scanf("%s" , tempstr) ;

for(i = 1 ; i < chrcount ; i++)
{
printf("%s %d %d %d %d\n", chromosomes[i], lineidx[i-1], chrinfo[i][0], lineidx[i] -1, chrinfo[i][1]) ;
for(j = lineidx[i-1] ; j < lineidx[i]; j++)
{
locarray[j] = chrinfo[i][0] + (j - lineidx[i-1]) * wsize ;  
data[j] = 0 ;
}
}

//printf("going to load file") ; scanf("%s" , tempstr) ;

infile = fopen(filename , "r") ;
sprintf( prevchr, "chr0") ; 

while(!feof(infile))
{
line = line + 1 ;
fscanf(infile, "%s", currchr) ;
fscanf(infile, "%s", instr) ; chrloc = atoi(instr) ;
fscanf(infile, "%s", instr) ; chrloc1 = atoi(instr) ; 
fscanf(infile, "%s", instr) ; ltags = atof(instr) ;


/*
printf("%s %d %d %d\n", currchr, line, minloc, maxloc ) ;
*/

if(strcmp(currchr , prevchr) != 0)
{
// printf("%s %d %d %d\n", prevchr, line, minloc, maxloc ) ;
count = 50 ;
for(j = 1 ; j < chrcount ; j++)
{
if(strcmp(chromosomes[j], currchr) == 0) { count = j ; }
}
strcpy(prevchr,  currchr) ;
}

if(count < chrcount)
{

if(strcmp(sign, "+") == 0)
{
chrloc1 = chrloc + taglength ;
}
else
{
chrloc = chrloc1 - taglength ;  
}

idx = lineidx[count-1] + (chrloc - chrinfo[count][0])/ wsize; 
idx1 = lineidx[count-1] + (chrloc1 - chrinfo[count][0])/ wsize; 

if((chrloc > chrinfo[count][0]) && (chrloc < chrinfo[count][1]) && (idx >= lineidx[count-1]) &&( idx < lineidx[count])) 
{
ratio = (locarray[idx] + wsize - chrloc)/ (taglength + 0.1) ;
data[idx] = data[idx] + ratio * ltags ;
}

if((chrloc1 > chrinfo[count][0]) && (chrloc1 < chrinfo[count][1])&& (idx1 >= lineidx[count-1]) &&( idx1 < lineidx[count])) 
{

if(idx != idx1) {
ratio = (chrloc1 - locarray[idx1])/ (taglength + 0.1) ;  
data[idx1] = data[idx1] + ratio *ltags;
}

if((idx1 - idx) > 1)
{
for(idxj = idx + 1 ; idxj < idx1 ; idxj++)
{ data[idxj] = data[idxj] + ltags / (idx1 - idx) ; }
}

}

}

}
fclose(infile) ;

/*
mean = 0 ;
for(j = 0 ; j < totalbins ;j++)
{
mean = mean + data[j] / totalbins ; 
}
*/


ofile = fopen( outfilename, "w") ;
//sprintf(trackname, "track-%s", argv[2]) ;
//ofile1 = fopen(trackname, "w") ;

//fprintf( ofile1, "track type=wiggle_0 name=\"%s\" description=\"%s\"\n", trackname, trackname ) ;

for(i = 1 ; i < chrcount ; i++)
{

//fprintf( ofile1, "variableStep chrom=%s span=200\n", chromosomes[i] ) ;
for(j = lineidx[i-1] ; j < lineidx[i]; j++)
{
//tagcounts = (int)(data[j] + 0.3) ;
//fprintf( ofile, "%d %d %.1f\n", chrnum, locarray[j], data[j] ) ;
fprintf( ofile, "%f\n",  data[j] ) ;
//fprintf( ofile1, "%d %d\n", locarray[j], tagcounts ) ;
}
}

fclose(ofile) ;
//fclose(ofile1) ;
}








