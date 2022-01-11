function makeMatrixtracks(id1 , id2)
 % the line below is not needed 
 %addpath('/home/vibhor/oldCodes/CD4/codes/') ;


dirs = { 'PLC'
'ILC'
'CGC'
'InsularCortex'
'NaCshell'
'NaCcore'
'LS'
'MS'
'BNST'
'MPOA'
'PVN'
'MDT'
'LH'
'VMH'
'BLA'
'CMA'
'dorDG'
'dorCA1'
'HB'
'venDG'
'venCA1'
'Arc'
'VTA'
'PAG'
'Raphe'
'LC'
'Snt'
} ;

cpath = pwd() ;
addpath(cpath) ;

[chrin chrin2 chrin4] = textread( 'rn5/bins-200-property-log' , '%s %d %d') ;

%fid1 = fopen('trackDb.txt' , 'w' ) ;
 data2save = load( 'rn5/bins-200-gc') ;
 gc = data2save(:, 3) ;
 
[pchr ploc1 ploc2] = textread('allpeaksPval5FC4.bed' , '%s %d %d') ; 
peaks =   peaks2posrn5( pchr, ploc1, ploc2, ones(size(ploc1))  ) ; 
size(peaks)
for kl = id1:id2
currdir = pwd() ;
cd(dirs{kl}) ;

vstepnames = textread( 'chipnames' , '%s') ;
cvstepnames = textread( 'inputnames' , '%s') ;

[bedro bedco] = size(vstepnames) ;


peakro = numel(peaks(:, 1)) ; 
cpeak = ceil(mean(peaks(:, 1:2), 2)) ;
peakscore = zeros(peakro , bedro) ;
peakscore1 = zeros(peakro , bedro) ;
meanwin =  ones([5, 1]) ;

bssum =  zeros(1 ,bedro ) ; 

for bi = 1:bedro
 nfilename{bi} = [ vstepnames{bi} '.vstep' ] ;
 datamod(:, bi) = load(nfilename{bi}) ;
 
 cnfilename{bi} = [ cvstepnames{bi} '.vstep' ] ;
 control = load(cnfilename{bi}) ;
 
  data = control ; calculateF = 0 ; GCcalc ;
  ratioGC = mhist ./ (mhist1 + 0.0001) ; calculateF = 1 ; GCcalc ;
  datamod(:, bi) = factorGC .* datamod(:, bi) ;
  gcbias(:, bi) = mhist ;
  

sdatamod(:, bi) = convn(datamod(:, bi) , meanwin , 'same') ;

for j = 1:peakro
peakscore(j,bi) =  sum( datamod(cpeak(j)-5: cpeak(j)+5, bi)) ;
peakscore1(j,bi) =  max( sdatamod(cpeak(j)-5: cpeak(j)+5, bi)) ;
end

[B sidx] = sort( peakscore(:, bi) , 'descend') ;
size(sidx)
for j = 1:10000
bssum(bi) = bssum(bi) + mean(datamod(cpeak(sidx(j))-1:cpeak(sidx(j)) +1 ,bi) ) ; 
end
bssum(bi) = bssum(bi) / (10000* (mean(datamod(:,bi)) + 0.01)) ;

amean = mean(peakscore(:,bi)) ;

sdatamod(:,bi) = sdatamod(:,bi) / amean ;

end
% colr = floor(kl * (255 /27)) ;
% colg = floor( (27 - kl) * (255 /27)) ;
%colb = floor (kl * (255 /27)) ;       
% agcbias{kl} = gcbias ;

%ploc2s = ceil( [ data2save(cpeak,1:2)  peakscore1 ] );

dlmwrite('signal2noise' , bssum , 'delimiter' , '\t') ;
dlmwrite('peakscore' , peakscore , 'delimiter' , '\t' , 'precision' , '%.2f' ) ;
dlmwrite( 'peakscore1' , peakscore1 , 'delimiter' , '\t' , 'precision' , '%d' )

%  for bi = 1:bedro
%  
%  
%  nfilename{bi} = [ 'GCcorr1-track-' vstepnames{bi} '.vstep' ] ;
% % 
% % fprintf(fid1, [ '\ntrack ' vstepnames{bi} '\n'] ) ;
% % fprintf(fid1, [ 'bigDataUrl  data/' nfilename{bi} '.bw\n'] ) ;
% % fprintf( fid1 , ['shortLabel ' vstepnames{bi} '\n'] ) ;
% % fprintf( fid1 , ['longLabel ' vstepnames{bi}  '\n'] ) ;
% % fprintf( fid1, 'type bigWig \n' ) ;
% % fprintf( fid1, ['color ' num2str(colr) ' ' num2str(colg) ' '  num2str(colb)  '\n' ]) ;
% % fprintf( fid1, 'maxHeightPixels 60:32:2 \n' ) ;
% % fprintf( fid1, 'autoScale on\n' ) ;
% % fprintf( fid1, 'alwaysZero on \n \n \n \n' ) ;
% % 
%  
%  fid = fopen(nfilename{bi} , 'w' ) ;
%  fprintf(fid, 'track type=wiggle_0 name="%s"description="%s"\n', vstepnames{bi}, vstepnames{bi}) ;   
%  
%  chr = chrin{data2save(1,1)} ;
%  fprintf(fid, 'variableStep chrom=%s span=199 \n',  chr ) ;
%  fprintf(fid, '%d %d\n' , data2save(1, 2) , sdatamod(1, bi)) ;
%  
%  for i = 2:numel(datamod(:, 1)) 
%      
%  if( data2save(i, 1) ~=  data2save((i-1), 1 ) )
%  chr = chrin{data2save(i ,1) } ; 
%  fprintf(fid, 'variableStep chrom=%s span=199 \n',  chr ) ;
%  end
%  fprintf(fid, '%d %d\n' , data2save(i, 2) , sdatamod(i, bi) ) ;   
%  end
%  
%  fclose(fid) ;
% % 
% command = ['perl ../convert.pl  ' nfilename{bi}  ] ;
% system(command) ;
% % 
% command = [' wigToBigWig ' nfilename{bi} '.wig   ../bins-200-property-log ' nfilename{bi}  '1.bw' ] ;
%  system(command) ;
%  end
% 
 cd (currdir) ;
 end

%fclose(fid1) ;




