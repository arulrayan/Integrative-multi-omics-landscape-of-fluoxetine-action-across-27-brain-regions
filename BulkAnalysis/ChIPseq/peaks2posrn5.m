function  posPeaks = peaks2posrn5(pchr , ploc1, ploc2, ptags )
codepath = '/home/vibhor/vibhor/programmes/DFilter2013' ;
binlimFile = [codepath '/pipeline/rn5/bins-200-property-log'] ;
[chrin chrin2 chrin4] = textread( binlimFile , '%s %d %d') ;
[totalchr cols] = size(chrin2) ;
wsize = 200 ;


chrin1(1) = 0 ;
chrin3(1) =  floor((chrin4(1) - chrin2(1))/wsize) +1 ;
for i = 2:totalchr
chrin1(i) = chrin3(i-1)+1 ;
chrin3(i) =  chrin3(i-1) + floor((chrin4(i) - chrin2(i))/wsize)+ 1  ;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% remove negative locations%%%%%%%%%%%%%%%%%%%%%
 [ro col]  =  size(pchr) ;
count = 0 ;
peakPos = zeros(ro, 2) ;

for i = 1:ro
    chrname = pchr{i} ;
    chr = 0 ;
    for j = 1:totalchr
        
        if(strcmp(chrin(j, 1) , pchr{i} ) == 1)
       chr = j ;
       break ;
        end
    end
    
    if(chr > 0)
    if(((ploc1(i) - 1000) > chrin2(chr) ) & ((ploc2(i) + 1000) < chrin4(chr) ))
        
       chrloc1 = chrin1(chr) + ceil((ploc1(i) - chrin2(chr)) / wsize)  ;
       chrloc2 = chrin1(chr) + ceil((ploc2(i) - chrin2(chr)) / wsize)  ;
       
    count = count + 1 ; 
    peakPos(count, 1) = chrloc1 ;  peakPos(count, 2) = chrloc2 ;
   
    peakPos(count, 3) = ptags(i) ;
   
    end
    
    end
end

if(count < ro)
    peakPos(count + 1:ro,:) = [] ;
end
    
posPeaks = peakPos ;



 
