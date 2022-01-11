
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mi = 0 ; ma = 1 ;
numbins = 100 ;
gap = (ma - mi) / numbins ;
loch = [] ;
mhist = [] ;
mhist1 = [] ;
count = 1 ;
datap = [] ;
 for m = 1: numbins
     pos = find(( gc >= (mi + (m-1) * gap) ) & ( gc <  (mi + m *gap + 0.01) ) ) ; 
    [rows1, cols1] = size(pos) ;
     mhist(count) =  rows1 * cols1 ;
     mhist1(count) = sum(data(pos)) ;
     loch(count) = mi + (m -1) * gap ;
  
  if (calculateF == 1)
       factorGC(pos,1) = ratioGC(count) ;
      if(ratioGC(count) < 0.2)
      factorGC(pos,1) = 0.2 ;
      else
       
        if(ratioGC(count) > 5)
        factorGC(pos,1) = 5 ;
        end
      end

  end
  count = count + 1 ;
  
 end
 
 mhist = mhist / sum(mhist) ;
 mhist1 = mhist1 / sum(mhist1) ;


