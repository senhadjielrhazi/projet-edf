function [D,q,L]= kstwo(x, y, flag,verb_flag) 
// flag = ["two.sided", "less", "greater"] 
  [lhs,rhs]=argn(0) ;
  if rhs <= 3 then verb_flag = %f ;end 
  if rhs <= 2 then flag = "two sided"; end 
  nx = size(x,'*'); 
  ny = size(y,'*'); 
  // rajouter des typeof ? 

  // we now compute S = S_Nx(u) - S_Ny(u) 
  // which is a piecewize function jumping at x or y elements 
  [Z,zk]=gsort( [x(:);y(:)]','g','i');
  // first compute a density 
  z= - 1/ny*ones(Z) ; 
  z( zk <= nx ) = 1/nx ; 

  // take care of redondant values in x and y 
  [Zu,Zui]=unique(Z) 
  if size(Zu,'*')<>size(Z,'*') then 
    // we have redundant values 
    i=0*Z;
    i(Zui)=1;
    ii=find(i==0); 
    // ii gives indices of redundant values not selected in Zu 
    // we move the associated values in the Zu selection 
    for vi=ii, 
      j=find( Zu == Z(vi)) 
      z(Zui(j))= z(Zui(j))+z(vi) 
    end 
    // we delete the points which have been moved 
    z(ii)=[];
  end 

  // now compute the values of jump of S 
  S = cumsum(z); 
  METHOD = "Two-sample Kolmogorov-Smirnov test"; 
  n =  nx * ny / (nx + ny)

  select flag 
   case "two sided" then 
    D = maxi(abs(S)) ; 
    T = "D" ; 
    q = cdfks(sqrt(n)*D)
   case "less" then 
    D= -min(S) ; 
    T = "D^-"
    q = exp(-2 *n* D^2)
   case "greater" then 
    D = max(S) ; 
    T = "D^+"
    q = exp(-2 *n* D^2)
  else error('flag should be ""two sided"" or ""less"" or ""greater"" '); 
     return 
  end 
  L= list(METHOD,T,q,D) 
  if verb_flag then 
    mprintf(METHOD+"\n") 
    mprintf(T+" = %5.2f, p-value = " + string(q) +"\n",D);
    mprintf("alternative hypothesis: "+ flag + "\n");
  end
endfunction 


function K_PQ = cdfks(x) 
stacksize(10000000);
I = 1:1000000;
K_PQ = 1+2*sum((-1)^I.*exp(-2*I^2*x^2));
endfunction