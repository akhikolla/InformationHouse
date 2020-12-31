zeros <-function (n,m=n,o=0) {
# zeros(n)     returns an n-by-n matrix of 0s. 
# zeros(n,1)   returns a vector of 0s 
# zeros(n,m)   returns an n-by-m matrix of zeros
# zeros(n,m,o) returns an 3D matrix  of zeros

# ALU

if (m==1) { # vector wird zurueckgegeben
   return(c(1:n)*0) ;   
}else{      # return n-by-m matrix of ones.         
  if  (o==0){
     return(matrix(0,n,m));
   }else{   #3D matrix
     nullen = rep(0, m*n*o);  # soviele nullen wie in die 3D matrix pasen
     return(array(nullen,c(n,m,o)));
    
   } # end  if  (o==0)
  } # end if (m==1) 
} # end  function  zeros
