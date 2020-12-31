
Imputation <- function(XX){
  w = which(is.na(XX), arr.ind=TRUE)
  r=w[,1]
  c=w[,2]
  for (i in 1:length(r)){
    XX[r[i],c[i]]=0
  }
  XX=matrix(as.numeric(XX),ncol=1)
  if (length(r)!=0){
    for (i in 1:length(r)){
      if (r[i]==1){
        if (XX[r[i]+1,1]==0){
          if (XX[r[i]+2,1]==0){
            if (XX[r[i]+3,1]==0){
              temp=XX[r[i]+4,1]
            }else{
              temp=XX[r[i]+3,1]
            }
          }else{
            temp=XX[r[i]+2,1]
          }
        }else{
          temp=XX[r[i]+1,1]
        }
      }else if(r[i]==nrow(XX)){
        if (XX[r[i]-1,1]==0){
          if (XX[r[i]-2,1]==0){
            if (XX[r[i]-3,1]==0){
              temp=XX[r[i]-4,1]
            }else{
              temp=XX[r[i]-3,1]
            }
          }else{
            temp=XX[r[i]-2,1]
          }
        }else{
          temp=XX[r[i]-1,1]
        }
      }else if(XX[r[i]+1,1]==0){
        temp=XX[r[i]-1,1]
      }else{
        temp=0.5*(XX[r[i]-1,1]+XX[r[i]+1,1])
      }
      XX[r[i],1]=temp
    }
  }
  return(XX)
}
