umatrixForEsom <- function(Weights, Lines=NULL, Columns=NULL, Toroid=T){
# Umatrix=umatrixForEsom(wts)
# Calculate the Umatrix for given EsomNeurons projection
# INPUT
# EsomNeurons[Lines,Columns,weights]		neuronen aus EsomNeurons
#       or: List containing EsomNeurons (result of esom())
# OPTIONAL
# Toroid				planar=F
# OUTPUT
# Umatrix[Lines,Columns]

# Editor: FL - 16.08.2016 - Statt EsomNeurons kann auch Liste uebergeben werden
#############################################
  
 # intern wird mit EsomNeurons gearbeitet (3 dimensionale arrays)
 EsomNeurons = ListAsEsomNeurons(Weights, Lines, Columns)
  
## Nachbarn()
nachbarn <- function(k, EsomNeurons, Toroid=FALSE){
# INPUT
# k Gitterpunkt
# EsomNeurons[Lines,Columns,weights]
# Toroid
# OUTPUT
# nb
M <- dim(EsomNeurons)[1]
N <- dim(EsomNeurons)[2]
if(Toroid){
	pos1 = c(k[1]-1,k[1]-1,k[1]-1,k[1],k[1],k[1]+1,k[1]+1,k[1]+1) %% M
	pos1[which(pos1==0)] = M
	pos2 = c(k[2]-1,k[2],k[2]+1,k[2]-1,k[2]+1,k[2]-1,k[2],k[2]+1) %% N
	pos2[which(pos2==0)] = N
	nb = cbind(pos1,pos2)
}else{# planar
	 if(k[1] == 1){
		if (k[2] == 1){
			nb = rbind(c(1,2), c(2,2), c(2,1))
		}else{
			if (k[2] == N){
				nb = rbind(c(1,(N-1)), c(2,(N-1)), c(2,N))
			}else{
				nb = rbind(c(1,(k[2]-1)), c(1,(k[2]+1)), c(2,(k[2]-1)), c(2,k[2]), c(2,(k[2]+1)))
			}
		}
	}#end fall 1 fuer planar
	if(k[1] == M){
		if(k[2] == 1){
			nb = rbind(c((M-1),1), c((M-1),2), c(M,2))
		}else{
			if(k[2] == N){
				nb = rbind(c((M-1),(N-1)), c((M-1),N), c(M,(N-1)))
			}else{
				nb = rbind(c((M-1),(k[2]-1)), c((M-1),k[2]), c((M-1),(k[2]+1)), c(M,(k[2]-1)), c(M,(k[2]+1)))
			}
		}
	}#end fall 2 fuer planar
	if(k[1] != 1 && k[1] != M){
		if(k[2] == 1){
							nb = rbind(c((k[1]-1),1), c((k[1]-1),2), c(k[1],2),c((k[1]+1),1), c((k[1]+1),2))
		}else{
			if(k[2] == N){
					nb = rbind(c((k[1]-1),(N-1)), c((k[1]-1),N), c(k[1],(N-1)), c((k[1]+1),(N-1)), c((k[1]+1),N))
			}else{
					nb = rbind(c((k[1]-1),(k[2]-1)), c((k[1]-1),k[2]), c((k[1]-1),(k[2]+1)), c(k[1],(k[2]-1)),c(k[1],(k[2]+1)), c((k[1]+1),(k[2]-1)), c((k[1]+1),k[2]), c((k[1]+1),(k[2]+1)))
			}
	 }
	}#end fall 3 fuer planar
}# end if Toroid
return(nb)
}
############################################
if(is.list(EsomNeurons)) EsomNeurons = EsomNeurons$EsomNeurons  # FL

k = dim(EsomNeurons)[1]
m = dim(EsomNeurons)[2]
Umatrix = matrix(0,k,m)

d=dim(EsomNeurons)[3]
if(is.null(d)){#wts als liste
  stop('EsomNeurons wts has to be an array[1:Lines,1:Columns,1:Weights], use ListAsEsomNeurons')
}

for(i in 1:k){
	for(j in 1:m){
		nbs=nachbarn(c(i,j),EsomNeurons,Toroid)
		wij=EsomNeurons[i,j,]
		n.nbs=dim(nbs)[1]
		for(l in 1:n.nbs){
			nij=EsomNeurons[nbs[l,1],nbs[l,2],]
			Umatrix[i,j]=Umatrix[i,j]+sqrt(sum((wij-nij)^2))
		}
		Umatrix[i,j]=Umatrix[i,j]/n.nbs
	}
}
return(Umatrix)
}

