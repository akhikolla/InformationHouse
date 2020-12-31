CalcUstarmatrix=CalcUstarmatrix=function(Umatrix,Pmatrix){
# ustarMatrix=calcUstarmatrix(Umatrix,Pmatrix)
# Berechnet die U*-Matrix aus der U-Matrix und der P-Matrix.
#
# INPUT
# Umatrix(1:Lines,1:Columns)                       Umatrix d.h. lokale Distanzen an jedem Punkt der ESOM
# Pmatrix(1:Lines,1:Columns)                       Pmatrix d.h. lokale Dichten   an jedem Punkt der ESOM
#
# OUTPUT
# list with
# UStarMatrix(1:Lines,1:Columns)                   U*matrix == Umatrix.* MagnificationFactor(Pmatrix)
#
# author: MT 04/2015
# Korrektur am 11/2015
# Example

# calcUstarmatrix(Umatrix,Pmatrix)
if(!prod(dim(Umatrix)==dim(Pmatrix))){
	stop('Dimensions of Pmatrix and Umatrix are not equal')
}
# # Bestimmung der Steigung A
Pmean = median(Pmatrix, na.rm = T)
P95   = quantile(as.vector(Pmatrix),c(0.95))
denominator = P95-Pmean
if(denominator < 0.01){
	denominator = 1
}
A = -1/denominator
# Additionsfaktor
B = -A*P95
scaleFactor = A*Pmatrix+B
# Make sure the intervall is between 0 and 1
scaleFactor = pmax(scaleFactor,0, na.rm = T)
scaleFactor = pmin(scaleFactor,1, na.rm = T)

# U-stern Matrix Skalieren
Ustar = Umatrix * scaleFactor
# Dort wo Hohe Umatrix-hoehen sind, geringere Reskalierung durch Dichte-Skalierungsfaktor
meanUheigths <- mean(Umatrix, na.rm = T)
hiInd <- which(Umatrix > meanUheigths)
Ustar[hiInd] <- 0.5*Umatrix[hiInd]+0.5*Ustar[hiInd]
# Sicherheitscheck
# Normierung, sodass maximale Hohe nicht die Umatrix-Maximal-Hoehe uebersteigt
uMax = max(Umatrix)
rows = nrow(Umatrix)
cols = ncol(Umatrix)
Ustar = matrix(pmin(Ustar,uMax),rows,cols)



## folgende Version geht nicht
#   # CONSTANTS
# MAXFACTOR = 2;      # U-Matrix Heights are at the most augmented by this Factor
# MINFACTOR = 0;      # U-Matrix Heights are at least multiplied by this Factor

# LINEARES MODELL:
# U*Matrix = U-Matrix .* Factor
# die Faktor Matrix wird wie folgt bestimmt
# Factor(p) == 1 <=> p =  mean(UHeights);
# Factor(p) == 1 <=> p =  mprctile(Pmatrix(:),95);
# Pmean = median(Pmatrix, na.rm = T)          # median der P-Matrix Hoehen
# P95   = quantile(as.vector(Pmatrix),c(0.95))         # 95 percentil der  P-Matrix Hoehen
#
# # lineare Funktion:
# # Faktor(x) = A*x + B;
# # Bestimmung der Steigung A
# Nenner =  Pmean - P95  #
# if(abs(Nenner)  < 0.01 ){ # zu wenig Unterschied zwischen median und p95
#     A = 1;  # keine Aenderung der
#     B = 0;
# }else{
#     A = 1/Nenner
#     B =  -A* P95
# } # if Nenner  < 0.01
#
# # damit wird der Skalierungsfaktor bestimmt:
# Factor = A * Pmatrix + B
# # Begrenzen auf das intervall [MINFACTOR   MAXFACTOR]
# #Factor = max(Factor,MINFACTOR)
# #Factor = min(Factor,MAXFACTOR)
# Factor[Factor<=MINFACTOR]=MINFACTOR
# Factor[Factor>=MAXFACTOR]=MAXFACTOR



return(Ustar)
}
