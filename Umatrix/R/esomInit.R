esomInit<-function(Data,Lines=50, Columns=82, method = "uni_min_max"){
# res <- esomInit(Data,50,82)
# Initializes a grid of weightvectors based on given Data (callled by esom function)
#
# INPUT
# Data(1:m,1:n)		      sampleData that will be used for initialization
#				      n datapoints with m attributes each
# OPTIONAL
# Lines				      Lines of the grid
# Columns			      Columns of the grid
#
# OUTPUT
# return: WeightVectors(1:m,1:n)              matrix that contains n weights with m components each
# author: Florian Lerch
# details: Reimplemented from Databionic Esom Tools (http://databionic-esom.sourceforge.net)
# esomInit(Data)
  initMethod = method
# initMethod			      name of the method that will be used to choose initializations
#				      Valid Inputs: uni_min_max  <- uniform distribution with minimum and maximum
#								    from Data
#						    uni_mean_2std <- uniform distribution based on mean and standard
#								     deviation of Data
#						    norm_mean_2std <- normal distribuation based on mean and standard
#								     deviation of Data

  # initialize the weight-vectors
  if(initMethod == "uni_min_max"){
    l = c()
    for(i in 1:ncol(Data)){
      l <- cbind(l,runif(Lines*Columns, min(Data[,i]), max(Data[,i])))
    }
    WeightVectors <- l
  }

  # else if(initMethod == "uni_min_max+ica"){
  #   l = c()
  #   for(i in 1:ncol(Data)){
  #     l <- cbind(l,runif(Lines*Columns, min(Data[,i]), max(Data[,i])))
  #   }
  #   WeightVectors <- l
  #
  #   res = ICA(Data)$ProjectedPoints
  #
  #   # Projektion auf Gitter anpassen
  #   res[,1] = res[,1] + (1-min(res[,1])) # beginn bei 1
  #   res[,1] = res[,1] * (Lines / max(res[,1]))
  #   res[,2] = res[,2] + (1-min(res[,2])) # beginn bei 1
  #   res[,2] = res[,2] * (Columns / max(res[,2]))
  #   res=round(res)
  #
  #   # Gitterpunkte mit hochdimensionalen Werten aus Data initialisieren
  #   for(i in 1:nrow(res)){
  #     WeightVectors[ (res[i,1]-1)*Columns + res[i,2], ] = Data[i,]
  #   }
  #
  # }

  else if(initMethod == "uni_mean_2std"){
    l = c()
    for(i in 1:ncol(Data)){
      l <- cbind(l,runif(Lines*Columns,
	mean(Data[,i]) - (2*sd(Data[,i])),
	mean(Data[,i]) + (2*sd(Data[,i]))))
    }
    WeightVectors <- l
  }

  # choose values out of normal distribution with mean and standard deviation from Data
  else if(initMethod == "norm_mean_2std"){
    l = c()
    for(i in 1:ncol(Data)){
      l <- cbind(l,rnorm(Lines*Columns,
	mean(Data[,i]), sd(Data[,i])))
    }
    WeightVectors <- l
  }

  # fill the weightvectors only with zeros
  else if(initMethod == "zero"){
    WeightVectors <- matrix(0,ncol=ncol(Data),nrow=Lines*Columns)
  }

  else stop("The initialization method is not recognized.")

  WeightVectors
}

