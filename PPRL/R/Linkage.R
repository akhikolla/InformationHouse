# Call for DeterministicLinkage and ProbabilisticLinkage

DeterministicLinkage <- function(IDA, dataA, IDB, dataB, blocking = NULL, similarity){
  blA <- vector('character')
  blB <- vector('character')
  vars1 <- vector('character')
  vars2 <- vector('character')
  res <- NULL
  method <- character()
  ind_c0 <- vector('logical')
  ind_c1 <- vector('logical')
  lenNgram <- vector('integer')
  i <- 1
  if (!is.null(blocking)){
  blA <- dataA[,match(blocking('variable1'),colnames( dataA ))] #gets the number of the column in dataA
  blB <- dataB[,match(blocking('variable2'),colnames( dataB ))]
  blockingMethod_ <- blocking('method')

  }
  else {
    blA <- NULL
    blB <- NULL
    blockingMethod_ <- "0"
  }
  # just one similarity option
  if (class(similarity) != "list"){
    vars1 <- c(vars1, dataA[,match(similarity('variable1'),colnames( dataA ))])
    vars2 <- c(vars2, dataB[,match(similarity('variable2'),colnames( dataB ))])
    method <- c(method, similarity('method'))
    ind_c0 <- c(ind_c0,  similarity('ind_c0'))
    ind_c1 <- c(ind_c1,  similarity('ind_c1'))
    lenNgram <- c(lenNgram,  similarity('lenNgram'))
    res<-.DeterministicLinkagec(IDA_ = IDA, dataA_ = dataA[,match(similarity('variable1'),colnames( dataA ))], blockingdataA_ = blA,
                                IDB_ = IDB, dataB_ = dataB[,match(similarity('variable2'),colnames( dataB ))], blockingdataB_ = blB,
                                method_ = similarity('method'), blocking_ = blockingMethod_, threshold_ = similarity('threshold'), lenNgram_ = lenNgram,
                                ind_c0_= ind_c0, ind_c1_= ind_c1, counterSim = i)
    res <- as.data.frame(res,stringsAsFactors = FALSE)
    names(res)[3] <- paste(similarity('variable1'), "vs", similarity('variable2'), names(res)[3]) #rename columns
  }
  # if similarity is a list
  else if (class(similarity) == "list"){
     for(l in similarity){
    #   # (match(similarity('variable1'),colnames( dataA ))) gets the number of the column in dataA
       vars1 <-NULL
       vars2 <-NULL
       vars1 <- c(vars1, dataA[,match(l('variable1'),colnames( dataA ))])
       vars2 <- c(vars2, dataB[,match(l('variable2'),colnames( dataB ))])
       method <- c(method, l('method'))
       ind_c0 <- c(ind_c0,  l('ind_c0'))
       ind_c1 <- c(ind_c1,  l('ind_c1'))
       lenNgram <- c(lenNgram,  l('lenNgram'))
       res<-c(res,.DeterministicLinkagec(IDA_ = IDA, dataA_ = vars1,  blockingdataA_ = blA,
                 IDB_ = IDB,  dataB_ = vars2, blockingdataB_ = blB,
                 method_ = method, blocking_ = blockingMethod_, threshold_ = l('threshold'), lenNgram_ = lenNgram,
                 ind_c0_= ind_c0, ind_c1_= ind_c1, counterSim = i))
       i <- i+1
     }
    res <- as.data.frame(res, stringsAsFactors = FALSE)
    i<-0
    for(l in similarity){
      names(res)[3+i] <- paste(l('variable1'), "vs" ,l('variable2'),"\n", names(res)[3+i] )
     i <- i+1
       }
    }
  return(res)
}

ProbabilisticLinkage <- function(IDA, dataA, IDB, dataB,  blocking = NULL , similarity){
  blA <- vector('character')
  blB <- vector('character')
  res <- NULL
  em <- NULL
  method <- character()
  ind_c0 <- vector('logical')
  ind_c1 <- vector('logical')
  lenNgram <- vector('integer')
  m <- double()
  u <- double()
  p <- double()
  epsilon <- double()
  if (!is.null(blocking)){
    blA <- dataA[,match(blocking('variable1'),colnames( dataA ))] #gets the number of the column in dataA
    blB <- dataB[,match(blocking('variable2'),colnames( dataB ))]
    blockingMethod_ <- blocking('method')
    #print(blA)
  }
  else {
    blA <- NULL
    blB <- NULL
    blockingMethod_ <- "0"
  }
  # if similarity is a list
  if (class(similarity) == "list"){
    vars1 <- NULL
    vars2 <- NULL
    for(l in similarity){
      # (match(similarity('variable1'),colnames( dataA ))) gets the number of the column in dataA
      vars1 <- c(vars1, list(dataA[,match(l('variable1'),colnames( dataA ))]))
      vars2 <- c(vars2, list(dataB[,match(l('variable2'),colnames( dataB ))]))
      method <- c(method, l('method') )
      ind_c0 <- c(ind_c0,  l('ind_c0'))
      ind_c1 <- c(ind_c1,  l('ind_c1'))
      lenNgram <- c(lenNgram,  l('lenNgram'))
      m <- c(m, l('m'))
      u <- c(u, l('u'))
      p <- c(p, l('p'))
      epsilon <- c(epsilon, l('epsilon'))

    }
      res<-.ProbabilisticLinkagec(IDA_ = IDA, dataA_ = vars1,  blockingdataA_ = blA,
                                 IDB_ = IDB,  dataB_ = vars2, blockingdataB_ = blB,
                                 method_ = method, blocking_ = blockingMethod_,
                                 threshold_ = l('threshold'), lenNgram_ = lenNgram,
                                 ind_c0_= ind_c0, ind_c1_= ind_c1,
                                 m_ = m,  u_ = u, p_ = p, e = epsilon,
                                 upper = l('upper'), lower = l('lower'),
                                 jaroWeightFactor = l('jaroWeightFactor'))


  }
  # just one similarity option
  else{
    method <- c(method, similarity('method') )
    ind_c0 <- c(ind_c0,  similarity('ind_c0'))
    ind_c1 <- c(ind_c1,  similarity('ind_c1'))
    m <- c(m, similarity('m'))
    u <- c(u, similarity('u'))
    #print(as.vector(table(dataA[,match(similarity('variable1'),colnames( dataA ))])))
    res<-.ProbabilisticLinkagec(IDA_ = IDA, dataA_ = dataA[,match(similarity('variable1'),colnames( dataA ))], blockingdataA_ = blA,
              IDB_ = IDB, dataB_ = dataB[,match(similarity('variable2'),colnames( dataB ))], blockingdataB_ = blB,
              method_ = similarity('method'), blocking_ = blockingMethod_, threshold_ = similarity('threshold'),
              lenNgram_ = similarity('lenNgram'),
              ind_c0_= similarity('ind_c0'), ind_c1_= similarity('ind_c1'),
              m_ = similarity('m'),  u_ = similarity('u'), p_ = similarity('p'), e = similarity('epsilon'),
              upper = similarity('upper'), lower = similarity('lower'),
              jaroWeightFactor =similarity('jaroWeightFactor'))
  }
  return(res)
}
