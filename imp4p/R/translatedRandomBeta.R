
translatedRandomBeta=function (n, min, max, param1 = 3, param2 = 1){
  scale <- max - min
  simu <- rbeta(n, param1, param2)
  res <- (simu * scale) + min
  return(res)
}
