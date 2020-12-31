getnodes =  function(qrf, data) {
  class(qrf) = "ranger"
  nodes = stats::predict(qrf, data = data, predict.all = TRUE)$predictions
  return(nodes)
}
