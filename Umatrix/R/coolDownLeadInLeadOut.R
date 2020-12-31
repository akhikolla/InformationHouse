coolDownLeadInLeadOut<-function(start, end, steps, step){
# res <- coolDownLeadInLeadOut(start, end, steps, step)
# divides range from start to end into equally sized steps and returns current value
#  the first 10% of the steps (or epochs) always return the start value
#  and the last 5% of the steps (or epochs) always return the end value
# INPUT
# start		  start value of linear cooling
# end		  end value of linear cooling
# steps		  number of steps in which the range from start to end should be split
# step		  the current step
# OUTPUT
# cooledDownVal	  value after current step
# author: Florian Lerch


  pctg = step/steps

  if(pctg <= 0.1) return(start)
  else if(pctg >= 0.95) return(end)
  else return(start-((start-end)/steps)*(step-1))


}
