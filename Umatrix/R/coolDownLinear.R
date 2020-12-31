coolDownLinear<-function(start, end, steps, step){
# res <- coolDownLinear(start, end, steps, step)
# divides range from start to end into equally sized steps and returns current value
# INPUT
# start		  start value of linear cooling
# end		  end value of linear cooling
# steps		  number of steps in which the range from start to end should be split
# step		  the current step
# OUTPUT
# cooledDownVal	  value after current step
# author: Florian Lerch
# coolDownLinear(start, end, steps, step)
  start-((start-end)/steps)*(step-1)
}
