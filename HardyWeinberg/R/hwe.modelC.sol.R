hwe.modelC.sol <- function(x,parinit=c(0.50,0.00,0.00),verbose=FALSE,tracing=1) {
  hwe.out <- solnp(parinit, fun = hwe.modelC.obj, eqfun = NULL, eqB = NULL,  
                   ineqfun = eqn3, ineqLB = c(0,0,0,-Inf,0,-Inf),
                   ineqUB = c(Inf,Inf,Inf,0,Inf,0), LB = c(0,-1,-1), UB = c(1,1,1), 
                   control = list(trace = tracing, tol = 1e-8, inner.iter = 1000), x)
  k <- hwe.out$pars
  H <- hwe.out$hessian
  return(list(k=k,H=H,hwe.out=hwe.out))
}
