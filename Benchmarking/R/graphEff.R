# $Id: graphEff.R 229 2020-07-04 13:39:18Z lao $

# Funktion til beregning af graf efficiens.  Beregning sker via
# bisection hvor der itereres mellem mulige og ikke-mulige loesninger
# et LP problem hvor venstreside er som in- og output orienteret
# efficiens. blot er foerste soejle erstattet af rene nuller, og G*X og
# (1/G)*Y optraeder paa hoejresiden. Minimering af 0 saa der blot soeges om
# der er en mulig loesning. Da soejlen for efficiens er bar 0'er vil
# justering af efficiens ud over G ikke ske, dvs. det er kun lambdaer
# der tilpasses for at se om der er en mulig loesning.

graphEff <- function(lps, X, Y, XREF, YREF, RTS, FRONT.IDX, rlamb, oKr, 
         param=param, TRANSPOSE=FALSE, SLACK=FALSE, FAST=FALSE, LP=FALSE,
         CONTROL=CONTROL) 
{
   m = dim(X)[2]  # number of inputs
   n = dim(Y)[2]  # number of outputs
   K = dim(X)[1]  # number of units, firms, DMUs
   Kr = dim(YREF)[1]  # number of units, firms, DMUs

if (LP)  {
    cat("m=",m, ", n=",n, ", K=",K, ", Kr=", Kr, "\n", sep="")
    flush.console()
}

   objval <- rep(NA,K)   # vector for the final efficiencies
   if ( FAST ) {
     lambda <- NULL
   } else {
      lambda <- matrix(NA, nrow=K, ncol=Kr) # lambdas one column per unit
   }
   set.column(lps, 1, rep(0,dim(lps)[1]))
   lpcontr <- lp.control(lps)
   tol <- lpcontr$epsilon["epsint"]
   lp.control(lps, timeout=5, verbose="severe")
   if (!missing(CONTROL)) set_control(lps, CONTROL)
   for ( k in 1:K)  {
      if ( LP )  { print(paste("Firm",k), quote=FALSE); flush.console()}
      # Lav bisection
      a <- 0
      b <- 2  # medfoerer start med G=1
      nIter <- 0
      gFundet <- FALSE
      xIset <- TRUE

      # Er G==1 en mulighed?
      G <- 1
      set.rhs(lps, c(-G*X[k,],Y[k,]/G), 1:(m+n))
      set.basis(lps, default=TRUE)
      status <- solve(lps)
      #if (status==7) { print(paste("For G=1, status =",status)); flush.console()}
      if (LP) { print(paste("For G=1, status =",status)); flush.console()}
      nIter <- nIter + 1
      if ( status == 0 )  {
        # G=1 er mulig; hvis G=1-tol ikke er mulig, er G=1 optimal loesning
         G <- 1 - sqrt(tol)
        if (LP) { print(paste("G korrigeret til", G, "; sqrt(tol) =", sqrt(tol))); flush.console()}
         set.rhs(lps, c(-G*X[k,],Y[k,]/G), 1:(m+n))
         #if (LP) write.lp(lps, filename="graphEff.lp")
         # Uden set.basis gaer 'solve' en sjaelden gang i en uendelig loekke
         set.basis(lps, default=TRUE)
         status <- solve(lps)
         #if (status==7) { print(paste("For G=1-eps, status =",status)); flush.console()}
         if (LP) { print(paste("For G=1-eps, status =",status)); flush.console()}
         nIter <- nIter + 1
         if ( status != 0 )  {
            # G=1 er mulig og G=1-eps er ikke-mulig ==> G==1
            G <- b <- 1
            gFundet <- TRUE
         }
      } else {
         if (LP) {warning("Firm outside technology set, firm =",k); flush.console()}
         # G=1 er ikke mulig; firm uden for teknology set saa G > 1; eller
         # der er skeet en fejl i solve.
         xIset <- FALSE
         # Bestem oevre graense
         b <- 2
         while (status != 0 && nIter < 50)  {
            set.rhs(lps, c(-b*X[k,], Y[k,]/b), 1:(m+n))
            set.basis(lps, default=TRUE)
            status <- solve(lps)
            if ( status==5 )  {
               set.basis(lps, default=TRUE)
               status <- solve(lps)
            }
            #if (status==7) { print(paste0("b G = ", G, ", status =",status)); flush.console()}
            nIter <- nIter + 1
            b <- b^2
         }
         # nedre graense
         if ( b > 2 )  { a <- sqrt(b) # 'b' blev oeget anden potens og forrige
                                              # vaerdi var sqrt(b) som saa er en mulig
                                              # nedre graense
          }  else  { a <- 1 }
      }

      status <- 0 # status kunne godt have en anden vaerdi og saa 
                  # ville naeste loekke ikke blive gennemloebet
      if ( !gFundet && xIset && status == 0 )  {
        # Find en nedre graense; find et interval under 1 som kan
        # bruges ved start af bisection.
         # G==1 er mulig og G er ikke fundet endnu
         dif <- .1
         i <- 1
         while ( status==0 && i < 10 )  {
            # Saet G til en mindre vaerdi saalaenge det er en mulig loesning.
            G <- 1 - i*dif
            set.rhs(lps, c(-G*X[k,],Y[k,]/G), 1:(m+n))
            set.basis(lps, default=TRUE)
            status <- solve(lps)
            if ( status==5 )  {
               set.basis(lps, default=TRUE)
               status <- solve(lps)
            }
            #if (status==7) { print(paste("Graense: G = ",G,"; status = ",status)); flush.console()}
            if (LP) { print(paste("G = ",G,"; status = ",status)); flush.console()}
            nIter <- nIter + 1
            i <- i+1
         }
         # enten er i==10 eller ogsaa er status!=0
         if ( i==10 )  {
            a <- 0
            b <- dif
         } else {
            a <- 1 - (i-1)*dif
            b <- 1 - (i-2)*dif
         }
      }


      # bisection loekke
      if (LP) { print(paste("Bisection interval: [",a,",",b,"]")); flush.console()}
      while ( !gFundet && b-a > tol && nIter < 50 )  {
        # if (LP) {cat("nIter =", nIter, "\n"); flush.console()}
         G <- (a+b)/2
         # if (LP) { print(paste("Bisect: G = ",G,"(",k,")")); flush.console()}
         set.rhs(lps, c(-G*X[k,],Y[k,]/G), 1:(m+n))
            set.basis(lps, default=TRUE)
            status <- solve(lps)
         #if (status==7) { print(paste("Bisect G = ",G,"(",k,"); status =",status)); flush.console()}
         if (LP) { print(paste("G = ",G,"(",k,"); status =",status)); flush.console()}
         if ( status == 0 ) {
            # loesning findes
            b <- G
         } else {
            a <- G
         }
         nIter <- nIter + 1
      }
      if (LP) {print(paste("nIter =",nIter,"; status =",status)); flush.console()}

      if ( status != 0 )  {
         # Hvis den sidste vaerdi af G ikke var mulig bruger vi den
         # oevre graense. Det er noedvendigt med en mulig loesning for at
         # kunne faa lambdaer og duale vaerdier.
         G <- b
          set.rhs(lps, c(-G*X[k,],Y[k,]/G), 1:(m+n))
          set.basis(lps, default=TRUE)
         status <- solve(lps)
         #if (status==7) { print(paste("Sidste G = ",G,"; status =",status)); flush.console()}
       }
      if (LP)  {
         print(paste("G = ",G,"(",k,"); status =",status))
         # print(rlamb)
         # print("Solution")
         # print(get.variables(lps))
         # print(lps)
         flush.console()
      }
      objval[k] <- G
      if ( LP && k == 1 )  print(lps)
      if ( !FAST )  {
         sol <- get.variables(lps)
         lambda[k,] <- sol[2:(1+Kr)]
      }
    if (LP && status==0) {
         print(paste("Objval, firm",k))
         print(get.objective(lps))
         # print("Solution/varaibles")
         # print(get.variables(lps))
         # print("Primal solution")
         # print(get.primal.solution(lps))
         # print("Dual solution:")
         # print(get.dual.solution(lps))
         flush.console()
      }
   }  # loop for each firm
    lp.control(lps, timeout=0, verbose="neutral")


   e <- objval
   e[abs(e-1) < sqrt(tol)] <- 1
   lambda[abs(lambda-1) < sqrt(tol)] <- 1   # taet ved 1
   lambda[abs(lambda) < sqrt(tol)] <- 0     # taet ved 0

   if ( FAST ) { 
      return(e)
      stop("Her skulle vi ikke kunne komme i 'dea'")
   }

   if ( length(FRONT.IDX)>0 )  {
      colnames(lambda) <- paste("L",(1:oKr)[FRONT.IDX],sep="")
   } else {
      colnames(lambda) <- paste("L",1:Kr,sep="")
   }

   primal <- dual <- NULL
   ux <- vy <- NULL

   if ( TRANSPOSE ) {
      lambda <- t(lambda)
   }

   oe <- list(eff=e, lambda=lambda, objval=objval, RTS=RTS,
              primal=primal, dual=dual, ux=ux, vy=vy, gamma=gamma,
              ORIENTATION="graph", TRANSPOSE=TRANSPOSE
              # ,slack=slack_, sx=sx, sy=sy
              , param=param 
              )
   class(oe) <- "Farrell"

   if ( SLACK ) {
      if ( TRANSPOSE )  { # Transponer tilbage hvis de blev transponeret
         X <- t(X)
         Y <- t(Y)
         XREF <- t(XREF)
         YREF <- t(YREF)
      }
      sl <- slack(X, Y, oe, XREF, YREF, FRONT.IDX, LP=LP)
      oe$slack <- sl$slack
      oe$sx <- sl$sx
      oe$sy <- sl$sy
      oe$lambda <- sl$lambda
      if (LP)  {
         print("slack fra slack:")
         print(sl$slack)
         print("slack efter slack:")
         print(oe$slack)
         flush.console()
      }
   }

   return(oe)
}

