sort_aux <- function(A) {
  output <- rep(0,13)

  for (k in seq_along(output)) {
    start <- 2^(k - 1)
    end <- -1 + 2^(k)
    if(end > length(A)) {
      end <- length(A)
      X <- sum(A[start:end])
      output[k] <- X
      break
    }
    X <- sum(A[start:end])
    output[k] <- X
  }
  return(output)
}

expected.SAD <- function(theta, m, J) {
  if (theta < 1) {
    stop("expected.SAD: ",
         "theta can not be below one")
  }
  if (m < 0) {
    stop("expected.SAD: ",
         "m can not be below zero")
  }
  if (m > 1) {
    stop("expected.SAD: ",
         "m can not be above 1")
  }
  if (J < 0) {
    stop("expected.SAD: ",
         "J can not be below zero")
  }

  I <- (J - 1) * m / (1 - m)
  aux <- pm_sad(theta, I, J)
  sad <- sort_aux(aux)
  return(sad)
}

expected.SAD.Guilds <- function(theta, alpha_x, alpha_y,
                                J, n_replicates = 100) {
  if (theta < 1) {
    stop("expected.SAD.Guilds: ",
         "theta can not be below one")
  }
  if (alpha_x < 0) {
    stop("expected.SAD.Guilds: ",
         "alpha_x can not be below zero")
  }
  if (alpha_x > 1) {
    stop("expected.SAD.Guilds: ",
         "alpha_x can not be above 1")
  }
  if (alpha_y < 0) {
    stop("expected.SAD.Guilds: ",
         "alpha_y can not be below zero")
  }
  if (alpha_y > 1) {
    stop("expected.SAD.Guilds: ",
         "alpha_y can not be above 1")
  }
  if (J < 1) {
    stop("expected.SAD: ",
         "J can not be below one")
  }


  meanx <- rep(0, J)
  meany <- rep(0, J)

  for (r in seq_len(n_replicates)) {
		M <- draw_local(theta, alpha_x, alpha_y, J)
		for (m in seq_along(M$guildX)) {
			meanx[m] <- meanx[m] + M$guildX[m]
		}
		for (m in seq_along(M$guildY)) {
			meany[m] <- meany[m] + M$guildY[m]
		}
  }
  meanx <- meanx / n_replicates
  meany <- meany / n_replicates

  gx <- sort_aux(meanx)
  gy <- sort_aux(meany)
  
  output <- list( guildX = gx, guildY = gy)
  return(output)
}

expected.SAD.Guilds.Conditional <- function(theta,
                                            alpha_x,
                                            alpha_y,
                                            Jx,
                                            Jy,
                                            n_replicates = 100) {
  if (theta < 1) {
    stop("expected.SAD.Guilds.Conditional: ",
         "theta can not be below one")
  }
  if (alpha_x < 0) {
    stop("expected.SAD.Guilds.Conditional: ",
         "alpha_x can not be below zero")
  }
  if (alpha_x > 1) {
    stop("expected.SAD.Guilds.Conditional: ",
         "alpha_x can not be above 1")
  }
  if (alpha_y < 0) {
    stop("expected.SAD.Guilds.Conditional: ",
         "alpha_y can not be below zero")
  }
  if (alpha_y > 1) {
    stop("expected.SAD.Guilds.Conditional: ",
         "alpha_y can not be above 1")
  }
  if (Jx < 1) {
    stop("expected.SAD.Guilds.Conditional: ",
         "Jx can not be below one")
  }
  if (Jy < 1) {
    stop("expected.SAD.Guilds.Conditional: ",
         "Jy can not be below one")
  }

  meanx <- rep(0,Jx)
  meany <- rep(0,Jy)

  for (r in seq_len(n_replicates)) {
		M <- draw_local_cond(theta, alpha_x, alpha_y, Jx, Jy)
		for (m in seq_along(M$guildX)) {
			meanx[m] <- meanx[m] + M$guildX[m]
		}
		for (m in seq_along(M$guildY)) {
			meany[m] <- meany[m] + M$guildY[m]
		}
  }
  meanx <- meanx / n_replicates
  meany <- meany / n_replicates

  gx <- sort_aux(meanx)
  gy <- sort_aux(meany)
  
  output <- list( guildX = gx, guildY = gy)
  return(output)
}