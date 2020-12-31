gif.generate<-function (M.save, M.source, M.target, K, file.name, fps, 
                        new.l, gif_type,out.col= grey(0:1000/1000),width=800,height=800) 
{
  if (!requireNamespace("animation", quietly = TRUE)) {
    stop("Package 'animation' required. Please install it.")
  }
  dirname <- "temp"
  dirno <- 1
  gif.length <- K
  while (!dir.create(dirname, showWarnings = FALSE)) {
    dirname <- paste("temp", dirno, sep = "")
    dirno <- dirno + 1
  }
  gif.pot <- 2
  while (gif.length > 10) {
    gif.pot <- gif.pot + 1
    gif.length <- gif.length/10
  }
  count <- 1
  zero.drop <- 0
  pot.check <- 10
  zero.string <- ""
  for (z in 1:(gif.pot)) {
    zero.string <- paste("0", zero.string, sep = "")
  }
  for (i in seq(1, 10)) {
    if (count == pot.check) {
      zero.drop <- zero.drop + 1
      pot.check <- pot.check * 10
      zero.string <- ""
      for (z in 1:(gif.pot - zero.drop)) {
        zero.string <- paste("0", zero.string, sep = "")
      }
    }
    name.number <- paste(zero.string, count, sep = "")
    name <- paste(dirname, "/output", name.number, ".png", 
                  sep = "")
    png(filename = name, width = width, height = height)
    image(M.source, col = out.col, axes = FALSE, 
          zlim = c(min(M.save), max(M.save)))
    dev.off()
    count <- count + 1
  }
  for (i in c(1:K)) {
    if (count == pot.check) {
      zero.drop <- zero.drop + 1
      pot.check <- pot.check * 10
      zero.string <- ""
      for (z in 1:(gif.pot - zero.drop)) {
        zero.string <- paste("0", zero.string, sep = "")
      }
    }
    name.number <- paste(zero.string, count, sep = "")
    name <- paste(dirname, "/output", name.number, ".png", 
                  sep = "")
    png(filename = name, width = width, height = height)
    image(M.save[, , i], col = out.col, axes = FALSE, 
          zlim = c(min(M.save), max(M.save)))
    dev.off()
    count <- count + 1
  }
  for (i in seq(1, 20)) {
    if (count == pot.check) {
      zero.drop <- zero.drop + 1
      pot.check <- pot.check * 10
      zero.string <- ""
      for (z in 1:(gif.pot - zero.drop)) {
        zero.string <- paste("0", zero.string, sep = "")
      }
    }
    name.number <- paste(zero.string, count, sep = "")
    name <- paste(dirname, "/output", name.number, ".png", 
                  sep = "")
    png(filename = name, width = width, height = height)
    image(M.target, col = out.col, axes = FALSE, 
          zlim = c(min(M.save), max(M.save)))
    dev.off()
    count <- count + 1
  }
  oopt = animation::ani.options(interval = 1/fps,ani.width=width,ani.height=height)
  #animation::ani.options(oopt)
  files_tmp<-list.files(path=dirname,pattern = "*.png")
  for (k in 1:length(files_tmp)){
    files_tmp[k]<-paste(dirname,"/",files_tmp[k],sep="")
  }
  if (gif_type == "gif_im") {
    animation::im.convert(files = files_tmp, output = file.name)
  }
  #animation::ani.options(oopt)
  unlink(dirname, recursive = TRUE)
  return(invisible())
}
transport_track<-function (source, target, tplan, K = 50, scmult = 1, smooth = FALSE,
                           H = matrix(c(1,0,0,1),2,2), create.file = c("none",
                           "gif_im"), file.name = "Rtransport.gif", fps = 20,cut = FALSE, col=grey((0:1000)/1000),width=800,height=800){
  if (!requireNamespace("ks", quietly = TRUE)) {
    stop("Package 'ks' required. Please install it.")
  }
  create.file <- match.arg(create.file)
  M.source.n <- source$mass
  M.target.n <- target$mass
  l <- sqrt(length(M.source.n))
  new.l <- l * scmult
  M.source <- matrix(0, new.l, new.l)
  M.target <- matrix(0, new.l, new.l)
  M.save <- array(0, c(new.l, new.l, (K)))
  for (i in 1:l) {
    for (j in 1:l) {
      M.source[((i - 1) * scmult + 1):(i * scmult), ((j - 1) * scmult + 
        1):(j * scmult)] <- matrix(M.source.n[i,j]/scmult/scmult, scmult, scmult)
      M.target[((i - 1) * scmult + 1):(i * scmult), ((j - 1) * scmult +
        1):(j * scmult)] <- matrix(M.target.n[i,j]/scmult/scmult, scmult, scmult)
    }
  }
  m <- length(tplan$from)
  arrows.from <- rep(0, m * scmult * scmult)
  arrows.to <- rep(0, m * scmult * scmult)
  arrows.mass <- rep(0, m * scmult * scmult)
  for (k in 1:m) {
    count <- 1
    from <- tplan$from[k]
    to <- tplan$to[k]
    mass <- tplan$mass[k]
    arrows.from.insert <- rep(0, scmult * scmult)
    arrows.to.insert <- rep(0, scmult * scmult)
    to.col <- (ceiling(to/l) - 1) * scmult
    r <- to%%l
    if (r == 0) {
      r <- l
    }
    to.row <- (r - 1) * scmult + 1
    to.corner <- new.l * to.col + to.row
    from.col <- (ceiling(from/l) - 1) * scmult
    r <- from%%l
    if (r == 0) {
      r <- l
    }
    from.row <- (r - 1) * scmult + 1
    from.corner <- new.l * from.col + from.row
    for (i in 1:scmult) {
      for (j in 1:scmult) {
        arrows.from.insert[count] <- from.corner + (j -
                                                      1) + (i - 1) * new.l
        arrows.to.insert[count] <- to.corner + (j - 1) +
          (i - 1) * new.l
        count <- count + 1
      }
    }
    arrows.to[((k - 1) * scmult * scmult + 1):(k * scmult *
      scmult)] <- arrows.to.insert
    arrows.from[((k - 1) * scmult * scmult + 1):(k * scmult *
      scmult)] <- arrows.from.insert
    arrows.mass[((k - 1) * scmult * scmult + 1):(k * scmult *
      scmult)] <- rep(mass/scmult/scmult, scmult * scmult)
  }
  m.new <- m * scmult * scmult
  alpha <- 1/(K+1)
  source <- matrix(0, m.new, 2)
  end <- matrix(0, m.new, 2)
  dist <- matrix(0, m.new, 2)
  for (i in (1:m.new)) {
    s.col <- ceiling(arrows.from[i]/new.l)
    s.row <- arrows.from[i]%%new.l
    if (s.row == 0) {
      s.row <- new.l
    }
    e.col <- ceiling(arrows.to[i]/new.l)
    e.row <- arrows.to[i]%%new.l
    if (e.row == 0) {
      e.row <- new.l
    }
    source[i, 1] <- s.col
    source[i, 2] <- s.row
    end[i, 1] <- e.col
    end[i, 2] <- e.row
    dist[i, 1] <- e.col - s.col
    dist[i, 2] <- e.row - s.row
  }
  eval <- matrix(0, new.l * new.l, 2)
  for (i in 1:new.l) {
    eval[((i - 1) * new.l + 1):(i * new.l), 1] <- t(1:new.l) -
      0.5
    eval[((i - 1) * new.l + 1):(i * new.l), 2] <- rep(i -
                                                        0.5, new.l)
  }
  M.setup <- M.source
  for (i in (1:m.new)) {
    M.setup[arrows.from[i]] <- M.setup[arrows.from[i]] -
      arrows.mass[i]
  }
  weight <- arrows.mass
  weight<-rbind(matrix(weight,length(weight),1),matrix(M.setup,length(M.setup),1))
  v.temp<-rep(0,new.l*new.l)
  for (k in 1:new.l){
    v.temp[((k-1)*new.l+1):(new.l*k)]<-rep(k,new.l)
  }
  grid.m<-cbind(v.temp,rep(1:new.l,new.l))
  grid.m<-grid.m-0.5
  if (smooth == FALSE) {
    source.kernel <- ks::binning(grid.m/new.l, w = matrix(M.source,length(M.source),1),
                             xmin = c(0,0), xmax = c(1, 1), bgridsize = c(new.l, new.l),H=H)
    est <- t(matrix(source.kernel$counts, new.l))
    M.source.kernel <-est
    end.kernel <- ks::binning(grid.m/new.l, w = matrix(M.target,length(M.source),1),
                          xmin = c(0,0), xmax = c(1, 1), bgridsize = c(new.l, new.l),H=H)
    est <- t(matrix(end.kernel$counts, new.l))
    M.end.kernel <- est
    for (k in (1:K)) {
      cat(k, ",", sep = "")
      m1 <- source[, 1] + (k * alpha * dist[, 1])
      m2 <- source[, 2] + (k * alpha * dist[, 2])
      mean.new <- cbind(m1, m2) - 0.5
      mean.new <- rbind(mean.new,grid.m)
      kernel <- ks::binning(mean.new/new.l, w = weight, xmin = c(0,
                                                             0), xmax = c(1, 1), bgridsize = c(new.l, new.l),H=H)
      p <- sum(kernel$counts)
      est <- t(matrix(kernel$counts, new.l))
      M.save[, , k] <- est
    }
  }
  else {
     weight<-weight/sum(weight)*(length(source[,1])+length(grid.m[,1]))
     source.kernel <- ks::kde(grid.m, w = (matrix(M.source,length(M.source),1)/sum(M.source)*length(M.source)),
                          eval.points = eval,H=H)
     est <- t(matrix(source.kernel$estimate, new.l))
     M.source.kernel <- est
     end.kernel <- ks::kde(grid.m, w = (matrix(M.target,length(M.target),1)/sum(M.source)*length(M.source)),
                       eval.points = eval,H=H)
     est <- t(matrix(end.kernel$estimate, new.l))
     M.end.kernel <- est
    for (k in (1:K)) {
      cat(k, ",", sep = "")
      m1 <- source[, 1] + (k * alpha * dist[, 1])
      m2 <- source[, 2] + (k * alpha * dist[, 2])
      mean.new <- cbind(m1, m2) - 0.5
      mean.new <- rbind(mean.new,grid.m)
      kernel <- ks::kde(mean.new, w = weight, eval.points = eval,H=H)
      p <- sum(kernel$estimate)
      est <- t(matrix(kernel$estimate, new.l))
      M.save[, , k] <- est
    }
  }
  cat("\n")
  if (cut==TRUE){
    M.save<-M.save[2:(new.l-1),2:(new.l-1),]
    M.source.kernel<-M.source.kernel[2:(new.l-1),2:(new.l-1)]
    M.end.kernel<-M.end.kernel[2:(new.l-1),2:(new.l-1)]
    new.l<-new.l-1
  }
   M.source <- M.source.kernel
   M.target <- M.end.kernel
  if ((create.file == "gif_im")){
  	if (!requireNamespace("animation", quietly = TRUE)) {
      stop("Package 'animation' required for creation of animated gif. Please install it or use option create.file='none'")
    }
    gif.generate(M.save,M.source,M.target,K,file.name,fps,new.l,create.file,out.col=col,width=width,height=height)
    return(invisible(M.save))
  } else {
    return(M.save)
  }
}
