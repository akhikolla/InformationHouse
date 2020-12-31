
#impute the data to get a gaussian distribution in each column

impute.igcda=function (tab, tab.imp, conditions, q = 0.95) 
{
  dataSet.imputed = tab
  nb_cond = length(levels(conditions))
  nb_rep = rep(0, nb_cond)
  k = 1
  j = 1
  jind = NULL
  for (i in 1:nb_cond) {
    nb_rep[i] = sum((conditions == levels(conditions)[i]))
    dataSet.mvs = tab[, (k:(k + nb_rep[i] - 1))]
    mim = suppressWarnings(apply(dataSet.mvs, 1, min, na.rm = TRUE))
    mam = suppressWarnings(apply(dataSet.mvs, 1, max, na.rm = TRUE))
    rmammim = mam - mim
    rangem = quantile(rmammim, probs = q)
    mini = min(mim, na.rm = T)
    while (j < (k + nb_rep[i])) {
      curr.sample = tab[, j]
      if (sum(is.na(curr.sample)) != 0) {
        Fobs = ecdf(curr.sample)
        if (sum(tab.imp[which(is.na(curr.sample)), j],na.rm = T) != 0) {
          sta = Fobs(quantile(tab.imp[which(is.na(curr.sample)),j], 0.9, na.rm = T))
          gamma = sum(is.na(curr.sample))/length(curr.sample)
          if (sta>0.99){sta=0.99;}
          interv = gamma + (1 - gamma) * seq(sta, 0.999, 
                                             length.out = 100)
          q.normal = qnorm(interv, mean = 0, sd = 1)
          q.curr.sample = quantile(curr.sample, probs = seq(sta, 
                                                            0.999, length.out = 100), na.rm = T)
          q.curr.sample[which(!is.finite(q.curr.sample))] = NA
          q.normal[which(!is.finite(q.normal))] = NA
          temp.QR = lm(q.curr.sample ~ q.normal, na.action=na.omit)
          m = temp.QR$coefficients[1]
          v = (as.numeric(temp.QR$coefficients[2]))^2
          dx = seq(mini - rangem, max(curr.sample, na.rm = T), 
                   length.out = 1000)
          Fn = pnorm(dx, mean = m, sd = sqrt(v))
          Fna = (Fn - (1 - gamma) * Fobs(dx))/gamma
          Fna[Fna > 1] = 1
          Fna = pava(Fna)
          pna = diff(Fna)
          gen.sample = NULL
          for (ll in 1:length(pna)) {
            gen.sample = c(gen.sample, runif(floor(pna[ll] * 
                                                     max(10000, sum(is.na(curr.sample)) + 1)), 
                                             dx[ll], dx[ll + 1]))
          }
          curr.sample.imputed = curr.sample
          gss = gen.sample[floor(runif(sum(is.na(curr.sample)), 
                                       1, length(gen.sample) + 1))]
          curr.sample.imputed[which(is.na(curr.sample))][order(mim[which(is.na(curr.sample))])] = sort(gss)
          dataSet.imputed[, j] = curr.sample.imputed
        }
        else {
          cat(paste("Warning: in sample ", j, ": the imputed sample is equal to the original sample. impute.slsa has been used."))
          jind = c(jind, j)
        }
      }
      j = j + 1
    }
    k = k + nb_rep[i]
  }
  if (!is.null(jind)) {
    dataSet.imputed = impute.slsa(dataSet.imputed, conditions)
  }
  return(dataSet.imputed)
}