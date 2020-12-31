rprmeddis2 <- function(I, sec, mat, ehat, location, scale, rprpair = "") {
    rprpair = tolower(rprpair)
    if (rprpair == "hl-disp") {
        location = scale = 2
    }
    if (rprpair == "med-mad") {
        location = scale = 1
    }
    if (location == 1) {
        matre = matrix(c(0, 0), ncol = 2)
        sch_size = as.vector(apply(mat, 1, sum))
        ahati = rep(0, I)
        epshatij = rep(0, sum(sch_size))
        ic = 0
        id = 0
        for (i in 1:I) {
            ehattmp = ehat[(ic + 1):(ic + sch_size[i])]
            ic = ic + sch_size[i]
            ahati[i] = median(ehattmp)
        }
        amed = median(ahati)
        ahati = ahati - amed
        for (j in 1:I) {
            for (k in 1:sch_size[j]) {
                epshatij[(id + k)] = ehat[(id + k)] - ahati[j]
            }
            id = id + sch_size[j]
        }
        epshatij = epshatij - median(epshatij)
    }
    if (location == 2) {
        matre = matrix(c(0, 0), ncol = 2)
        sch_size = as.vector(apply(mat, 1, sum))
        ahati = rep(0, I)
        epshatij = rep(0, sum(sch_size))
        ic = 0
        id = 0
        for (i in 1:I) {
            ehattmp = ehat[(ic + 1):(ic + sch_size[i])]
            ic = ic + sch_size[i]
            if (length(ehattmp) == 1) {
                ahati[i] = ehattmp
            }
            else {
                ahati[i] = onesampwil(ehattmp, maktable = F, 
                  plotb = F)$est
            }
        }
        amed = onesampwil(ahati, maktable = F, plotb = F)$est
        ahati = ahati - amed
        for (j in 1:I) {
            for (k in 1:sch_size[j]) {
                epshatij[(id + k)] = ehat[(id + k)] - ahati[j]
            }
            id = id + sch_size[j]
        }
        epshatij = epshatij - onesampwil(epshatij, maktable = F, 
            plotb = F)$est
    }
    if (scale == 1) {
        siga2 = mad(ahati)^2
        sigmae2 = mad(epshatij)^2
    }
    if (scale == 2) {
        siga2 = dispvar(ahati)^2
        sigmae2 = dispvar(epshatij)^2
    }
    list(frei = ahati, free = epshatij, siga2 = siga2, sigmae2 = sigmae2)
}
