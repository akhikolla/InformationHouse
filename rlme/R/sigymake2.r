sigymake2 <- function(I, sec, mat, siga2, sige2) {
    ss = apply(mat, 1, sum)
    n = sum(ss)
    iflag = 0
    for (i in 1:I) {
        for (j in 1:sec[i]) {
            mattmp2 = jmake(mat[i, j])
            if (j == 1) {
                secpart = mattmp2
            }
            else {
                secpart = bdmake(secpart, mattmp2)
            }
        }
        mattmp = siga2 * jmake(ss[i]) + sige2 * diag(rep(1, ss[i]))
        
        eigs = eigen(mattmp, symmetric = T)
        
        if (min(eigs$values) < 0) {
            iflag = 1
        }
        
        #mattmp3 = sigma12i(mattmp)
        
        mattmp3 = eigs$vectors %*% diag(1/eigs$values^0.5) %*% 
          t(eigs$vectors)
        
        mattmp4 = eigs$vectors %*% diag(eigs$values^0.5) %*% 
          t(eigs$vectors)
        
        if (i == 1) {
            schpart = mattmp
            schpart12i = mattmp3
            siggma = mattmp4
        }
        else {
            schpart = bdmake(schpart, mattmp)
            schpart12i = bdmake(schpart12i, mattmp3)
            siggma = bdmake(siggma, mattmp4)
        }
    }
    list(sigy2 = schpart, sigy12i = schpart12i, iflag = iflag, siggma = siggma)
}
