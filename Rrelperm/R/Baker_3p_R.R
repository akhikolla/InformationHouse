func <- function(s, SW3t, SG3t, SWCON, SORW, SORG, SGCON, KROCW, KROGCG, NOW, NOG) {
   f <- vector(length = length(s))
   SW <- s[1]^2
   SG <- s[2]^2
   f[1] <- krow2p_BC(SW, SWCON, SORW, KROCW, NOW) - krgl2p_BC(SG, SWCON, SORG, SGCON, KROGCG, NOG)
   f[2] <- (SW3t - SWCON) * (SG3t - SGCON) - (SW3t - SW) *  (SG3t - SG)
   # f[2] <- (SWCON - SW) * (SG3t - SGCON) - (SW3t - SW) *  (SG - SGCON)
   return(f)
}


kr3p_Baker_R <- function(SWCON, SWCRIT, SOIRW, SORW, SOIRG, SORG, SGCON, SGCRIT, KRWIRO, KROCW, KRGCL, NW, NOW, NG, NOG, NP) {

   threshold <- 1e-6
   len <- NP
   if (len > 501) len <- 501
   lent <- len * (len + 1) / 2
   KROGCG <- KROCW
   counter_1 <- 0
   counter_2 <- 0
   SW <- seq(0, 1, by = 1 / (len - 1))
   SG <- seq(0, 1, by = 1 / (len - 1))
   SLRG <- SWCON + SORG
   mat2 <- matrix(NA, nrow = lent, ncol = 4)
   m <- (SORG - (SORW - SGCON)) / (SWCON - (1 - SORW))
   b <- (SORG) - m * SWCON
   for (i in 1:len) {
      for (j in 1:(len - i + 1)) {
         mat2[counter_2 + j,1] <- SW[i]
         mat2[counter_2 + j,2] <- SG[j]
         mat2[counter_2 + j,3] <- 1 - SW[i] - SG[j]
         SW3t <- SW[i]
         SG3t <- SG[j]
         SO3t <- 1 - SW[i] - SG[j]
         if ((SW3t > SWCON) & (SW3t < (1 - SORW)) & (SG3t < (1 - SORG - SWCON)) & (SG3t > SGCON)) {
            SOr <- m * SW3t + b
            if (SO3t <= SOr) {
               mat2[counter_2 + j,4] <- 0
            } else {
               s <- sqrt(c(0.66 * (SWCON + (1 - SORW)), 0.66 * (SGCON + (1 - SORG))))
               res <- nleqslv::nleqslv(x = s, fn = func, SW3t = SW3t, SG3t = SG3t, SWCON = SWCON,
                                       SORW = SORW, SORG = SORG, SGCON = SGCON, KROCW = KROCW,
                                       KROGCG = KROGCG, NOW = NOW, NOG = NOG, method = c("Broyden"))
               SW2 <- res$x[1] * res$x[1]
               mat2[counter_2 + j,4] <- krow2p_BC(SW2, SWCON, SORW, KROCW, NOW)
            }
         } else if (SW3t <= SWCON) {
            mat2[counter_2 + j,4] <- krgl2p_BC(SG3t, SWCON, SORG, SGCON, KROGCG, NOG)
         } else if (SW3t >= (1 - SORW)) {
            mat2[counter_2 + j,4] <- krow2p_BC(SW3t, SWCON, SORW, KROCW, NOW)
         } else if (SG3t <= SGCON) {
            mat2[counter_2 + j,4] <- krow2p_BC(SW3t, SWCON, SORW, KROCW, NOW)
         } else {
            mat2[counter_2 + j,4] <- krgl2p_BC(SG3t, SWCON, SORG, SGCON, KROGCG, NOG)
         }
         if (mat2[counter_2 + j,4] <= threshold) {
            mat2[counter_2 + j,4] <- 0
         }
      }
      counter_1 <- len - i + 1
      counter_2 <- counter_2 + counter_1
   }
   colnames(mat2) <- c("Sw", "Sg", "So", "Kro")
   return(mat2);
}
