context("check_motif")
### Testing if the check_motif function works the way we want it to work #####

all_mot_perm <- gen_all_mot_perm()

for(i in 1:17) {
  mot <- motifs[[i]]
  N <- nrow(mot)
  M <- ncol(mot)
  pr <- gtools::permutations(N, N, 1:N, repeats = FALSE)
  pc <- gtools::permutations(M, M, 1:M, repeats = FALSE)
  # pick a random permutation for rows and cols
  rr <- floor(stats::runif(1, min=1, max= nrow(pr) + 1))
  rc <- floor(stats::runif(1, min=1, max= nrow(pc) + 1))

  mot2 <- mot[pr[rr, ], ,drop = FALSE]
  mot3 <- mot2[, pc[rc, ], drop = FALSE]

  numb <- check_motif(mot3, all_mot_perm)

  expect(numb == i, failure_message = "failed")
}
