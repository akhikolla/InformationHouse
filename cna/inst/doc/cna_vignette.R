### R code from vignette source 'cna_vignette.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: data examples
###################################################
library(cna)
cna(d.educate)
cna(d.women)


###################################################
### code chunk number 3: data type (eval = FALSE)
###################################################
## cna(d.jobsecurity, type = "fs")
## fscna(d.jobsecurity)
## cna(d.pban, type = "mv") 
## mvcna(d.pban)


###################################################
### code chunk number 4: configTable1 (eval = FALSE)
###################################################
## configTable(d.women)
## mvct(d.pban) 


###################################################
### code chunk number 5: configTable2 (eval = FALSE)
###################################################
## pact.ct <- configTable(d.pacts, type = "fs", case.cutoff = 2)
## ct2df(pact.ct)


###################################################
### code chunk number 6: simul1 (eval = FALSE)
###################################################
## allCombs(c(2, 2, 2)) - 1 
## allCombs(c(3, 4, 5))
## full.ct("A + B*c")
## full.ct(6)
## full.ct(list(A = 1:2, B = 0:1, C = 1:4)) 


###################################################
### code chunk number 7: simul1 (eval = FALSE)
###################################################
## dat1 <- allCombs(c(2, 2, 2)) - 1 
## selectCases("A + B <-> C", dat1)
## selectCases("(h*F + B*C*k + T*r <-> G)*(A*b + H*I*K <-> E)")
## target <- randomCsf(full.ct(6))
## selectCases(target)


###################################################
### code chunk number 8: simul2 (eval = FALSE)
###################################################
## dat2 <- full.ct(list(EN = 0:2, TE = 0:4, RU = 1:4)) 
## selectCases1("EN=1*TE=3 + EN=2*TE=0 <-> RU=2", dat2, con = .75, cov = .75)


###################################################
### code chunk number 9: simul3 (eval = FALSE)
###################################################
## makeFuzzy(selectCases("Hunger + Heat <-> Run"), 
##           fuzzvalues = seq(0, 0.4, 0.05))


###################################################
### code chunk number 10: simul4 (eval = FALSE)
###################################################
## dat3 <- allCombs(c(3, 4, 5))
## dat4 <- selectCases("A=1*B=3 + A=3 <-> C=2", mvct(dat3))
## some(dat4, n = 10, replace = FALSE)
## some(dat4, n = 1000, replace = TRUE)


###################################################
### code chunk number 11: cons1 (eval = FALSE)
###################################################
## fscna(d.jobsecurity, ordering = list("JSR"), strict = TRUE,
##       con = 1, cov = 1)
## fscna(d.jobsecurity, ordering = list("JSR"), strict = TRUE,
##       con = .9, cov = .9)


###################################################
### code chunk number 12: cons2 (eval = FALSE)
###################################################
## fscna(d.jobsecurity, ordering = list("JSR"), strict = TRUE,
##       con = .75, cov = .75)


###################################################
### code chunk number 13: odering1
###################################################
dat.aut.1 <- d.autonomy[15:30, c("AU","EM","SP","CO")]
ana.aut.1 <- fscna(dat.aut.1, ordering = list(c("EM","SP","CO"), "AU"), 
  strict = TRUE, con = .9, cov = .9)
printCols <- c("condition", "consistency", "coverage")
csf(ana.aut.1)[printCols]


###################################################
### code chunk number 14: odering2
###################################################
ana.aut.2 <- fscna(dat.aut.1, ordering = list(c("EM","SP","CO"), "AU"), 
  strict = FALSE, con = .9, cov = .9)
csf(ana.aut.2)[printCols]


###################################################
### code chunk number 15: maxstep1 (eval = FALSE)
###################################################
## cna(d.volatile, ordering = list("VO2"), maxstep = c(3,4,10))
## vol1 <- cna(d.volatile, ordering = list("VO2"), maxstep = c(4,4,10))
## csf(vol1, n.init = 3000)


###################################################
### code chunk number 16: maxstep2 (eval = FALSE)
###################################################
## cna(d.volatile, ordering = list("VO2"), maxstep = c(8,10,40), 
##   suff.only = TRUE)


###################################################
### code chunk number 17: maxstep3
###################################################
ana.jsc.1 <- fscna(d.jobsecurity, ordering = list("JSR"), con=.9, cov=.85)
csf(ana.jsc.1)[printCols]
ana.jsc.2 <- fscna(d.jobsecurity, ordering = list("JSR"), con=.9, cov=.85,
  maxstep = c(3,5,12))
csf(ana.jsc.2)[printCols]


###################################################
### code chunk number 18: notcols
###################################################
ana.aut.3 <- fscna(dat.aut.1, ordering = list(c("EM","SP","CO"), "AU"),
  strict = FALSE, con = .88, cov = .82, notcols = c("AU", "EM")) 
csf(ana.aut.3)[printCols]


###################################################
### code chunk number 19: tab2c (eval = FALSE)
###################################################
## dat5 <- allCombs(c(2, 2, 2, 2, 2)) -1
## dat6 <- selectCases("(A + B <-> C)*(A*B + D <-> E)", dat5)
## set.seed(3) 
## tab2c <- makeFuzzy(dat6, fuzzvalues = seq(0, 0.4, 0.01))
## fscna(tab2c, con = .8, cov = .8, what = "mac")


###################################################
### code chunk number 20: details (eval = FALSE)
###################################################
## cna(d.educate, details = TRUE)
## cna(d.educate, details = c( "co", "cy"))


###################################################
### code chunk number 21: what (eval = FALSE)
###################################################
## cna(d.educate, what = "tm")
## cna(d.educate, what = "mac")
## cna(d.educate, what = "all")


###################################################
### code chunk number 22: vol1 (eval = FALSE)
###################################################
## vol2 <- cna(d.volatile, ordering = list("VO2"))
## msc(vol2)
## asf(vol2)
## print(asf(vol2), Inf)


###################################################
### code chunk number 23: vol1 (eval = FALSE)
###################################################
## csf(vol2, n.init = 2000)
## csf(vol2, n.init = 100)


###################################################
### code chunk number 24: vol1 (eval = FALSE)
###################################################
## csf(vol2, verbose = TRUE)


###################################################
### code chunk number 25: inus1
###################################################
dat.inu.1 <- allCombs(c(2, 2, 2)) -1
dat.inu.2 <- some(dat.inu.1, 40, replace = TRUE)
dat.inu.3 <- selectCases1("A + a*B <-> C",  con = 1, cov = 1, dat.inu.2)
asf(cna(dat.inu.3, con = 1, cov = 1, inus.only = FALSE))


###################################################
### code chunk number 26: inus2
###################################################
set.seed(26)
dat.inu.4 <- some(dat.inu.1, 40, replace = TRUE)
dat.inu.5 <- selectCases1("A + a*B <-> C", con = .8, cov = .8, dat.inu.4)
asf(cna(dat.inu.5, con = .8, cov = .8, inus.only = FALSE))


###################################################
### code chunk number 27: inus4
###################################################
asf(cna(dat.inu.5, con = .8, cov = .8, inus.only = TRUE))


###################################################
### code chunk number 28: d.edu1
###################################################
printCols <- c("condition", "consistency", "coverage", "exhaustiveness")
csf(cna(d.educate, details = "exhaust"))[printCols]


###################################################
### code chunk number 29: d.edu2
###################################################
csf(cna(d.educate[-1,], details = "exhaust"))[printCols]


###################################################
### code chunk number 30: d.edu3
###################################################
printCols <- c("condition", "consistency", "coverage", "faithfulness")
csf(cna(d.educate, details = "faithful"))[printCols]


###################################################
### code chunk number 31: d.edu4
###################################################
csf(cna(rbind(d.educate,c(1,1,0,1,0)), con = .8, details = "f"))[printCols]


###################################################
### code chunk number 32: rownames
###################################################
rownames(d.educate) <- 1:8


###################################################
### code chunk number 33: coherence
###################################################
d.edu.exp1 <- rbind(d.educate, c(1,0,1,0,0))
printCols <- c("condition", "consistency", "coverage", "coherence")
csf(cna(d.edu.exp1, con = .8, details = "cohere"))[printCols]


###################################################
### code chunk number 34: redundant1
###################################################
(dat.redun <- ct2df(selectCases("(A*B + C <-> D)*(a + c <-> E)")))


###################################################
### code chunk number 35: redundant2
###################################################
printCols <- c("condition", "consistency", "coverage", "inus", "redundant")
csf(cna(dat.redun, details = "r"))[printCols]


###################################################
### code chunk number 36: redundant3
###################################################
csf(cna(dat.redun, details = "r"), inus.only = FALSE, 
    minimalizeCsf = FALSE)[printCols]


###################################################
### code chunk number 37: redundant4a
###################################################
options(width=80)


###################################################
### code chunk number 38: redundant4
###################################################
printCols <- c("condition", "consistency", "coverage", "inus")
csf(fscna(d.autonomy, ordering = list("AU"), con = .9, cov = .94, 
   maxstep = c(2, 2, 8), inus.only = FALSE))[printCols]


###################################################
### code chunk number 39: redundant5
###################################################
condTbl("EM -> SP", fsct(d.autonomy))


###################################################
### code chunk number 40: redundant6
###################################################
csf(fscna(d.autonomy, ordering = list("AU"), con = .9, cov = .94, 
   maxstep = c(2, 2, 8), inus.only = TRUE))[printCols]


###################################################
### code chunk number 41: cycle1
###################################################
printCols <- c("condition", "cyclic")
csf(fscna(d.pacts, con = .8, cov = .8, acyclic.only = FALSE, 
   details = "cy"))[printCols]
csf(fscna(d.pacts, con = .8, cov = .8, acyclic.only = TRUE, 
   details = "cy"))[printCols]


###################################################
### code chunk number 42: ambigu1
###################################################
dat7 <- selectCases("a*B + A*b + B*C <-> D")
printCols <- c("condition", "consistency", "coverage", "inus",
  "exhaustiveness")
csf(cna(dat7, details = c("exhaust", "inus")))[printCols]


###################################################
### code chunk number 43: ambigu2
###################################################
csf(mvcna(d.pban, cov = .95, maxstep = c(3, 5, 10)))["condition"]


###################################################
### code chunk number 44: ambigu3
###################################################
ana.pact.1 <- fscna(d.pacts, ordering = list("PACT"), con = .8, cov = .8,
  maxstep = c(4, 5, 15), details = TRUE)
csf.pact.1 <- csf(ana.pact.1, Inf)
length(csf.pact.1$condition)


###################################################
### code chunk number 45: ambigu4
###################################################
csf.pact.1.ex <- subset(csf.pact.1, exhaustiveness >= .85)
length(csf.pact.1.ex$condition)


###################################################
### code chunk number 46: ambigu5
###################################################
csf.pact.1.ex.co <- subset(csf.pact.1.ex, 
  coherence == max(csf.pact.1.ex$coherence))
length(csf.pact.1.ex.co$condition)


###################################################
### code chunk number 47: back1
###################################################
dat.aut.2 <- d.autonomy[15:30, c("AU","EM","SP","CO","RE","DE")]
ana.aut.3 <- fscna(dat.aut.2, con = .91, cov = .91,
  ordering = list(c("RE", "DE","SP","CO"),"EM","AU"),
  strict = TRUE)
fscond(csf(ana.aut.3)$condition, dat.aut.2)


###################################################
### code chunk number 48: back2
###################################################
group.by.outcome(fscond(asf(ana.aut.3)$condition, dat.aut.2))$AU


###################################################
### code chunk number 49: details (eval = FALSE)
###################################################
## # Draw a ground truth.
## fullData <- mvct(allCombs(c(4,4,4,4,4)))
## groundTruth <- randomCsf(fullData, n.asf = 2, compl = 2)
## # Generate ideal data for groundTruth.
## idealData <- ct2df(selectCases(groundTruth, fullData))
## # Introduce 20% fragmentation.
## fragData <- idealData[-sample(1:nrow(idealData), nrow(idealData)*0.2), ] 
## # Introduce 10% random noise (cases incompatible with ground truth).
## incompCases <- dplyr::setdiff(ct2df(fullData), idealData)
## x <- rbind(ct2df(fullData[sample(1:nrow(fullData), 
##    nrow(incompCases) * 0.1), ]), fragData)  
## # Run CNA without an ordering.
## csfs <- csf(mvcna(x, con = .7, cov = .7, maxstep = c(3, 3, 12)))
## # Check whether no causal fallacy (no false positive) is returned.
## if(length(csfs$condition)==0) {
##   TRUE } else {any(is.submodel(csfs$condition, groundTruth))}


