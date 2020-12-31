motif_info <- function(m, node = TRUE, link = TRUE) {
# m is the number of the motif
# the function returns which node and link positions are in that motif
# result is a list, first component: vector of node positions
  # second component: vector of link positions
# if node = FALSE, only second component of the list returned
  # likewise if link = FALSE, only first component returned

  if (m == 1) {
    l <- list(c(1,2), c(1))
  }
  if (m == 2) {
    l <- list(c(3,4), c(2))
  }
  if (m == 3) {
    l <- list(c(5,6), c(3))
  }
  if (m == 4) {
    l <- list(c(7,8), c(4))
  }
  if (m == 5) {
    l <- list(9:12, 5:7)
  }
  if (m == 6) {
    l <- list(13:14, 8)
  }
  if (m == 7) {
    l <- list(15:16, 9)
  }
  if (m == 8) {
    l <- list(17:18, 10)
  }
  if (m == 9) {
    l <- list(19:22, 11:13)
  }
  if (m == 10) {
    l <- list(23:25, 14:15)
  }
  if (m == 11) {
    l <- list(26:29, 16:18)
  }
  if (m == 12) {
    l <- list(30:31, 19)
  }
  if (m == 13) {
    l <- list(32:35, 20:22)
  }
  if (m == 14) {
    l <- list(36:38, 23:24)
  }
  if (m == 15) {
    l <- list(39:42, 25:27)
  }
  if (m == 16) {
    l <- list(43:44, 28)
  }
  if (m == 17) {
    l <- list(45:46, 29)
  }
  if (m == 18) {
    l <- list(47:48, 30)
  }
  if (m == 19) {
    l <- list(49:52, 31:33)
  }
  if (m == 20) {
    l <- list(53:57, 34:37)
  }
  if (m == 21) {
    l <- list(58:61, 38:40)
  }
  if (m == 22) {
    l <- list(62:64, 41:42)
  }
  if (m == 23) {
    l <- list(65:68, 43:45)
  }
  if (m == 24) {
    l <- list(69:70, 46)
  }
  if (m == 25) {
    l <- list(71:74, 47:49)
  }
  if (m == 26) {
    l <- list(75:78, 50:52)
  }
  if (m == 27) {
    l <- list(79:82, 53:55)
  }
  if (m == 28) {
    l <- list(83:88, 56:60)
  }
  if (m == 29) {
    l <- list(89:94, 61:66)
  }
  if (m == 30) {
    l <- list(95:99, 67:70)
  }
  if (m == 31) {
    l <- list(100:103, 71:73)
  }
  if (m == 32) {
    l <- list(104:108, 74:77)
  }
  if (m == 33) {
    l <- list(109:112, 78:80)
  }
  if (m == 34) {
    l <- list(113:114, 81)
  }
  if (m == 35) {
    l <- list(115:118, 82:85)
  }
  if (m == 36) {
    l <- list(119:122, 86:88)
  }
  if (m == 37) {
    l <- list(123:124, 89)
  }
  if (m == 38) {
    l <- list(125:128, 90:92)
  }
  if (m == 39) {
    l <- list(129:133, 93:96)
  }
  if (m == 40) {
    l <- list(134:137, 97:99)
  }
  if (m == 41) {
    l <- list(138:140, 100:101)
  }
  if (m == 42) {
    l <- list(141:144, 102:104)
  }
  if (m == 43) {
    l <- list(145:146, 105)
  }
  if (m == 44) {
    l <- list(147:148, 106)
  }

  if(node & link) {
    return(l)
  }
  else if (node & !link) {
    return(l[[1]])
  }
  else if (!node & link) {
    return(l[[2]])
  }
  else {
    print('error in motif_info: must have at least one true')
  }
}
