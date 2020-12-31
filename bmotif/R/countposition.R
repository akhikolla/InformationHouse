countposition <- function(M, p, jZ, jP, JP, JZ, JP3 = NULL, JZ3 = NULL, KP3 = NULL, KZ3 = NULL, MT, N, NT, dZ, dP, Z, Y, X, P, Q, R, MTA = NULL, MTB = NULL, MTC = NULL, MTD = NULL, MA = NULL, MB = NULL, MC = NULL, MD = NULL, NTA = NULL, NTB = NULL, NTC = NULL, Na = NULL, NB = NULL, NC = NULL){
  if (p == 1) {a <- dP}
  else if (p == 2) {a <- dZ}
  else if (p == 3) {a <- P %*% jP - dP}
  else if (p == 4) {a <- dZ * (dZ - jZ) / 2}
  else if (p == 5) {a <- dP * (dP - jP) / 2}
  else if (p == 6) {a <- Z %*% jZ - dZ}
  else if (p == 7) {a <- dP * (dP - jP) * (dP - 2 * jP) / 6}
  else if (p == 8) {a <- M %*% ((dP - jP) * (dP - 2 * jP)) / 2}
  else if (p == 9) {a <- (P * R) %*% jP}
  else if (p == 10) {a <- (P * Q) %*% jP}
  else if (p == 11) {a <- (X * Z) %*% jZ}
  else if (p == 12) {a <- (Y * Z) %*% jZ}
  else if (p == 13) {a <- (P * (P - JP)) %*% jP / 2 - dP * (dP - jP) / 2}
  else if (p == 14) {a <- (Z * (Z - JZ)) %*% jZ / 2 - dZ * (dZ - jZ) / 2}
  else if (p == 15) {a <- MT %*% ((dZ - jZ) * (dZ - 2 * jZ)) / 2}
  else if (p == 16) {a <- dZ * (dZ - jZ) * (dZ - 2 * jZ) / 6}
  else if (p == 17) {a <- dP * (dP - jP) * (dP - 2 * jP) * (dP - 3 * jP) / 24}
  else if (p == 18) {a <- M %*% ((dP - jP) * (dP - 2 * jP) * (dP - 3 * jP)) / 6}
  else if (p == 19) {a <- (P * R * (R - JP)) %*% jP / 2}
  else if (p == 20) {a <- (P * Q * (Q - JP)) %*% jP / 2}
  else if (p == 21) {a <- (N * (M %*% ((Q - JP) * P))) %*% jP}
  else if (p == 22) {a <- (M * (M %*% (Q * (Q - JP)))) %*% jP / 2}
  else if (p == 23) {a <- (P * Q * R) %*% jP}
  else if (p == 24) {a <- (N * (M %*% (P * R))) %*% jP}
  else if (p == 25) {a <- (M * (M %*% (Q * R))) %*% jP / 2}
  else if (p == 26) {a <- (P * (P - JP) * R) %*% jP / 2}
  else if (p == 27) {a <- (P * (P - JP) * Q) %*% jP / 2}
  else if (p == 28) {a <- (N * (M %*% (P * (P - JP)))) %*% jP / 2}
  else if (p == 29) {a <- (M * (M %*% (Q * (P - JP)))) %*% jP}
  else if (p == 30) {a <- (P * (P - JP) * (P - 2 * JP)) %*% jP / 6 - dP * (dP - jP) * (dP - 2 * jP) / 6}
  else if (p == 31) {a <- (M * M %*% ((P - JP) * (P - 2 * JP))) %*% jP / 4 - M %*% ((dP - jP) * (dP - 2 * jP)) / 4}
  else if (p == 32) {a <- (NT * (MT %*% ((Y - JZ) * Z))) %*% jZ}
  else if (p == 33) {a <- (MT * (MT %*% (Y * (Y - JZ)))) %*% jZ / 2}
  else if (p == 34) {a <- (Z * X * (X - JZ)) %*% jZ / 2}
  else if (p == 35) {a <- (Z * Y * (Y - JZ)) %*% jZ / 2}
  else if (p == 36) {a <- (NT * (MT %*% (Z*X))) %*% jZ}
  else if (p == 37) {a <- (MT * (MT %*% (Y * X))) %*% jZ / 2}
  else if (p == 38) {a <- (Z * Y * X) %*% jZ}
  else if (p == 39) {a <- (NT * (MT %*% (Z * (Z - JZ)))) %*% jZ / 2}
  else if (p == 40) {a <- (MT * (MT %*% (Y * (Z - JZ)))) %*% jZ}
  else if (p == 41) {a <- (Z * (Z - JZ) * X) %*% jZ / 2}
  else if (p == 42) {a <- (Z * (Z - JZ) * Y) %*% jZ / 2}
  else if (p == 43) {a <- (MT * (MT %*% ((Z - JZ) * (Z - 2 * JZ)))) %*% jZ / 4 - MT %*% ((dZ - jZ) * (dZ - 2 * jZ)) / 4}
  else if (p == 44) {a <- (Z * (Z - JZ) * (Z - 2 * JZ)) %*% jZ / 6 - dZ * (dZ - jZ) * (dZ - 2 * jZ) / 6}
  else if (p == 45) {a <- MT %*% ((dZ - jZ) * (dZ - 2 * jZ) * (dZ - 3 * jZ)) / 6}
  else if (p == 46) {a <- dZ * (dZ - jZ) * (dZ - 2 * jZ) * (dZ - 3 * jZ) / 24}

  # 6 node
  else if (p == 47) {a <- dP * (dP - jP) * (dP - 2 * jP) * (dP - 3 * jP) * (dP - 4 * jP) / 120}
  else if (p == 48) {a <- M %*% ((dP - jP) * (dP - 2 * jP) * (dP - 3 * jP) * (dP - 4 * jP)) / 24}
  else if (p == 49) {a <- (P * R * (R - JP) * (R - 2 * JP)) %*% jP / 6}
  else if (p == 50) {a <- (P * Q * (Q - JP) * (Q - 2 * JP)) %*% jP / 6}
  else if (p == 51) {a <- (N * (M %*% (P * (Q - JP) * (Q - 2 * JP)))) %*% jP / 2}
  else if (p == 52) {a <- (M * (M %*% (Q * (Q - JP) * (Q - 2 * JP)))) %*% jP / 6}
  else if (p == 53) {a <- (P * Q * R * (R - JP)) %*% jP / 2}
  else if (p == 54) {a <- (P * Q * (Q - JP) * R) %*% jP / 2}
  else if (p == 55) {a <- (M * (N %*% (P * Q * (Q - JP)))) %*% jP / 2}
  else if (p == 56) {a <- (M * (N %*% (P * Q * (R - JP)))) %*% jP}
  else if (p == 57) {a <- (M * (M %*% (Q * (Q - JP) * R))) %*% jP / 2}
  else if (p == 58) {a <- (P * (P - JP) * R * (R - JP)) %*% jP / 4}
  else if (p == 59) {a <- (P * (P - JP) * Q * (Q - JP)) %*% jP / 4}
  else if (p == 60) {a <- (N * (M %*% ((Q - JP) * P * (P - JP)))) %*% jP / 2}
  else if (p == 61) {a <- (M * (M %*% ((P - JP) * Q * (Q - JP)))) %*% jP / 2}
  else if (p == 62) {a <- (P * (P - JP) * Q * R) %*% jP / 2}
  else if (p == 63) {a <- (N * (M %*% (P * (P - JP) * R))) %*% jP / 2}
  else if (p == 64) {a <- (M * (M %*% ((P - JP) * Q * R))) %*% jP / 2}
  else if (p == 65) {a <- (P * (P - JP) * (P - 2 * JP) * R) %*% jP / 6}
  else if (p == 66) {a <- (P * (P - JP) * (P - 2 * JP) * Q) %*% jP / 6}
  else if (p == 67) {a <- (N * (M %*% (P * (P - JP) * (P - 2 * JP)))) %*% jP / 6}
  else if (p == 68) {a <- (M * (M %*% ((P - JP) * (P - 2 * JP) * Q))) %*% jP / 2}
  else if (p == 69) {a <- (P * (P - JP) * (P - 2 * JP) * (P - 3 * JP)) %*% jP / 24 - dP * (dP - jP) * (dP - 2 * jP) * (dP - 3 * jP) / 24}
  else if (p == 70) {a <- (M * (M %*% ((P - JP) * (P - 2 * JP) * (P - 3 * JP)))) %*% jP / 12 - M %*% ((dP - jP) * (dP - 2 * jP) * (dP - 3 * jP)) / 12}
  else if (p == 71) {a <- apply((NTB * (NTB - JP3) * MTA * KP3), MARGIN = 1, sum) / 2}
  else if (p == 72) {a <- apply((MTD * (MTD - JP3) * MTA * KP3), MARGIN = 1, sum) / 4}
  else if (p == 73) {a <- apply((NB * (NB - JZ3) * MA * KZ3), MARGIN = 1, sum) / 2}
  else if (p == 74) {a <- apply((MD * (MD - JZ3) * MA * KZ3), MARGIN = 1, sum) / 4}
  else if (p == 75) {a <- apply(MTD * MTA * NTB * KP3, MARGIN = 1, sum)}
  else if (p == 76) {a <- apply(MTA * NTB * NTC * KP3, MARGIN = 1, sum) / 2}
  else if (p == 77) {a <- apply(MD * MB * MC * KZ3, MARGIN = 1, sum) / 2}
  else if (p == 78) {a <- apply(MB * NB * Na * KZ3, MARGIN = 1, sum)}
  else if (p == 79) {a <- apply(MTB * MTC * MTD * KP3, MARGIN = 1, sum) / 2}
  else if (p == 80) {a <- apply(MTB * NTB * NTA * KP3, MARGIN = 1, sum)}
  else if (p == 81) {a <- apply(MA * MD * NB * KZ3, MARGIN = 1, sum)}
  else if (p == 82) {a <- apply(MA * NB * NC * KZ3, MARGIN = 1, sum) / 2}
  else if (p == 83) {a <- apply(MTB * NTA * NTC * KP3, MARGIN = 1, sum)}
  else if (p == 84) {a <- apply(NTA * MTB * MTD * KP3, MARGIN = 1, sum)}
  else if (p == 85) {a <- apply(MTB * MTC * NTC * KP3, MARGIN = 1, sum)}
  else if (p == 86) {a <- apply(MB * Na * NC * KZ3, MARGIN = 1, sum)}
  else if (p == 87) {a <- apply(MD * MB * Na * KZ3, MARGIN = 1, sum)}
  else if (p == 88) {a <- apply(NB * MB * MC * KZ3, MARGIN = 1, sum)}
  else if (p == 89) {a <- apply(MTA * NTA * NTB * KP3, MARGIN = 1, sum)}
  else if (p == 90) {a <- apply(MTA * MTB * NTB * KP3, MARGIN = 1, sum)}
  else if (p == 91) {a <- apply(MTA * MTB * MTD * KP3, MARGIN = 1, sum)}
  else if (p == 92) {a <- apply(MA * Na * NB * KZ3, MARGIN = 1, sum)}
  else if (p == 93) {a <- apply(MB * NB * MA * KZ3, MARGIN = 1, sum)}
  else if (p == 94) {a <- apply(MA * MB * MD * KZ3, MARGIN = 1, sum)}
  else if (p == 95) {a <- apply(MTB * NTA * (NTA - JP3) * KP3, MARGIN = 1, sum) / 2}
  else if (p == 96) {a <- apply(MTB * (MTB - JP3) * NTA * KP3, MARGIN = 1, sum) / 2}
  else if (p == 97) {a <- apply(MTB * MTC * (MTC - JP3) * KP3, MARGIN = 1, sum) / 2}
  else if (p == 98) {a <- apply(MD * MA * Na * KZ3, MARGIN = 1, sum) / 2}
  else if (p == 99) {a <- apply(MA * NB * MC * KZ3, MARGIN = 1, sum)}
  else if (p == 100) {a <- apply(MTA * NTA * (NTA - JP3) * KP3, MARGIN = 1, sum) / 4}
  else if (p == 101) {a <- apply(MTA * MTB * (MTB - JP3) * KP3, MARGIN = 1, sum) / 2}
  else if (p == 102) {a <- apply(MA * (MA - JZ3) * NB * KZ3, MARGIN = 1, sum) / 2}
  else if (p == 103) {a <- apply(MD * MA * (MA - JZ3) * KZ3, MARGIN = 1, sum) / 4}
  else if (p == 104) {a <- apply(MTD * MTA * NTA * KP3, MARGIN = 1, sum) / 2}
  else if (p == 105) {a <- apply(MTB * MTA * NTC * KP3, MARGIN = 1, sum)}
  else if (p == 106) {a <- apply(MB * Na * (Na - JZ3) * KZ3, MARGIN = 1, sum) / 2}
  else if (p == 107) {a <- apply(MB * (MB - JZ3) * Na * KZ3, MARGIN = 1, sum) / 2}
  else if (p == 108) {a <- apply(MB * (MB - JZ3) * MC * KZ3, MARGIN = 1, sum) / 2}
  else if (p == 109) {a <- apply(MTA * (MTA - JP3) * NTB * KP3, MARGIN = 1, sum) / 2}
  else if (p == 110) {a <- apply(MTA * (MTA - JP3) * MTD * KP3, MARGIN = 1, sum) / 4}
  else if (p == 111) {a <- apply(MA * Na * (Na - JZ3) * KZ3, MARGIN = 1, sum) / 4}
  else if (p == 112) {a <- apply(MA * MB * (MB - JZ3) * KZ3, MARGIN = 1, sum) / 2}
  else if (p == 113) {a <- apply(MTB * MTC * NTA * KP3, MARGIN = 1, sum) / 2}
  else if (p == 114) {a <- apply(MB * MC * Na * KZ3, MARGIN = 1, sum) / 2}
  else if (p == 115) {a <- apply(MTA * MTB * NTA * KP3, MARGIN = 1, sum)}
  else if (p == 116) {a <- apply(MTA * MTB * MTC * KP3, MARGIN = 1, sum) / 2}
  else if (p == 117) {a <- apply(MA * Na * MB * KZ3, MARGIN = 1, sum)}
  else if (p == 118) {a <- apply(MA * MB * MC * KZ3, MARGIN = 1, sum) / 2}
  else if (p == 119) {a <- apply(MTA * (MTA - JP3) * NTA * KP3, MARGIN = 1, sum) / 4}
  else if (p == 120) {a <- apply(MTA * (MTA - JP3) * MTB * KP3, MARGIN = 1, sum) / 2}
  else if (p == 121) {a <- apply(MA * (MA - JZ3) * Na * KZ3, MARGIN = 1, sum) / 4}
  else if (p == 122) {a <- apply(MB * MA * (MA - JZ3) * KZ3, MARGIN = 1, sum) / 2}
  else if (p == 123) {a <- apply(MTA * (MTA - JP3) * (MTA - 2 * JP3) * KP3, MARGIN = 1, sum) / 12}
  else if (p == 124) {a <- apply(MA * (MA - JZ3) * (MA - 2 * JZ3) * KZ3, MARGIN = 1, sum) / 12}
  else if (p == 125) {a <- (NT * (MT %*% (Z * (Y - JZ) * (Y - 2 * JZ)))) %*% jZ / 2}
  else if (p == 126) {a <- (MT * (MT %*% (Y * (Y - JZ) * (Y - 2 * JZ)))) %*% jZ / 6}
  else if (p == 127) {a <- (Z * X * (X - JZ) * (X - 2 * JZ)) %*% jZ / 6}
  else if (p == 128) {a <- (Z * Y * (Y - JZ) * (Y - 2 * JZ)) %*% jZ / 6}
  else if (p == 129) {a <- (MT * (NT %*% (Z * Y * (Y - JZ)))) %*% jZ / 2}
  else if (p == 130) {a <- (MT * (NT %*% (Z * Y * (X - JZ)))) %*% jZ}
  else if (p == 131) {a <- (MT * (MT %*% (Y * (Y - JZ) * X))) %*% jZ / 2}
  else if (p == 132) {a <- (Z * Y * X * (X - JZ)) %*% jZ / 2}
  else if (p == 133) {a <- (Z * Y * (Y - JZ) * X) %*% jZ / 2}
  else if (p == 134) {a <- (NT * (MT %*% ((Y - JZ) * Z * (Z - JZ)))) %*% jZ / 2}
  else if (p == 135) {a <- (MT * (MT %*% ((Z - JZ) * Y * (Y - JZ)))) %*% jZ / 2}
  else if (p == 136) {a <- (Z * (Z - JZ) * X * (X - JZ)) %*% jZ / 4}
  else if (p == 137) {a <- (Z * (Z - JZ) * Y * (Y - JZ)) %*% jZ / 4}
  else if (p == 138) {a <- (NT * (MT %*% (Z * (Z - JZ) * X))) %*% jZ / 2}
  else if (p == 139) {a <- (MT * (MT %*% ((Z - JZ) * Y * X))) %*% jZ / 2}
  else if (p == 140) {a <- (Z * (Z - JZ) * Y * X) %*% jZ / 2}
  else if (p == 141) {a <- (NT * (MT %*% (Z * (Z - JZ) * (Z - 2 * JZ)))) %*% jZ / 6}
  else if (p == 142) {a <- (MT * (MT %*% ((Z - JZ) * (Z - 2 * JZ) * Y))) %*% jZ / 2}
  else if (p == 143) {a <- (Z * (Z - JZ) * (Z - 2 * JZ) * X) %*% jZ / 6}
  else if (p == 144) {a <- (Z * (Z - JZ) * (Z - 2 * JZ) * Y) %*% jZ / 6}
  else if (p == 145) {a <- (MT * (MT %*% ((Z - JZ) * (Z - 2 * JZ) * (Z - 3 * JZ)))) %*% jZ / 12 - MT %*% ((dZ - jZ) * (dZ - 2 * jZ) * (dZ - 3 * jZ)) / 12}
  else if (p == 146) {a <- (Z * (Z - JZ) * (Z - 2 * JZ) * (Z - 3 * JZ)) %*% jZ / 24 - dZ * (dZ - jZ) * (dZ - 2 * jZ) * (dZ - 3 * jZ) / 24}
  else if (p == 147) {a <- MT %*% ((dZ - jZ) * (dZ - 2 * jZ) * (dZ - 3 * jZ) * (dZ - 4 * jZ)) / 24}
  else if (p == 148) {a <- dZ * (dZ - jZ) * (dZ - 2 * jZ) * (dZ - 3 * jZ) * (dZ - 4 * jZ) / 120}

  else {a <- -1}
  matrix(a)
}
