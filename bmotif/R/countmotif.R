countmotif <- function(x, motif, z, p, JP, JZ, P, Q, R, Z, Y, X, dP, jP, dZ, jZ, J3 = NULL, MA = NULL, MB = NULL, MC = NULL, MD = NULL, Na = NULL, NB = NULL, NC = NULL){
  if(motif == 1){
    sum(x)
  } else if(motif == 2){
    sum(dZ * (dZ - jZ)) / 2
  } else if(motif == 3){
    sum(dP * (dP - jP)) / 2
  } else if(motif == 4){
    sum(dP * (dP - jP) * (dP - 2 * jP)) / 6
  } else if(motif == 5){
    sum(Z * Y)
  } else if(motif == 6){
    sum(Z * (Z - JZ)) / 4 - sum(dZ * (dZ - jZ)) / 4
  } else if(motif == 7){
    sum(dZ * (dZ - jZ) * (dZ - 2 * jZ)) / 6
  } else if(motif == 8){
    sum(dP * (dP - jP) * (dP - 2 * jP) * (dP - 3 * jP)) / 24
  } else if(motif == 9){
    sum(P * Q * (Q - JP)) / 2
  } else if(motif == 10){
    sum(P * Q * R) / 2
  } else if(motif == 11){
    sum(P * (P - JP) * Q) / 2
  } else if(motif == 12){
    sum(P * (P - JP) * (P - 2 * JP)) / 12 - sum(dP * (dP - jP) * (dP - 2 * jP)) / 12
  } else if(motif == 13){
    sum(Z * Y * (Y - JZ)) / 2
  } else if(motif == 14){
    sum(Z * Y * X) / 2
  } else if(motif == 15){
    sum(Z * (Z - JZ) * Y) / 2
  } else if(motif == 16){
    sum(Z * (Z - JZ) * (Z - 2 * JZ)) / 12 - sum(dZ * (dZ - jZ) * (dZ - 2 * jZ)) / 12
  } else if(motif == 17){
    sum(dZ * (dZ - jZ) * (dZ - 2 * jZ) * (dZ - 3 * jZ)) / 24
  } else if(motif == 18){
    sum(dP * (dP - jP) * (dP - 2 * jP) * (dP - 3 * jP) * (dP - 4 * jP)) / 120
  } else if(motif == 19){
    sum(P * Q * (Q - JP) * (Q - 2 * JP)) / 6
  } else if(motif == 20){
    sum(Q * P * R * (Q - JP)) / 2
  } else if(motif == 21){
    sum(Q * (Q - JP) * P * (P - JP)) / 4
  } else if(motif == 22){
    sum(Q * P * (P - JP) * R) / 4
  } else if(motif == 23){
    sum(P * (P - JP) * (P - 2 * JP) * Q) / 6
  } else if(motif == 24){
    sum(P * (P - JP) * (P - 2 * JP) * (P - 3 * JP)) / 48 - sum(dP * (dP - jP) * (dP - 2 * jP) * (dP - 3 * jP)) / 48
  } else if(motif == 25){
    sum(MA * NB * (NB - J3)) / 4
  } else if(motif == 26){
    if(z <= p){
      sum(MD * MB * MC) / 2
    } else {
      sum(NB * NC * MA) / 2
    }
  } else if(motif == 27){
    if(z <= p){
      sum(NB * NC * MA) / 2
    } else {
      sum(MD * MB * MC) / 2
    }
  } else if(motif == 28){
    sum(MB * MC * NC)
  } else if(motif == 29){
    sum(MA * MB * MD)
  } else if(motif == 30){
    if(z <= p){
      sum(MA * MB * NC) / 2
    } else {
      sum(MB * (MB - J3) * MC) / 2
    }
  } else if(motif == 31){
    if(z <= p){
      sum(MA * NB * (MA - J3)) / 4
    } else {
      sum(MB * (MB - J3) * MA) / 4
    }
  } else if(motif == 32){
    if(z <= p){
      sum(MB * (MB - J3) * MC) / 2
    } else {
      sum(MA * MB * NC) / 2
    }
  } else if(motif == 33){
    if(z <= p){
      sum(MB * (MB - J3) * MA) / 4
    } else {
      sum(MA * NB * (MA - J3)) / 4
    }
  } else if(motif == 34){
    sum(Na * MB * MC) / 6
  } else if(motif == 35){
    sum(MA * MB * MC) / 2
  } else if(motif == 36){
    sum(MA * (MA - J3) * MB) / 4
  } else if(motif == 37){
    sum(MA * (MA - J3) * (MA - 2 * J3)) / 36
  } else if(motif == 38){
    sum(Z * Y * (Y - JZ) * (Y - 2 * JZ)) / 6
  } else if(motif == 39){
    sum(Z * (Y - JZ) * X * Y) / 2
  } else if(motif == 40){
    sum(Z * (Z - JZ) * Y * (Y - JZ)) / 4
  } else if(motif == 41){
    sum(Z * (Z - JZ) * X * Y) / 4
  } else if(motif == 42){
    sum(Z * (Z - JZ) * (Z - 2 * JZ) * Y) / 6
  } else if(motif == 43){
    sum(Z * (Z - JZ) * (Z - 2 * JZ) * (Z - 3 * JZ)) / 48 - sum(dZ * (dZ - jZ) * (dZ - 2 * jZ) * (dZ - 3 * jZ)) / 48
  } else if(motif == 44){
    sum(dZ * (dZ - jZ) * (dZ - 2 * jZ) * (dZ - 3 * jZ) * (dZ - 4 * jZ)) / 120
  }
}
