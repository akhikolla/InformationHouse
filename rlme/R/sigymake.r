sigymake <- function (I, sec, mat, siga2, sigw2, sige2) {
    students.per.school = apply(mat, 1, sum)
    n = sum(students.per.school)
    
    sigma.matrix = matrix(0, nrow=n, ncol=n)
    sigma.matrix.i = matrix(0, nrow=n, ncol=n)
    siggma.matrix  = matrix(0, nrow=n, ncol=n)
    
    iflag = 0
    
    position = 1
    
    for(school.index in 1:I) {
      school.range = position:(position + students.per.school[school.index] - 1)
      
      sigma.matrix[school.range, school.range] = siga2 + diag(students.per.school[school.index]) * sige2
      
      for(section.index in 1:sec[school.index]) {
        students = mat[school.index, section.index]
        
        section.range = position:(position+students-1)
        
        sigma.matrix[section.range, section.range] = sigma.matrix[section.range, section.range] + sigw2
        
        position = position + students
      }
      
      eigs = eigen(sigma.matrix[school.range, school.range], symmetric = T)
      if (min(eigs$values) < 0) {
        iflag = 1
      }
      
      sigma.matrix.i[school.range, school.range] = eigs$vectors %*% diag(1/eigs$values^0.5) %*% t(eigs$vectors)
      siggma.matrix[school.range, school.range] = eigs$vectors %*% diag(eigs$values^0.5) %*% t(eigs$vectors)
    }
    
    return(list(
      sigy2 = sigma.matrix,
      sigy12i = sigma.matrix.i,
      siggma = siggma.matrix,
      iflag = iflag
    ))
  }