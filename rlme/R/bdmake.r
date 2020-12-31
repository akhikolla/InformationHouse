bdmake <- function(mat1, mat2) {
    r1 = length(mat1[, 1])
    r2 = length(mat2[, 1])
    r3 = r1 + r2
    res = matrix(rep(0, r3^2), ncol = r3)
    res[1:r1, 1:r1] = mat1
    res[(r1 + 1):r3, (r1 + 1):r3] = mat2
    return(res)
}
