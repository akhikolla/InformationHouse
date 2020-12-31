ListAsEsomNeurons <- function(WtsList, Lines, Columns){
# NeuronWeights <- ListAsEsomNeurons(WtsList, Lines, Columns)
# INPUT
# WtsList[ind,Weight]	Matrix with weights in the 2nd dimension(not list() like in R)
# Lines                        Lines/Height of the desired grid
# Columns		       Columns/Width of the desired grid
# OUTPUT
# EsomNeurons           3dimensional array containing the weights of the grid


    nr_of_arguments = ncol(WtsList)

    # first row will be filled with nans as placeholder
    newMatrix <- array(NaN, dim=c(Lines,Columns,nr_of_arguments));

    # list is filled linewise
    for(xline in 1:Lines){
        for(xcol in 1:Columns){
            newMatrix[xline,xcol,] = WtsList[(xline-1)*Columns + xcol,]
        }
    }

    newMatrix
}

