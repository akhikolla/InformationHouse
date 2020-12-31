EsomNeuronsAsList <- function(EsomNeurons){
# EsomNeurons_list <- EsomNeuronsAsList(EsomNeurons)
# INPUT
# EsomNeurons[Row,Col,Weight]         3 dimensional Array with weights in the 3rd dimension, like EsomNeurons DataIO
# OUTPUT
# EsomNeurons_list         List of Weights as matrix (not list() like in R)
# Author: Florian Lerch

    # first row will be filled with nans as placeholder
    newMatrix <- matrix(nrow=(dim(EsomNeurons)[2]*dim(EsomNeurons)[1]),ncol=dim(EsomNeurons)[3])

    i=1
    # list gets filled linewise
    for(xline in 1:nrow(EsomNeurons)){ # Height
        for(xcol in 1:ncol(EsomNeurons)){ # Width
            weight_vector = EsomNeurons[xline,xcol,]
						newMatrix[i,] = weight_vector
						i=i+1
        }
    }
    return(newMatrix)
}
