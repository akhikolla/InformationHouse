
#' Export results into a template file
#'
#' @param numbers    vector of results
#' @param template.f template file name
#' @param out.f      output file name
#' @param sub.str    pattern of string to be replaced
#'
#'
#' @export
#'
tkExpRst <- function(numbers, template.f,  out.f="rst.txt", sub.str="AA") {
    if (!file.exists(template.f)) {
        return(NULL)
    }
    ##read
    tpla <- readChar(template.f, file.info(template.f)$size);

    ##substitute
    for (i in 1:length(numbers)) {
        tpla <- sub(sub.str, numbers[i], tpla);
    }

    ##write out
    write(tpla, file = out.f);
}


#' Import objects in a list into a designated environment
#'
#' @param alist list of objects
#' @param dest.env designated environment
#'
#' @export
#'
tkMakeLocal <- function(alist, dest.env) {
    for (i in 1:length(alist)) {
        assign(names(alist[i]), alist[[i]], dest.env);
    }
}


#' Call function by its name organized as a vector
#'
#' @param vec function names as a vector
#' @param ... Parameters needed for the actual function
#'
#' @export
#'
tkCallFun <- function(vec, ...) {
    rst <- NULL
    eval(parse(text = paste("rst <- ",
                            paste(vec, collapse = ""),
                            "(...)",
                            sep = "")
               )
         )
    rst
}
