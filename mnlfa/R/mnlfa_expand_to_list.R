## File Name: mnlfa_expand_to_list.R
## File Version: 0.03

mnlfa_expand_to_list <- function(x, names_list)
{
    I <- length(names_list)
    res <- mnlfa_convert_to_list(x=x, I=I, names_list=names_list,
                create_list=TRUE)
    return(res)
}
