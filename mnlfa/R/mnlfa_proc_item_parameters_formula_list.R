## File Name: mnlfa_proc_item_parameters_formula_list.R
## File Version: 0.04


mnlfa_proc_item_parameters_formula_list <- function(formula_parm, items)
{
    I <- length(items)
    if ( ! is.list(formula_parm) ){
        formula_parm <- mnlfa_convert_to_list(x=formula_parm, I=I, names_list=items)
    }
    return(formula_parm)
}
