## File Name: immer_score_person_adjusted.R
## File Version: 0.04


immer_score_person_adjusted <- function( sum_score, max_pers, eps)
{
    score_pers <- sum_score
    score_pers <- ifelse( sum_score==0, eps, score_pers )
    score_pers <- ifelse( sum_score==max_pers, max_pers - eps, score_pers )
    return(score_pers)
}
