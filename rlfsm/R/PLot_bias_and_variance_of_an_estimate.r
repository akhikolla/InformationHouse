
#### PLot bias and variance ####

#' A function to plot variance- and bias dependencies of estimators on the lengths of sample paths. Works in conjunction with \code{\link{MCestimLFSM}} function.
#' @param data a list created by \code{\link{MCestimLFSM}}
#' @return The function returns a ggplot2 graph.
#' @seealso \code{\link{Plot_dens}}
#' @examples
#'
#' # Light weight computaions
#'
#' m<-25; M<-50
#' alpha<-1.8; H<-0.8; sigma<-0.3
#' S<-c(1:3)*1e2
#' p<-.4; p_prime<-.2; t1<-1; t2<-2
#' k<-2; NmonteC<-50
#'
#' # Here is the continuous H-1/alpha inference procedure
#' theor_3_1_H_clt<-MCestimLFSM(s=S,fr='H',Nmc=NmonteC,
#'                      m=m,M=M,alpha=alpha,H=H,
#'                      sigma=sigma,ContinEstim,
#'                      t1=t1,t2=t2,p=p,k=k)
#' Plot_vb(theor_3_1_H_clt)
#'
#' \donttest{
#' # More demanding example (it is better to use multicore setup)
#' # General low frequency inference
#'
#' m<-45; M<-60
#' alpha<-0.8; H<-0.8; sigma<-0.3
#' S<-c(1:15)*1e2
#' p<-.4; t1<-1; t2<-2
#' NmonteC<-50
#'
#' # Here is the continuous H-1/alpha inference procedure
#' theor_4_1_H_clt<-MCestimLFSM(s=S,fr='H',Nmc=NmonteC,
#'                      m=m,M=M,alpha=alpha,H=H,
#'                      sigma=sigma,GenLowEstim,
#'                      t1=t1,t2=t2,p=p)
#' Plot_vb(theor_4_1_H_clt)
#' }
#'
#' @export
Plot_vb<-function(data){

    s<-NULL # avoids NOTEs when being builded
    value<-NULL # avoids NOTEs when being builded

    bs <- data$biases
    ms <- data$means
    sds <- data$sds

    part_bs <- reshape2::melt(bs, id.vars="s")
    part_bs$type='bias'

    part_ms <- reshape2::melt(ms, id.vars="s")
    part_ms$type='mean'

    part_sds <- reshape2::melt(sds, id.vars="s")
    part_sds$type='sd'

    data_to_plot <- rbind(part_bs, part_ms, part_sds)

    pl <- ggplot2::ggplot(data_to_plot, aes(x=s, y=value)) + ggplot2::geom_point(colour="cyan4") #geom_line(size = 0.5, colour = "blue") #
    pl<-pl + ggplot2::facet_wrap(variable~type, scales = "free", labeller = label_both) #facet_grid(alpha ~ H, scales = "free", labeller = label_both)
    pl<-pl + ggplot2::theme_bw()
    pl<-pl + ggplot2::geom_hline(yintercept=0, colour="brown")
    pl<-pl + ggplot2::geom_smooth() # local regression

    pl
}




# A special version for MP paper (black and white palette)
Plot_vb_paper<-function(data){

    n<-NULL # avoids NOTEs when being built
    value<-NULL # avoids NOTEs when being built

    data<-as.data.frame(data)
    data_to_plot <- melt(data, id.vars="n")

    pl <- ggplot2::ggplot(data_to_plot, aes(x=n, y=value)) + ggplot2::geom_point(colour="black") #geom_line(size = 0.5, colour = "blue") #
    pl<-pl + ggplot2::facet_wrap(~variable, scales = "free", labeller = label_both) #facet_grid(alpha ~ H, scales = "free", labeller = label_both)
    pl<-pl + ggplot2::theme_bw()
    pl<-pl + ggplot2::geom_hline(yintercept=0, colour="gray40")
    pl<-pl+ggplot2::geom_smooth() # local regression

    ggplot2::ggsave("Variance and bias dependence on n.pdf", width = 10, height = 10)
    pl
}


