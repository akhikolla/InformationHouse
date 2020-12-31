

#' (alpha,H,sigma)- density plot
#'
#' Plots the densities of the parameters (alpha,H,sigma) estimated in Monte-Carlo experiment.
#' Works in conjunction with \code{\link{MCestimLFSM}} function.
#' @param par_vec vector of parameters which are to be plotted
#' @param MC_data a list created by \code{\link{MCestimLFSM}}
#' @param Nnorm number of point sampled from standard normal distribution
#' @seealso \code{\link{Plot_vb}}  to plot variance- and bias dependencies on n.
#' @examples
#' \donttest{
#' m<-45; M<-60
#'
#' p<-.4; p_prime<-.2
#' t1<-1; t2<-2; k<-2
#'
#' NmonteC<-5e2
#' S<-c(1e3,1e4)
#' alpha<-.8; H<-0.8; sigma<-0.3
#' theor_4_1_clt_new<-MCestimLFSM(s=S,fr='L',Nmc=NmonteC,
#'                        m=m,M=M,
#'                        alpha=alpha,H=H,sigma=sigma,
#'                        GenLowEstim,t1=t1,t2=t2,p=p)
#' l_plot<-Plot_dens(par_vec=c('sigma','alpha','H'), MC_data=theor_4_1_clt_new, Nnorm=1e7)
#' l_plot
#'
#' }
#' @export
Plot_dens<-function(par_vec=c('alpha','H','sigma'), MC_data, Nnorm=1e7) {

    N<-NULL # avoids NOTEs when being builded
    values<-NULL # avoids NOTEs when being builded
    ff<-NULL
    s<-NULL
    index<-NULL

    fr=MC_data$freq

    # Generate standard normal
    norm_dens<-rnorm(Nnorm)
    nor.data<-cbind(N="true_normal", as.data.frame(norm_dens))
    names(nor.data)<-c('s','values')

    #MC_dataset<-MC_data$MC_dataset
    MC_dataset<-MC_data$data_nor

    # Inscriptions on plots
    if(fr=='L') freq_text<-'Low' else freq_text<-'High'
    #if(MC_data$Inference==) freq_text<-'Low' else freq_text<-'High'

    # Plot a parameter
    ff<-function(param){

        #chart_name<-paste("Distribution of", param, "estimate.", "Theorem 3.1.", freq_text, "frequency", sep = " ")
        index  <-which(names(MC_dataset)==param)
        index_s<-which(names(MC_dataset)=='s')

        data_to_plot<-MC_dataset[, c(index_s,index)]
        names(data_to_plot)<-c('s','values')
        data_to_plot<-rbind(data_to_plot,nor.data)

        ggp <- ggplot(data_to_plot, aes(x=values)) +

            geom_line(aes(colour=s, alpha=s), stat="density", size=1)+
            scale_alpha_discrete(range = c(0.8,0.3)) +
            # ggtitle(expression(chart_name)) +
            xlim(-3.5, 3.5) +
            ylim(0, 0.7) +
            theme_bw() +
            theme(strip.placement = "outside", strip.background = element_blank()) + xlab(param) +
            theme(legend.position="top") +
            theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"))

        ggp
    }

    l<-lapply(par_vec,ff)
    l
}
