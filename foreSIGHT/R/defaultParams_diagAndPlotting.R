# This file contains the collated diagnostics and plotting parameters used by different parts of the code
# Include a function which can be used to read and modify the default settings in the package for reasonable parameters


tag_text <- c(expression(paste(italic("fore"), "SIGHT")))
tag_textCol <- "royalblue"

# Limits used to create traffic plots
#---------------------------------------------------

trafficLim <- list()
trafficLim[["pc.lim"]] <- c(5,10)             #SET FAIR AT 5-10 AND POOR 10+
trafficLim[["diff.lim"]] <- c(0.5,1)

# Colours used in traffic plots
traffic.col <- c("chartreuse3", "gold1", "red1")

foreSIGHT.colmap <- viridisLite::viridis

# Performance Space settings
# both
perfSpace_nlevel <- 256
perfSpace_contours <- TRUE
perfSpace_threshCol <- "white"
perfSpace_climDataBg <- "lightgray"
# Fill colour is changed if there are more than 300 points
perfSpace_climDataBg2 <- "black"
perfSpace_climDataCol <- "black"
perfSpace_alpha <- 0.9
# desired number of contours on the plot
perfSpace_nContour <- 3
# ggplot
threshLineSize <- 3
threshLabelSize <- 4.5
threshLabelText <- "Threshold"

# settings used by plotOptions - to compare two performance metrics & movement of contours
foreSIGHT.divColmap <- RColorBrewer::brewer.pal(n = 11, name = "BrBG")
option1_threshCol <- "white"
option1_lineSize <- 3
option1_lineAlpha <- 0.80
option2_threshCol <- "white"
option2_lineSize <- 3
option2_lineAlpha <- 1
# Way to create a diverging colour palette from the foreSIGHT.colMap
# endCol <- foreSIGHT.colmap(2)
# divColFn <- colorRampPalette(c(endCol[1], "white", endCol[2]))
# divCol <- divColFn(11)

# Base R - not used anymore
# perfSpace_nContour <- 5
threshLwd <- 2
threshLty <- 1
perfSpace_climDataPch <- 21


# Performance OAT settings

OATplot_lineSize <- 1
OATplot_fillAlpha <- 0.3
OATplot_textSize <- 11

heatPlot_textSize <- 12.5
heatPlot_margins <- c(0.5, 0.5, 0.5, 0.5)

theme_heatPlot <- function(textSize = heatPlot_textSize) {
  
  theme(axis.line = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  
  theme(panel.border = element_rect(colour = "black", size = 1, linetype = "solid", fill = NA),
        plot.margin = unit(heatPlot_margins, "cm"),
        plot.tag.position = c(0.95, 0.0)) +
    
  theme(plot.title = element_text(color = "black", size = textSize + 2, face = "plain", hjust = 0.5),
        axis.title.x = element_text(color = "black", size = textSize, angle = 0, hjust = .5, vjust = 0, face = "plain"),   # vjust = 2
        axis.title.y = element_text(color = "black", size = textSize, angle = 90, hjust = .5, vjust = 0, face = "plain", margin = margin(r = 10, unit = "pt")),  # vjust = -0.5, angle = 90
        axis.text.x = element_text(color = "black", size = textSize, face = "plain", vjust = 0),                                  
        axis.text.y = element_text(color = "black", size = textSize, face = "plain", hjust = 0, margin = margin(r = 2, unit = "pt")),
        plot.tag = element_text(color = tag_textCol, size = textSize)) +
  
  
  theme(legend.text = element_text(color = "black", size = textSize*0.9, face = "plain", margin = margin(r = 0.5, l = 0, unit = "cm")),
        legend.position = "bottom", legend.justification = "center",
        legend.title = element_text(color = "black", size = textSize*0.9, face = "plain", vjust = 1),
        
        # applicable only to the colorbar, set it there
        # legend.key.width = unit(1.5, "cm"), legend.key.height = unit(0.3, "cm"),
        
        # to remove the grey background colours from the legend
        legend.key = element_rect(fill = NA))

}



# contour colour?

# viridis - is it possible to set the full colours in viridis? - check
# default.colmap <- c("#440154FF", "#1F968BFF", "#FDE725FF")
# Chose colours using
# library(scales)
# show_col(viridis_pal()(20))

traffic_tileOutline <- "white"
traffic_textSize <- 12.5
# applies for no.of attributes <= 80
traffic_margins <- c(0.5,1,1,0)
traffic_upperMarg <- c(0.5, 1, 0, 0)
traffic_lowerMarg <- c(0, 1, 1, 0)

# for attributes > 80
traffic_tightMargins <- c(0.3,0.3,0.3,0)

# *** NOT USED - may be used for traffic light plots of individual targets
# The labels of attribute names on the heatmap (green, yellow = dark labels; red = light labels)
traffic_labelDarkCol <- "black"
traffic_labelLightCol <- "white"

# Anjana: For heatmaps we may use this hack to crates secondary dicrete x-axis in ggplot2: https://github.com/tidyverse/ggplot2/issues/3171
#         Or may add in the axis when this feature is implemented.


plotTitleTraffic <- c(Mean = "Mean of Absolute Biases", 
                      SD = "Standard Deviation of Absolute Biases")
legendTitleTraffic <- c(PType = "Bias relative to the target (in %)",
                        TempType = paste0("Bias relative to the target (in \u00B0C)"))

theme_traffic <- function(traffic_textSize) {

  theme(axis.line = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.ticks = element_blank()) +
  
  theme(legend.position = "bottom",legend.justification = "left",
        plot.tag.position = c(0.95, 0.0))    +  #%+replace%
  
  theme(plot.title = element_text(color = "black", size = traffic_textSize + 2, face = "plain", hjust = 0.5),
        axis.text.x = element_text(color = "black", size = traffic_textSize*0.8, angle = 270, hjust = 0, vjust = 0.2, face = "plain"),  # for full definition as label w/o geom_text
        axis.text.y = element_text(color = "black", size = traffic_textSize*0.8, angle = 0, hjust = 1, vjust = 0.2, face = "plain"),  
        axis.title.y = element_text(color = "black", size = traffic_textSize, angle = 90, hjust = .5, vjust = 1, face = "plain"),
        axis.title.x = element_text(color = "black", size = traffic_textSize, face = "italic", hjust = 1, vjust = 0),
        legend.text = element_text(color = "black", size = traffic_textSize*0.8, face = "plain"),
        legend.title = element_text(color = "black", size = traffic_textSize*0.8, face = "plain", vjust = 1),
        plot.tag = element_text(color = tag_textCol, size = traffic_textSize))

}

# theme to use if there is a lower plot
theme_traffic_upper <- function(traffic_textSize) {
  
  theme(axis.line = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.ticks = element_blank(), 
        axis.title.x = element_blank(), axis.title.y = element_blank()) +
    
    theme(legend.position = "bottom",legend.justification = "left") +
           #)plot.tag.position = c(0.95, 0.0))    +  #%+replace%
    
    theme(plot.title = element_text(color = "black", size = traffic_textSize + 2, face = "plain", hjust = 0.5, vjust = 1),
          axis.text.x = element_text(color = "black", size = traffic_textSize*0.8, angle = 90, hjust = 0, vjust = 0.2, face = "plain"),
          axis.text.y = element_text(color = "black", size = traffic_textSize*0.8, angle = 0, hjust = 1, vjust = 0.2, face = "plain"),  
          # axis.title.y = element_text(color = "black", size = traffic_textSize, angle = 90, hjust = .5, vjust = 1, face = "plain"),
          # axis.title.x = element_text(color = "black", size = traffic_textSize, face = "italic", hjust = 0),
          legend.text = element_text(color = "black", size = traffic_textSize*0.8, face = "plain"),
          legend.title = element_text(color = "black", size = traffic_textSize*0.8, face = "plain", vjust = 1)
          # plot.tag = element_text(color = tag_textCol, size = traffic_textSize
          )
}

theme_traffic_lower <- function(traffic_textSize) {
  
  theme(axis.line = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.ticks = element_blank(), plot.title = element_blank()) +
    
    theme(legend.position = "bottom",legend.justification = "left",
          plot.tag.position = c(0.95, 0.0))    +  #%+replace%
    
    theme(axis.text.x = element_text(color = "black", size = traffic_textSize*0.8, angle = 270, hjust = 0, vjust = 0.2, face = "plain"),  # for full definition as label w/o geom_text
          axis.text.y = element_text(color = "black", size = traffic_textSize*0.8, angle = 0, hjust = 1, vjust = 0.2, face = "plain"),  
          axis.title.y = element_text(color = "black", size = traffic_textSize, angle = 90, hjust = .5, vjust = 1, face = "plain"),
          axis.title.x = element_text(color = "black", size = traffic_textSize, face = "italic", hjust = 1, vjust = 0),
          legend.text = element_text(color = "black", size = traffic_textSize*0.8, face = "plain"),
          legend.title = element_text(color = "black", size = traffic_textSize*0.8, face = "plain", vjust = 1),
          plot.tag = element_text(color = tag_textCol, size = traffic_textSize))
  
}
