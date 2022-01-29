#' theme_pub
#'
#' @param base_size font size
#' @param base_family font name
#'
#' @import grid
#' @import ggthemes


theme_pub <- function(base_size=14, base_family="helvetica") {

  (theme_foundation(base_size=base_size,
                    base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2),
                                      hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    )
   )

}

#' scale_fill_pub
#'
#' @param n number pf colours
#' @param color.theme color themes of RColorBrewer
#' @param ... ggplot2 extension
#' @import scales
scale_fill_pub <- function(n, color.theme, ...){
  discrete_scale("fill",
                 "Publication",
                 manual_pal(values = brewer.pal(n, color.theme)), ...)

}

#' scale_colour_pub
#'
#' @param n number pf colours
#' @param color.theme color themes of RColorBrewer
#' @param ... ggplot2 extension
#' @import scales
scale_colour_pub <- function(n, color.theme, ...){
  discrete_scale("colour",
                 "Publication",
                 manual_pal(values = brewer.pal(n, color.theme)), ...)

}
