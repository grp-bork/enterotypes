addSmallLegend <- function(myPlot, pointSize = 0.5, textSize = 3, spaceLegend = 0.1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

enterotype_colors  = c("Prevotella"= "#DF6929",
                       "Firmicutes" = "#197EA5",
                       "Bact_Phoc" = "#2B9348", 
                       "Firmicutes_1"="#197EA5", 
                       "Firmicutes_2" = "#46A0C3", 
                       "Firmicutes_3" = "#7BC4E1",
                       "Bact_Phoc_Firm" =  "#2BA187",
                       "dysbiosis" = "#aa3e98"
)

enterotype_shapes  = c("Prevotella"= 21,
                       "Firmicutes" = 22,
                       "Bact_Phoc" = 23, 
                       "Firmicutes_1"=22, 
                       "Firmicutes_2" = 22, 
                       "Firmicutes_3" = 22,
                       "Bact_Phoc_Firm" =  25
)

phyla_colors = c("Firmicutes" = "#197ea5",
                      "Bacteroidota" = "#588157",
                      "Actinobacteriota" = "#b5838d",
                      "Proteobacteria" = "#DF692A",
                      "Verrucomicrobiota" = "#ffd166",
                      "Firmicutes_C" = "#46a0c3",
                      "Firmicutes_A" = "#7bc4e1")

disease_colors = c("patients" = "#aa3e98", "healthy subjects" = "#ffce1f" , "unknown disease status" = "#9DA39A")#
disease_colors_poster = c("diseased" = "#193F90", "healthy" = "#B65417")
dysbiosis_colors = c("no dysbiosis" = "#193F90", "dysbiosis" = "#B65417")
#disease_colors = c("diseased" = "#C48FCC", "healthy" = "#95F5D3" )#9DA39A

#library(microshades)
#ms_colors_1 <- c(microshades_palette("micro_green", 1, lightest = FALSE), 
#                 microshades_palette("micro_blue", 1, lightest = FALSE), 
#                 microshades_palette("micro_purple", 1, lightest = FALSE), 
#                 microshades_palette("micro_orange", 1, lightest = FALSE), 
#                 microshades_palette("micro_brown", 1, lightest = FALSE), 
#                 microshades_palette("micro_gray", 1, lightest = FALSE))
#
#ms_colors_2 <- c(microshades_palette("micro_green", 2, lightest = FALSE), 
#                 microshades_palette("micro_blue", 2, lightest = FALSE), 
#                 microshades_palette("micro_purple", 2, lightest = FALSE), 
#                 microshades_palette("micro_orange", 2, lightest = FALSE), 
#                 microshades_palette("micro_brown", 2, lightest = FALSE), 
#                 microshades_palette("micro_gray", 2, lightest = FALSE))
#

#ms_colors_5 <- c(microshades_palette("micro_green", 5, lightest = FALSE), 
#              microshades_palette("micro_blue", 5, lightest = FALSE), 
#              microshades_palette("micro_purple", 5, lightest = FALSE), 
#              microshades_palette("micro_orange", 5, lightest = FALSE), 
#              microshades_palette("micro_brown", 5, lightest = FALSE), 
#              microshades_palette("micro_gray", 5, lightest = FALSE))
#

colors_blue_red <- c("#264653", "#287271", "#2A9D8F", "#8AB17D", "#BABB74", "#E9C46A", "#EFB366", "#EE8959", "#E76F51", "#345D58", "#677B61", "#7F805B",
	       "#A68D50", "#AD834C", "#B17748", "#A66649", "#9C5848")


color_values <-c("#E76F51","#E9C46A","#F4A261","#F3DFC1","#2A9D8F","#457B9D","#A8DADC","#CAD2C5","#84A98C","#52796F","#354F52","#99b31b", "#6D6875", "#B5838D","#E5989B","#FFB4A2","#FFCDB2","#9d0208","#264653","#e9c46a", "#e76f51", "#90be6d", "#723c70","#2a9d8f")
