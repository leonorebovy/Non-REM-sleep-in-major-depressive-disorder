## Script with common theme
## Project: Non-REM sleep in major depressive disorder
## Author: Leonore Bovy

options(warn = -1)
library(ggplot2)

raincloud_theme <- theme_classic() +  theme(
  text             = element_text(size = 10),
  axis.title.x     = element_text(size = 16),
  axis.title.y     = element_text(size = 16),
  axis.text.x      = element_text(size = 20, color="#000000"), 
  axis.text.y      = element_text(size = 20, color="#000000"),
  legend.title     = element_text(size = 16),
  legend.text      = element_text(size = 16),
  legend.position  = "left",
  plot.title       = element_text(lineheight = .8, face = "bold", size = 16),
  panel.border     = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x      = element_line(colour = "black", size = 1, linetype = "solid"),
  axis.line.y      = element_line(colour = "black", size = 1, linetype = "solid"))

my_theme <- theme_classic() + 
  theme(axis.title.x = element_text(size   = 16),
        axis.title.y = element_text(size   = 16),
        axis.text    = element_text(size   = 16, color="#000000"), 
        plot.title   = element_text(size   = 18, hjust = 0.5),
        axis.line.x  = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.line.y  = element_line(colour = "black", size = 1, linetype = "solid"))

my_theme_freq <- theme_classic() + 
  theme(axis.title.x = element_text(size   = 16),
        axis.title.y = element_text(size   = 16),
        axis.text    = element_text(size   = 16, color="#000000"), 
        plot.title   = element_text(size   = 18, hjust = 0.5),
        axis.line.x  = element_line(colour = "black", size = 1.2, linetype = "solid"),
        axis.line.y  = element_line(colour = "black", size = 1.2, linetype = "solid"))

pretty_regressionplot <- function(dataframe, x, y, grouping) {
  ggplot(dataframe, aes(x = x, y = y, colour = grouping)) + 
    geom_point(data = dataframe, aes(y = y, group = grouping), size = 4) +
    geom_smooth(method = "lm", se = FALSE, fullrange=TRUE, size = 1.2) +   
    theme_classic() +                   
    theme(plot.title = element_text(hjust = 0.5)) +    
    theme(axis.text.x = element_text(size=20, color="#000000", family="serif" ), 
          axis.text.y = element_text(size=20, color="#000000", family="serif"), 
          plot.title = element_text(size = 15,family="serif"),
          axis.title.x = element_text(size=17,family="serif", margin = margin(b=10)),
          axis.title.y = element_text(size=17,family="serif"),
          axis.line.x  = element_line(colour = "black", size = 1),
          axis.line.y  = element_line(colour = "black", size = 1))}

color_dataset_a     = (c("#646464", "#b20000"))
color_dataset_b     = (c("#646464", "#66c4ff"))
color_dataset_b_pat = (c("#66c4ff", "#009dff"))
color_dataset_b_all = (c("#646464", "#66c4ff", "#009dff"))
color_dataset_c     = (c("#646464", "#7fcf7f"))
color_dataset_c_pat = (c("#7fcf7f", "#00a000"))
color_dataset_c_all = (c("#646464", "#7fcf7f", "#00a000"))
color_dataset_all   = (c("#646464", "#b20000", "#646464", "#66c4ff", "#009dff", "#646464", "#7fcf7f", "#00a000"))
color_dataset_three = (c("#646464", "#b20000", "#009dff", "#646464", "#7fcf7f", "#00a000"))
color_dataset_hamd  = (c("#b20000", "#66c4ff", "#7fcf7f"))
color_dataset_beh   = (c("#646464", "#b20000","#646464", "#66c4ff", "#009dff")) # A and B
color_dataset_comb  = (c("#646464", "#871F78")) #purple
