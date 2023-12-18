whittaker_biomes_plot <- function(data_MATp_MAPp) {


whittaker_base_plot()+
  # add the temperature - precipitation data points
  geom_point(data = data_MATp_MAPp, 
             aes(x = MATp, 
                 y = MAPp), 
             size   = 2,
             shape  = 21,
             colour = "gray95", 
             fill   = "black",
             stroke = 1,
             alpha  = 0.8) +
  theme_bw() +
    theme(
    #  legend.justification = c(0, 1), # pick the upper left corner of the legend box and
    #  legend.position = c(0, 1), # adjust the position of the corner as relative to axis
    #  legend.background = element_rect(fill = NA), # transparent legend background
    #  legend.box = "horizontal", # horizontal arrangement of multiple legends
    #  legend.spacing.x = unit(0.5, units = "cm"), # horizontal spacing between legends
      legend.position = "bottom",
      panel.grid = element_blank() # eliminate grids
    ) +
    facet_wrap(~ name)
}