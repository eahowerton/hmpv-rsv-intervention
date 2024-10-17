## plot conceptual diagram to illustrate possible mechanisms of interaciton
library(dplyr)
library(ggplot2)
library(RColorBrewer)

rsv_vl = data.frame(x = c(2, 5, 16), 
                    y = c(0, 50, 0), 
                    id = "rsv")
mpv_vl = data.frame(x = c(11, 13.5, 23), 
                     y = c(0, 40, 0), 
                    id = "mpv")
innate = data.frame(x = c(3, 6, 20), 
                    y = c(0, 58, 0), 
                    id = "innate")
mpv_winteraction = data.frame(x = c(12, 15, 26), 
                                 y = c(0, 32, 0), 
                                 id = "mpv_winteraction")
all_curves = bind_rows(rsv_vl, mpv_vl, innate, mpv_winteraction)

ggplot(data = all_curves, aes(x = x, y = y, color = id, linetype = id)) + 
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = c(2, 11), 
                     expand = c(0,0), 
                     labels = c("RSV\ninfection", "HMPV\ninfection"),
                     limits = c(1, max(all_curves$x)), name = "time") + 
  scale_y_continuous(expand = c(0,0), name = "log(viral load)") + 
  scale_color_manual(values = c("darkgray", brewer.pal(3, "Dark2")[1], brewer.pal(3, "Dark2")[1:2])) + 
  scale_linetype_manual(values = c("solid", "solid", "dashed", "solid")) + 
  theme_classic() +
  theme(axis.text.y = element_blank(), 
        axis.ticks = element_blank(), 
        legend.key.width = unit(1, "cm"),
        legend.position = c(0.9,0.8),
        legend.title = element_blank(),
        panel.grid = element_blank())
ggsave("figures/conceptual_figure.pdf", width = 6, height = 3)

