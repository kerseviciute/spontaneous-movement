saveRDS(snakemake, '.figure1.R.RDS')
# snakemake <- readRDS('.figure1.R.RDS')

library(data.table)
library(dplyr)
library(ggplot2)
library(wesanderson)
library(foreach)
library(ggpubr)
library(cowplot)
library(gridGraphics)

colors <- wes_palette('Zissou1')

data <- fread(snakemake@input$data) %>%
  .[ , Color := 'Reg' ] %>%
  .[ Movement == TRUE, Color := 'Movement' ] %>%
  .[ NoMovement == TRUE, Color := 'No Movement' ] %>%
  .[ , NextTime := shift(Time) ] %>%
  .[ , NextEMG := shift(EMG) ]

limit <- data[ , max(abs(EMG)) ]

plots <- foreach(channel = data[ , unique(Channel) ]) %do% {
  p <- data[ Channel == channel ] %>%
    ggplot() +
    geom_segment(aes(x = Time, xend = NextTime, y = EMG, yend = NextEMG, color = Color), linewidth = 0.1) +
    theme_void(base_size = 8) +
    scale_colour_manual(name = 'Event', values = c('Movement' = colors[ 4 ], 'No Movement' = colors[ 1 ])) +
    theme(legend.position = 'None') +
    scale_y_continuous(breaks = 0, labels = channel, limits = c(-12, 12)) +
    xlim(0, 10) +
    theme(axis.text.y = element_text())
}

p1 <- ggarrange(plotlist = plots, ncol = 1)

p1_legend <- get_legend(data %>%
                          ggplot() +
                          geom_segment(aes(x = Time, xend = NextTime, y = EMG, yend = NextEMG, color = Color), linewidth = 1) +
                          theme_void(base_size = 8) +
                          scale_colour_manual(name = '', values = c('Movement' = colors[ 4 ], 'No Movement' = colors[ 1 ])) +
                          theme(legend.position = 'None') +
                          scale_y_continuous(breaks = 0, labels = channel, limits = c(-12, 12)) +
                          xlim(0, 10) +
                          theme(axis.text.y = element_text()) +
                          theme(legend.position = 'top', legend.text = element_text(size = 8)))

event_data <- fread(snakemake@input$event_data)

average_data <- event_data %>%
  .[ , list(TKEO = mean(TKEO), VM = mean(VM), EMG = mean(EMG)), by = Time ] %>%
  .[ , Time := Time - 0.4 ]

p2 <- average_data %>%
  ggplot() +
  geom_line(aes(x = Time, y = scale(EMG)), linewidth = 0.1, color = 'gray40') +
  geom_line(aes(x = Time, y = scale(TKEO, center = FALSE)), color = colors[ 5 ], linewidth = 0.35) +
  theme_light(base_size = 8) +
  ylab('EMG, scaled') +
  xlab('Time, s')

p2_legend <- get_legend(data.table(Type = c('EMG', 'TKEO'), Point = c(1, 1, 2, 2)) %>%
                          ggplot() +
                          geom_line(aes(x = Point, y = Point, color =  Type), linewidth = 1) +
                          theme_light(base_size = 8) +
                          scale_colour_manual(name = '', values = c('EMG' = 'gray40', 'TKEO' = colors[ 5 ])) +
                          theme(legend.position = 'top'))

p3 <- data[ Channel == '2019_12_18t16I0_18' ] %>%
  ggplot() +
  geom_line(aes(x = Time, y = scale(EMG)), linewidth = 0.1, color = 'gray40') +
  geom_line(aes(x = Time, y = scale(VM, scale = FALSE)), color = colors[ 5 ], linewidth = 0.1) +
  scale_y_continuous(
    name = 'EMG, scaled',
    sec.axis = sec_axis(~. + data[ Channel == channel ][ , mean(VM) ], name = 'Voltage, V (mV)')
  ) +
  theme_light(base_size = 8) +
  xlab('Time, s')

p3_legend <- get_legend(data.table(Type = c('EMG', 'V'), Point = c(1, 1, 2, 2)) %>%
                          ggplot() +
                          geom_line(aes(x = Point, y = Point, color = Type), linewidth = 1) +
                          theme_light(base_size = 8) +
                          scale_colour_manual(name = '', values = c('EMG' = 'gray40', 'V' = colors[ 5 ])) +
                          theme(legend.position = 'top'))

final <- ggdraw() +
  draw_plot(p1, 0.0, 0, 0.6, 0.95) +
  draw_plot(p1_legend, 0.0, 0.95, 0.6, 0.05) +
  draw_plot(p2, 0.6, 0.5, 0.4, 0.45) +
  draw_plot(p2_legend, 0.6, 0.96, 0.4, 0.04) +
  draw_plot(p3, 0.6, -0.01, 0.4, 0.45) +
  draw_plot(p3_legend, 0.6, 0.45, 0.4, 0.04) +
  draw_plot_label(letters[ 1:3 ], x = c(0, 0.6, 0.6), y = c(1, 1, 0.5), size = 9)

ggsave(final, filename = snakemake@output$png, height = 4, width = 8, bg = 'white')

# TODO: add time scale in A
