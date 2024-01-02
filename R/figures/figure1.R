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
library(scales)

colors <- wes_palette('Zissou1')

##################
# p1 (a)
##################

emg <- fread(snakemake@input$emg_data)
vm <- fread(snakemake@input$vm_data)

channel <- '2019_12_09t6I0_0'

channelData <- emg[ Channel == channel ] %>%
  merge(vm[ , list(VM, ChannelID, Time) ], by = c('ChannelID', 'Time'))

p1 <- channelData %>%
  ggplot() +
  geom_line(aes(x = Time, y = EMG), linewidth = 0.1, color = 'gray40') +
  geom_line(aes(x = Time, y = scale(VM, scale = FALSE) + 10), color = colors[ 5 ], linewidth = 0.1) +
  scale_y_continuous(
    name = 'EMG (mV)',
    sec.axis = sec_axis(~. + channelData[ , mean(VM) ] - 10, name = 'Membrane potential, V (mV)')
  ) +
  theme_light(base_size = 8) +
  xlab('Time, t (s)') +
  xlim(4.5, 7.5)

p1_legend <- get_legend(data.table(Type = c('EMG', 'V'), Point = c(1, 1, 2, 2)) %>%
                          ggplot() +
                          geom_line(aes(x = Point, y = Point, color = Type), linewidth = 1) +
                          theme_light(base_size = 8) +
                          scale_colour_manual(name = '', values = c('EMG' = 'gray40', 'V' = colors[ 5 ])) +
                          theme(legend.position = 'top'))

##################
# p2 (b)
##################

emgLimit <- channelData[ , TKEO ] %>% abs() %>% max()

p2 <- channelData %>%
  .[ , ScaledEMG := rescale(EMG, to = c(-emgLimit, emgLimit)) ] %>%
  ggplot() +
  geom_line(aes(x = Time, y = scale(ScaledEMG, scale = FALSE)), linewidth = 0.1, color = 'gray40') +
  geom_line(aes(x = Time, y = TKEO), linewidth = 0.5, color = colors[ 4 ]) +
  geom_hline(yintercept = 0.15, linewidth = 0.25, color = 'black') +
  theme_light(base_size = 8) +
  ylab('TKEO (a.u.)') +
  xlim(4.5, 7.5) +
  xlab('Time, t (s)') +
  annotate('text', x = 5, y = 0.3, label = 'Threshold 0.15 a.u.', size = 2.5)

p2_legend <- get_legend(data.table(Type = c('EMG', 'TKEO'), Point = c(1, 1, 2, 2)) %>%
                          ggplot() +
                          geom_line(aes(x = Point, y = Point, color = Type), linewidth = 1) +
                          theme_light(base_size = 8) +
                          scale_colour_manual(name = '', values = c('EMG' = 'gray40', 'TKEO' = colors[ 4 ])) +
                          theme(legend.position = 'top'))

##################
# p3 (c)
##################

movement_info <- fread(snakemake@input$movement_info)
data <- emg[ , Movement := 'None' ] %>%
  .[ , NextEMG := shift(EMG, type = 'lead') ] %>%
  .[ , NextTime := shift(Time, type = 'lead') ]

for (i in seq_len(nrow(movement_info))) {
  channel <- movement_info[ i, Channel ]
  start <- movement_info[ i, Start ]
  end <- movement_info[ i, End ]

  data[ Channel == channel & Time > start & Time < end, Movement := 'Move' ]
}

no_movement <- fread(snakemake@input$no_movement_information)

for (i in seq_len(nrow(no_movement))) {
  channel <- no_movement[ i, Channel ]
  start <- no_movement[ i, Start ]
  end <- no_movement[ i, End ]

  data[ Channel == channel & Time > start & Time < end, Movement := 'Rest' ]
}

p3 <- data %>%
  .[ Movement == 'None' ] %>%
  ggplot() +
  facet_wrap(~Channel, ncol = 1) +
  geom_line(aes(x = Time, y = EMG), color = 'gray80') +
  geom_segment(data = data[ Movement == 'Move' ], aes(x = Time, xend = NextTime, y = EMG, yend = NextEMG), color = colors[ 4 ]) +
  geom_segment(data = data[ Movement == 'Rest' ], aes(x = Time, xend = NextTime, y = EMG, yend = NextEMG), color = colors[ 1 ]) +
  theme_light(base_size = 8) +
  theme(strip.background = element_blank()) +
  theme(strip.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  theme(panel.grid = element_blank()) +
  theme(panel.spacing = unit(0, 'lines')) +
  xlab('Time, t (s)')

p3_legend <- get_legend(data.table(Type = c('Move', 'Rest'), Point = c(1, 1, 2, 2)) %>%
                          ggplot() +
                          geom_line(aes(x = Point, y = Point, color = Type), linewidth = 1) +
                          theme_light(base_size = 8) +
                          scale_colour_manual(name = '', values = c('Move' = colors[ 4 ], 'Rest' = colors[ 1 ])) +
                          theme(legend.position = 'top'))

final <- ggdraw() +
  draw_plot(p1, 0, 0.5, 0.4, 0.45) +
  draw_plot(p1_legend, 0, 0.96, 0.4, 0.04) +
  draw_plot(p2, 0, 0, 0.372, 0.45) +
  draw_plot(p2_legend, 0, 0.46, 0.4, 0.04) +
  draw_plot(p3, 0.4, 0, 0.6, 0.96) +
  draw_plot(p3_legend, 0.4, 0.96, 0.6, 0.04) +
  draw_plot_label(letters[ 1:3 ], x = c(0, 0, 0.4), y = c(0.99, 0.49, 0.99), size = 9)

ggsave(final, filename = snakemake@output$png, height = 6, width = 8, bg = 'white', dpi = 1000)
