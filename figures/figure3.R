library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)

example_data <- fread('output/spontaneous-movement/figures/event_data.csv')
sids <- example_data[ , unique(SID) ][1:2]
event_id <- 0

example_data <- example_data %>%
  .[ , Time := Time - unique(EventStart) ] %>%
  .[ , EventStart := 0 ] %>%
  .[ , Type := 'None' ] %>%
  .[ Time >= -0.4 & Time <= -0.2, Type := 'Baseline' ] %>%
  .[ Time >= -0.1 & Time <= 0, Type := 'Pre-movement' ] %>%
  .[ Time > 0 & Time <= 0.2, Type := 'Movement onset' ] %>%
  .[ Time > 0.2 & Time <= 0.4, Type := 'Late movement' ]

average_data <- example_data %>%
    .[ , list(TKEO = mean(TKEO), VM = mean(VM), EMG = mean(EMG)), by = list(Region, Time) ]


plot <- list()

for (region in c('S1_L23', 'S1_L5', 'M1_L23', 'M1_L5')) {
  sids <- example_data[ Region == region ][ , unique(SID) ][ c(4, 3) ]

  p1 <- example_data %>%
    .[ SID == sids[ 1 ] ] %>%
    .[ EventId == event_id ] %>%
    ggplot() +
    geom_line(aes(x = Time, y = scale(EMG)), alpha = 0.25) +
    geom_line(aes(x = Time, y = scale(TKEO, center = FALSE)), color = 'violet') +
    theme_void() +
    geom_vline(xintercept = 0, linewidth = 1, color = 'gray', linetype = 2)

  p2 <- example_data %>%
    .[ SID == sids[ 1 ] ] %>%
    .[ EventId == event_id ] %>%
    ggplot() +
    geom_line(aes(x = Time, y = VM), color = 'black') +
    theme_void() +
    geom_vline(xintercept = 0, linewidth = 1, color = 'gray', linetype = 2)

  p3 <- example_data %>%
    .[ SID == sids[ 2 ] ] %>%
    .[ EventId == event_id ] %>%
    ggplot() +
    geom_line(aes(x = Time, y = scale(EMG)), alpha = 0.25) +
    geom_line(aes(x = Time, y = scale(TKEO, center = FALSE)), color = 'violet') +
    theme_void() +
    geom_vline(xintercept = 0, linewidth = 1, color = 'gray', linetype = 2)

  p4 <- example_data %>%
    .[ SID == sids[ 2 ] ] %>%
    .[ EventId == event_id ] %>%
    ggplot() +
    geom_line(aes(x = Time, y = VM), color = 'black') +
    theme_void() +
    geom_vline(xintercept = 0, linewidth = 1, color = 'gray', linetype = 2)

  p5 <- average_data %>%
    .[ Region == region ] %>%
    ggplot() +
    geom_line(aes(x = Time, y = scale(EMG)), alpha = 0.25) +
    geom_line(aes(x = Time, y = scale(TKEO, center = FALSE)), color = 'violet') +
    theme_void() +
    geom_vline(xintercept = 0, linewidth = 1, color = 'gray', linetype = 2)

  p6 <- average_data %>%
    .[ Region == region ] %>%
    ggplot() +
    geom_line(aes(x = Time, y = VM), color = 'black') +
    theme_void() +
    geom_vline(xintercept = 0, linewidth = 1, color = 'gray', linetype = 2)

  region_average <- average_data[ Region == region ]

  p7 <- region_average %>%
    ggplot() +
    annotate('rect', xmin = -0.4, xmax = -0.2, ymin = region_average[ , min(VM)], ymax = region_average[ , max(VM)], fill='gray', size=0.5, alpha =0.2) +
    annotate('rect', xmin = -0.1, xmax = 0, ymin = region_average[ , min(VM)], ymax = region_average[ , max(VM)], fill='yellow', size=0.5, alpha =0.2) +
    annotate('rect', xmin = 0, xmax = 0.1, ymin = region_average[ , min(VM)], ymax = region_average[ , max(VM)], fill='red', size=0.5, alpha =0.2) +
    annotate('rect', xmin = 0.2, xmax = 0.4, ymin = region_average[ , min(VM)], ymax = region_average[ , max(VM)], fill='green', size=0.5, alpha =0.2) +
    geom_line(aes(x = Time, y = VM), color = 'black') +
    theme_void() +
    geom_vline(xintercept = 0, linewidth = 1, color = 'gray', linetype = 2)
  
  p8 <- example_data[ Region == region ] %>%
    .[ Type != 'None' ] %>%
    .[ Type == 'Baseline', Type := 'B' ] %>%
    .[ Type == 'Pre-movement', Type := 'PM' ] %>%
    .[ Type == 'Movement onset', Type := 'O' ] %>%
    .[ Type == 'Late movement', Type := 'LM' ] %>%
    .[ , Type := factor(Type, levels = c('B', 'PM', 'O', 'LM'))] %>%
    ggplot(aes(x = Type, y = VM, color = Type)) +
    geom_boxplot(outlier.alpha = 0) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    xlab('') +
    ylab('Vm') +
    stat_compare_means(aes(label = after_stat(p.signif)),
      comparisons = list(c('B', 'PM'), c('PM', 'O'), c('O', 'LM')),
      method = 't.test', position = -20) +
    theme(legend.position = 'none')

  p8

  plot[[region]] <- ggarrange(
    p1, p2, p3, p4, p5, p6, p7, p8,
    heights = c(0.5, 1, 0.5, 1, 0.5, 1, 2, 4),
    nrow = 8,
    ncol = 1
  )
}

final <- ggarrange(plotlist = plot, ncol = 4, labels = letters[1:4], font.label = list(size = 12))

ggsave(final, filename = 'plot.pdf', height = 6, width = 8, dpi = 1000)
