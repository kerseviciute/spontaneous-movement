saveRDS(snakemake, '.figure3.R.RDS')
# snakemake <- readRDS('.figure3.R.RDS')

library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(wesanderson)
library(cowplot)

colors <- wes_palette('Zissou1')

event_data <- fread(snakemake@input$event_data)

event_data <- event_data %>%
  .[ , Time := Time - unique(EventStart) ] %>%
  .[ , EventStart := 0 ] %>%
  .[ , Type := 'None' ] %>%
  .[ Time >= -0.4 & Time <= -0.2, Type := 'Baseline' ] %>%
  .[ Time >= -0.1 & Time <= 0, Type := 'Pre-movement' ] %>%
  .[ Time > 0 & Time <= 0.2, Type := 'Movement onset' ] %>%
  .[ Time > 0.2 & Time <= 0.4, Type := 'Late movement' ]

sample_data <- event_data[ , list(SID, Region) ] %>%
  unique()

average_data <- event_data %>%
  .[ , list(TKEO = mean(TKEO), VM = mean(VM), EMG = mean(EMG)), by = list(Region, Time) ]

plot <- list()
for (region in c('S1_L23', 'S1_L5', 'M1_L23', 'M1_L5')) {
  region_average <- average_data[ Region == region ]

  region_label <- gsub(x = region, pattern = '_', replacement = ' ')
  region_label <- gsub(x = region_label, pattern = '23', replacement = '2/3')

  max_vm <- region_average[ , max(VM) ] + 0.25

  p1 <- region_average %>%
    ggplot() +
    annotate('rect', xmin = -0.4, xmax = -0.2, ymin = -Inf, ymax = Inf, fill = colors[ 1 ], size = 0.5, alpha = 0.5) +
    annotate('rect', xmin = -0.1, xmax = 0, ymin = -Inf, ymax = Inf, fill = colors[ 2 ], size = 0.5, alpha = 0.5) +
    annotate('rect', xmin = 0, xmax = 0.1, ymin = -Inf, ymax = Inf, fill = colors[ 3 ], size = 0.5, alpha = 0.5) +
    annotate('rect', xmin = 0.2, xmax = 0.4, ymin = -Inf, ymax = Inf, fill = colors[ 4 ], size = 0.5, alpha = 0.5) +
    geom_line(aes(x = Time, y = VM), color = colors[ 5 ], linewidth = 0.25) +
    theme_light(base_size = 8) +
    ylab('Membrane\npotential, V (mV)') +
    ggtitle(region_label) +
    theme(plot.title = element_text(hjust = 0.5, size = 9, face = 'bold')) +
    scale_x_continuous(name = 'Time, t (s)', breaks = c(-0.4, -0.2, 0, 0.2, 0.4), limits = c(-0.4, 0.4)) +
    annotate('text', x = -0.3, y = max_vm, label = 'B', size = 2.5) +
    annotate('text', x = -0.05, y = max_vm, label = 'P', size = 2.5) +
    annotate('text', x = 0.05, y = max_vm, label = 'O', size = 2.5) +
    annotate('text', x = 0.3, y = max_vm, label = 'L', size = 2.5)

  type_average <- event_data[ Region == region ] %>%
    .[ , list(VM = mean(VM)), by = list(SID, Type) ] %>%
    .[ Type != 'None' ] %>%
    .[ Type == 'Baseline', Type := 'B' ] %>%
    .[ Type == 'Pre-movement', Type := 'P' ] %>%
    .[ Type == 'Movement onset', Type := 'O' ] %>%
    .[ Type == 'Late movement', Type := 'L' ] %>%
    .[ , Type := factor(Type, levels = c('B', 'P', 'O', 'L')) ] %>%
    .[ , Animal := gsub(x = SID, pattern = '(W[0-9]).*', replacement = '\\1') ] %>%
    .[ , Animal := factor(Animal, levels = c('W1', 'W2', 'W3', 'W4')) ]

  p2 <- type_average %>%
    ggpaired(x = 'Type', y = 'VM', id = 'SID',
             line.color = 'gray', line.size = 0.2, point.size = 0, color = 'white') +
    geom_boxplot(aes(color = Type), outlier.alpha = 0, fill = NA, linewidth = 0.25) +
    geom_point(aes(x = Type, y = VM, color = Type), alpha = 0.75, size = 0.75) +
    stat_compare_means(
      paired = TRUE,
      comparisons = list(c('B', 'P'), c('P', 'O'), c('B', 'L')),
      method = 'wilcox', size = 2.5) +
    stat_compare_means(
      paired = TRUE,
      comparisons = list(c('O', 'L')),
      method = 'wilcox', size = 2.5) +
    theme_light(base_size = 8) +
    theme(legend.position = 'none') +
    ylab('Membrane potential, V (mV)') +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    xlab('') +
    scale_colour_manual(values = c('B' = colors[ 1 ], 'P' = colors[ 2 ], 'O' = colors[ 3 ], 'L' = colors[ 4 ]))

  ap_time <- fread(snakemake@input$ap_time) %>%
    .[ , SID := gsub(x = Trial, pattern = '.*(W[0-9]+_C[0-9]+)', replacement = '\\1') ] %>%
    merge(sample_data) %>%
    .[ , Time := Time - 0.4 ] %>%
    .[ Region == region ] %>%
    .[ , Trial := as.factor(Trial) %>% as.numeric() ]

  p3 <- ap_time %>%
    ggplot() +
    annotate('rect', xmin = -0.4, xmax = -0.2, ymin = -Inf, ymax = Inf, fill = colors[ 1 ], size = 0.5, alpha = 0.5) +
    annotate('rect', xmin = -0.1, xmax = 0, ymin = -Inf, ymax = Inf, fill = colors[ 2 ], size = 0.5, alpha = 0.5) +
    annotate('rect', xmin = 0, xmax = 0.1, ymin = -Inf, ymax = Inf, fill = colors[ 3 ], size = 0.5, alpha = 0.5) +
    annotate('rect', xmin = 0.2, xmax = 0.4, ymin = -Inf, ymax = Inf, fill = colors[ 4 ], size = 0.5, alpha = 0.5) +
    geom_histogram(aes(x = Time), fill = 'white', color = 'black', binwidth = 0.025) +
    theme_light(base_size = 8) +
    ylab('AP count') +
    theme(strip.background = element_rect(fill = 'white')) +
    theme(strip.text = element_text(colour = 'black')) +
    scale_x_continuous(name = 'Time, t (s)', breaks = c(-0.4, -0.2, 0, 0.2, 0.4)) +
    expand_limits(y = c(0, 70)) +
    annotate('text', x = -0.3, y = 69, label = 'B', size = 2.5) +
    annotate('text', x = -0.05, y = 69, label = 'P', size = 2.5) +
    annotate('text', x = 0.05, y = 69, label = 'O', size = 2.5) +
    annotate('text', x = 0.3, y = 69, label = 'L', size = 2.5)

  p5 <- ggdraw() +
    draw_plot(p1, 0, 0.7, 1, 0.3) +
    draw_plot(p2, 0.068, 0.25, 1 - 0.068, 0.45) +
    draw_plot(p3, 0.08, 0, 1 - 0.08, 0.25)

  # Only for the first in list
  if (region == 'S1_L23') {
    p5 <- p5 +
      draw_plot_label(letters[ 1:3 ], x = c(0, 0, 0), y = c(1, 0.7, 0.25), size = 9)
  }

  plot[[region]] <- p5
}

final <- ggarrange(plotlist = plot, ncol = 4)

ggsave(final, filename = snakemake@output$png, height = 5, width = 8, bg = 'white', dpi = 1000)
