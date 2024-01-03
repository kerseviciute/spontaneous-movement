saveRDS(snakemake, '.figure2.R.RDS')
# snakemake <- readRDS('.figure2.R.RDS')

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(foreach)
library(RColorBrewer)

samples <- fread(snakemake@input$sample_sheet)

movement <- foreach(file = snakemake@input$movement_information, .combine = rbind) %do% {
  fread(file)
} %>%
  merge(samples) %>%
  .[ , .N, by = list(SID, Region) ] %>%
  .[ , Region := gsub(x = Region, pattern = '_', replacement = ' ') ] %>%
  .[ , Region := gsub(x = Region, pattern = '23', replacement = '2/3') ] %>%
  .[ , Type := 'Move' ]

no_movement <- fread(snakemake@input$rest_information) %>%
  .[ , .N, by = list(SID, Region) ] %>%
  .[ , Region := gsub(x = Region, pattern = '_', replacement = ' ') ] %>%
  .[ , Region := gsub(x = Region, pattern = '23', replacement = '2/3') ] %>%
  .[ , Type := 'Rest' ]

###############
# P1 (a)
###############

dt <- rbind(movement, no_movement) %>%
  .[ , Region := factor(Region, levels = c('S1 L2/3', 'S1 L5', 'M1 L2/3', 'M1 L5')) ] %>%
  .[ , Animal := gsub(x = SID, pattern = '(W[0-9]).*', replacement = '\\1') ] %>%
  .[ , Animal := factor(Animal, levels = c('W1', 'W2', 'W3', 'W4')) ]

p1 <- dt %>%
  .[ , Type := factor(Type, levels = c('Rest', 'Move')) ] %>%
  ggplot(aes(x = Region, y = N)) +
  facet_wrap(~Type, scale = 'free_y') +
  geom_boxplot(outlier.alpha = 0, linewidth = 0.25) +
  geom_jitter(aes(color = Animal), height = 0, width = 0.2, alpha = 0.5) +
  stat_compare_means(comparisons = list(
    c('S1 L2/3', 'S1 L5'),
    c('S1 L5', 'M1 L2/3'),
    c('S1 L2/3', 'M1 L2/3'),
    c('S1 L5', 'M1 L5'),
    c('S1 L2/3', 'M1 L5')
  ), method = 'wilcox', size = 2.5) +
  stat_compare_means(comparisons = list(
    c('M1 L2/3', 'M1 L5')
  ), method = 'wilcox', size = 2.5) +
  theme_light(base_size = 8) +
  xlab('') +
  theme(legend.position = 'top') +
  ylab('Number of events') +
  theme(strip.background = element_rect(fill = 'white')) +
  theme(strip.text = element_text(colour = 'black')) +
  scale_color_manual(name = 'Animal', values = brewer.pal(11, 'RdBu')[ c(1, 2, 10, 11) ]) +
  scale_fill_manual(name = 'Animal', values = brewer.pal(11, 'RdBu')[ c(1, 2, 10, 11) ]) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

###############
# P2 (b)
###############

vm <- fread(snakemake@input$vm)

averages <- vm %>%
  .[ , list(Average = mean(Vm), SD = sd(Vm)), by = list(SID, Type, Region) ] %>%
  .[ , Region := gsub(x = Region, pattern = '_', replacement = ' ') ] %>%
  .[ , Region := gsub(x = Region, pattern = '23', replacement = '2/3') ] %>%
  .[ , Region := factor(Region, levels = c('S1 L2/3', 'S1 L5', 'M1 L2/3', 'M1 L5')) ] %>%
  .[ Type == 'Movement', Type := 'Move' ] %>%
  .[ Type == 'No movement', Type := 'Rest' ] %>%
  .[ , Type := factor(Type, c('Move', 'Rest')) ] %>%
  .[ , Animal := gsub(x = SID, pattern = '(W[0-9]).*', replacement = '\\1') ] %>%
  .[ , Animal := factor(Animal, levels = c('W1', 'W2', 'W3', 'W4')) ] %>%
  .[ , Type := factor(Type, levels = c('Rest', 'Move')) ]

p2 <- averages %>%
  ggpaired(x = 'Type', y = 'Average', id = 'SID',
           line.color = 'gray', line.size = 0.2, point.size = 0, color = 'white') +
  facet_wrap(~Region, ncol = 4) +
  geom_boxplot(outlier.alpha = 0, fill = NA, linewidth = 0.25) +
  geom_point(aes(x = Type, y = Average, color = Animal), alpha = 0.5) +
  theme_light(base_size = 8) +
  stat_compare_means(
    paired = TRUE,
    comparisons = list(c('Move', 'Rest')),
    method = 'wilcox',
    size = 2.5) +
  theme(strip.background = element_rect(fill = 'white')) +
  theme(strip.text = element_text(colour = 'black')) +
  theme(legend.position = 'top') +
  xlab('') +
  ylab('Membrane potential, V (mV)') +
  scale_color_manual(name = 'Animal', values = brewer.pal(11, 'RdBu')[ c(1, 2, 10, 11) ]) +
  scale_fill_manual(name = 'Animal', values = brewer.pal(11, 'RdBu')[ c(1, 2, 10, 11) ]) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

###############
# P3 (c)
###############

p3 <- averages %>%
  ggpaired(x = 'Type', y = 'SD', id = 'SID',
           line.color = 'gray', line.size = 0.2, point.size = 0, color = 'white') +
  facet_wrap(~Region, ncol = 4) +
  geom_boxplot(outlier.alpha = 0, fill = NA, linewidth = 0.25) +
  geom_point(aes(x = Type, y = SD, color = Animal), alpha = 0.5) +
  theme_light(base_size = 8) +
  stat_compare_means(
    paired = TRUE,
    comparisons = list(c('Move', 'Rest')),
    method = 'wilcox',
    size = 2.5) +
  theme(strip.background = element_rect(fill = 'white')) +
  theme(strip.text = element_text(colour = 'black')) +
  theme(legend.position = 'top') +
  xlab('') +
  ylab('Membrane potential SD, V (mV)') +
  scale_color_manual(name = 'Animal', values = brewer.pal(11, 'RdBu')[ c(1, 2, 10, 11) ]) +
  scale_fill_manual(name = 'Animal', values = brewer.pal(11, 'RdBu')[ c(1, 2, 10, 11) ]) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

###############
# P4 (d)
###############

p4 <- averages %>%
  ggplot(aes(x = Region, y = Average)) +
  facet_wrap(~Type) +
  geom_boxplot(outlier.alpha = 0, linewidth = 0.25) +
  geom_jitter(aes(color = Animal), height = 0, width = 0.2, alpha = 0.5) +
  theme(legend.position = 'none') +
  stat_compare_means(comparisons = list(
    c('S1 L2/3', 'S1 L5'),
    c('S1 L5', 'M1 L2/3'),
    c('S1 L2/3', 'M1 L2/3'),
    c('S1 L5', 'M1 L5'),
    c('S1 L2/3', 'M1 L5')
  ), method = 'wilcox', size = 2.5) +
  stat_compare_means(comparisons = list(
    c('M1 L2/3', 'M1 L5')
  ), method = 'wilcox', size = 2.5) +
  theme_light(base_size = 8) +
  xlab('') +
  theme(legend.position = 'none') +
  ylab('Membrane potential, V (mV)') +
  theme(strip.background = element_rect(fill = 'white')) +
  theme(strip.text = element_text(colour = 'black')) +
  scale_color_manual(name = '', values = brewer.pal(11, 'RdBu')[ c(1, 2, 10, 11) ]) +
  scale_fill_manual(name = '', values = brewer.pal(11, 'RdBu')[ c(1, 2, 10, 11) ]) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

###############
# P5 (e)
###############

ap_data <- fread(snakemake@input$ap_count) %>%
  .[ , Region := gsub(x = Region, pattern = '_', replacement = ' ') ] %>%
  .[ , Region := gsub(x = Region, pattern = '23', replacement = '2/3') ] %>%
  .[ , Region := factor(Region, levels = c('S1 L2/3', 'S1 L5', 'M1 L2/3', 'M1 L5')) ] %>%
  .[ Type == 'Movement', Type := 'Move' ] %>%
  .[ Type == 'No movement', Type := 'Rest' ] %>%
  .[ , Type := factor(Type, c('Move', 'Rest')) ] %>%
  .[ , list(CountAP = mean(CountAP)), by = list(SID, Type, Region) ] %>%
  .[ , Animal := gsub(x = SID, pattern = '(W[0-9]).*', replacement = '\\1') ] %>%
  .[ , Animal := factor(Animal, levels = c('W1', 'W2', 'W3', 'W4')) ] %>%
  .[ , Type := factor(Type, levels = c('Rest', 'Move')) ]

p5 <- ap_data %>%
  .[ , y := log(CountAP + 1) ] %>%
  ggpaired(x = 'Type', y = 'y', id = 'SID',
           line.color = 'gray', line.size = 0.2, point.size = 0, color = 'white') +
  facet_wrap(~Region, ncol = 4, scale = 'free_y') +
  geom_boxplot(outlier.alpha = 0, fill = NA, linewidth = 0.25) +
  geom_point(aes(x = Type, y = y, color = Animal), alpha = 0.5) +
  theme_light(base_size = 8) +
  stat_compare_means(
    paired = TRUE,
    comparisons = list(c('Move', 'Rest')),
    method = 'wilcox',
    size = 2.5) +
  theme(strip.background = element_rect(fill = 'white')) +
  theme(strip.text = element_text(colour = 'black')) +
  theme(legend.position = 'none') +
  xlab('') +
  ylab('Number of AP, log') +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_color_manual(name = 'Animal', values = brewer.pal(11, 'RdBu')[ c(1, 2, 10, 11) ]) +
  scale_fill_manual(name = 'Animal', values = brewer.pal(11, 'RdBu')[ c(1, 2, 10, 11) ])

p123 <- ggarrange(
  p1, p2, p3,
  widths = c(0.4, 0.3, 0.3),
  ncol = 3,
  nrow = 1,
  labels = letters[ 1:3 ],
  font.label = list(size = 9),
  common.legend = TRUE
)

p45 <- ggarrange(
  p4, p5,
  ncol = 2,
  nrow = 1,
  labels = letters[ 4:5 ],
  font.label = list(size = 9)
)

final <- ggarrange(
  p123, p45,
  heights = c(0.55, 0.45),
  ncol = 1,
  nrow = 2
)

ggsave(final, filename = snakemake@output$png, height = 6, width = 8, bg = 'white', dpi = 1000)
