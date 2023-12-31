library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

movement <- fread('figures/movement_info.csv') %>%
  .[ , .N, by = list(SID, Region) ] %>%
  .[ , Region := gsub(x = Region, pattern = '_', replacement = ' ') ] %>%
  .[ , Region := gsub(x = Region, pattern = '23', replacement = '2/3') ] %>%
  .[ , Type := 'Move' ]
no_movement <- fread('figures/no_movement_info.csv') %>%
  .[ , .N, by = list(SID, Region) ] %>%
  .[ , Region := gsub(x = Region, pattern = '_', replacement = ' ') ] %>%
  .[ , Region := gsub(x = Region, pattern = '23', replacement = '2/3') ] %>%
  .[ , Type := 'Rest' ]

dt <- rbind(movement, no_movement) %>%
  .[ , Region := factor(Region, levels = c('S1 L2/3', 'S1 L5', 'M1 L2/3', 'M1 L5')) ]

p1 <- dt %>%
  ggplot(aes(x = Region, y = N)) +
  facet_wrap(~Type, scale = 'free_y') +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(height = 0, width = 0.2, alpha = 0.5) +
  stat_compare_means(comparisons = list(
    c('S1 L2/3', 'S1 L5'),
    c('S1 L2/3', 'M1 L2/3'),
    c('S1 L2/3', 'M1 L5'),
    c('S1 L5', 'M1 L2/3'),
    c('S1 L5', 'M1 L5'),
    c('M1 L2/3', 'M1 L5')
  ), method = 't.test', size = 2.5) +
  theme_light(base_size = 8) +
  xlab('') +
  theme(legend.position = 'none') +
  ylab('Number of events') +
  theme(strip.background = element_rect(fill = 'white')) +
  theme(strip.text = element_text(colour = 'black'))

vm <- fread('figures/vm_data.csv')

averages <- vm %>%
 .[ , list(Average = mean(Vm), SD = sd(Vm)), by = list(SID, EventID, Type, Region) ] %>%
   .[ , Region := gsub(x = Region, pattern = '_', replacement = ' ') ] %>%
  .[ , Region := gsub(x = Region, pattern = '23', replacement = '2/3') ] %>%
  .[ , Region := factor(Region, levels = c('S1 L2/3', 'S1 L5', 'M1 L2/3', 'M1 L5')) ] %>%
  .[ Type == 'Movement', Type := 'Move' ] %>%
  .[ Type == 'No movement', Type := 'Rest'] %>%
  .[ , Type := factor(Type, c('Move', 'Rest')) ]

p2 <- averages %>%
  ggplot(aes(x = Type, y = Average)) +
  facet_wrap(~Region, ncol = 4) +
  geom_boxplot(outlier.alpha = 0.25) +
  theme_light(base_size = 8) +
  stat_compare_means(
    comparisons = list(c('Move', 'Rest')),
    method = 't.test',
    size = 2.5) +
  theme(strip.background = element_rect(fill = 'white')) +
  theme(strip.text = element_text(colour = 'black')) +
  xlab('') +
  ylab('Average membrane potential, V (mV)')

p3 <- averages %>%
  ggplot(aes(x = Type, y = SD)) +
  facet_wrap(~Region, ncol = 4) +
  geom_boxplot(outlier.alpha = 0.25) +
  theme_light(base_size = 8) +
  stat_compare_means(
    comparisons = list(c('Move', 'Rest')),
    method = 'wilcox.test',
    size = 2.5) +
  theme(strip.background = element_rect(fill = 'white')) +
  theme(strip.text = element_text(colour = 'black')) +
  xlab('') +
  ylab('Membrane potential standard deviation (mV)')


ap_data <- fread('figures/ap_data.csv') %>%
.[ , Region := gsub(x = Region, pattern = '_', replacement = ' ') ] %>%
  .[ , Region := gsub(x = Region, pattern = '23', replacement = '2/3') ] %>%
  .[ , Region := factor(Region, levels = c('S1 L2/3', 'S1 L5', 'M1 L2/3', 'M1 L5')) ] %>%
  .[ Type == 'Movement', Type := 'Move' ] %>%
  .[ Type == 'No movement', Type := 'Rest'] %>%
  .[ , Type := factor(Type, c('Move', 'Rest')) ]

p4 <- ap_data %>%
  ggplot(aes(x = Type, y = log(CountAP + 1))) +
  facet_wrap(~Region, ncol = 4, scale = 'free_y') +
  geom_boxplot(fill = 'white', alpha = 0.1, outlier.alpha = 0.2) +
  theme_light(base_size = 8) +
  stat_compare_means(
    comparisons = list(c('Move', 'Rest')),
    method = 'wilcox.test',
    size = 2.5) +
  theme(strip.background = element_rect(fill = 'white')) +
  theme(strip.text = element_text(colour = 'black')) +
  xlab('') +
  ylab('Number of AP, log')

final <- ggarrange(
  p1, p2, p3, p4,
  ncol = 2,
  nrow = 2,
  labels = letters[1:4],
  font.label = list(size = 9)
)

ggsave(final, filename = 'figure2.png', height = 6, width = 8, bg = 'white', dpi = 600)
