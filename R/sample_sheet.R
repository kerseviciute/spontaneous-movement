saveRDS(snakemake, '.sample_sheet.R.RDS')
# snakemake <- readRDS('.sample_sheet.R.RDS')

library(data.table)
library(dplyr)
library(openxlsx)

data <- openxlsx::read.xlsx(snakemake@input$data, sheet = 3, colNames = FALSE) %>% setDT()

data <- rbind(
  data[ 3:12, list(SID = X1, Region = data[ 1, X1 ], Date = X2, Depth = X3, Count = X4) ],
  data[ 3:12, list(SID = X5, Region = data[ 1, X5 ], Date = X6, Depth = X7, Count = X8) ],
  data[ 3:14, list(SID = X9, Region = data[ 1, X9 ], Date = X10, Depth = X11, Count = X12) ],
  data[ 3:12, list(SID = X13, Region = data[ 1, X13 ], Date = X14, Depth = X15, Count = X16) ]
) %>%
  .[ , AnimalID := gsub(x = SID, pattern = '(W[0-9]+)_.+', replacement = '\\1') ] %>%
  .[ , CellName := gsub(x = SID, pattern = 'W[0-9]+_(C[0-9]+)', replacement = '\\1') ] %>%
  data.table::setcolorder(c('SID', 'AnimalID', 'CellName')) %>%
  .[ , Region := gsub(x = Region, pattern = '/', replacement = '') ]

head(data)

fwrite(data, snakemake@output$sample_sheet)
saveRDS(data, snakemake@output$samples)
