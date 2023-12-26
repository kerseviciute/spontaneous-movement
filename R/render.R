snakemake@scriptdir <- dirname(snakemake@scriptdir)
saveRDS(snakemake, paste0('.', basename(snakemake@params$script), '.RDS'))

script <- file.path(
  snakemake@scriptdir,
  snakemake@params$script
)

rmarkdown::render(
  input = script,
  knit_root_dir = getwd(),
  output_dir = dirname(snakemake@output$report),
  output_file = snakemake@output$report
)
