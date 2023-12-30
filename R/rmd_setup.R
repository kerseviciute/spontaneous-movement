library(reticulate)
library(knitr)
library(DT)
reticulate::use_condaenv('spontaneous-movement-mne')
knitr::knit_engines$set(python = reticulate::eng_python)
setOption <- options
setOption(knitr.graphics.rel_path = FALSE)

knitr::opts_chunk$set(
  context = "render",
  echo = FALSE,
  include = FALSE,
  out.width = "100%",
  fig.width = 7.2,
  fig.height = 4.3
)

showTable <- function(dt) {
  scrollY <- 50 + nrow(dt) * 21
  scrollY <- min(scrollY, 400)

  DT::datatable(
    dt,
    rownames = FALSE,
    width = '100%',
    style = 'bootstrap',
    class = 'table-condensed table-hover',
    extensions = c('Scroller'),
    escape = FALSE,
    options = list(
      scrollX = TRUE,
      fixedColumnts = TRUE,
      scrollY = scrollY,
      scroller = TRUE,
      dom = 'frtip',
      searching = FALSE
    )
  )
}
