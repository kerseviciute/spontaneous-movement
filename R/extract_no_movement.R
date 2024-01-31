saveRDS(snakemake, '.extract_no_movement.R.RDS')
# snakemake <- readRDS('.extract_no_movement.R.RDS')

library(data.table, verbose = FALSE, quietly = TRUE, warn.conflicts = FALSE)
library(dplyr, verbose = FALSE, quietly = TRUE, warn.conflicts = FALSE)


data <- fread(snakemake@input$data)

tkeoThreshold <- snakemake@params$tkeoThreshold
maxTimeApart <- snakemake@params$maxTimeApart
minLength <- snakemake@params$minLength
expandBy <- snakemake@params$expandBy

channels <- data[ , unique(Channel) ]

events <- list()
for (channel in channels) {
  channelId <- data[ Channel == channel, unique(ChannelID) ]
  channelData <- data[ Channel == channel ] %>%
    .[ order(Time) ] %>%
    .[ , list(Time, TKEO, EMG) ] %>%
    .[ , Movement := TKEO < tkeoThreshold ] %>%
    .[ , MovementID := rleid(Movement) ] %>%
    .[ , MovementID := make.names(MovementID) ] %>%
    .[ , NextMovementID := shift(MovementID, type = 'lead') ]

  switchData <- channelData[ MovementID != NextMovementID ][ Movement == TRUE ][ order(Time) ]

  if (channelData[ , length(unique(MovementID)) ] > 2) {
    # Merge events that are very close in time
    for (i in 1:(nrow(switchData))) {
      end <- channelData[ MovementID == switchData[ i, MovementID ], max(Time) ]
      start <- channelData[ MovementID == switchData[ i, NextMovementID ], max(Time) ]

      if (end + maxTimeApart >= start) {
        mergeId <- switchData[ i, NextMovementID ]
        channelData[ MovementID == mergeId, Movement := TRUE ]
      }
    }
  }

  channelData <- channelData %>%
    .[ , MovementID := rleid(Movement) ] %>%
    .[ , MovementID := make.names(MovementID) ]

  channelEvents <- channelData %>%
    .[ Movement == TRUE, list(
      Start = min(Time),
      End = max(Time),
      Channel = channel,
      ChannelID = channelId
    ), by = list(MovementID) ] %>%
    # Expand the event start and end times
      .[ , Start := Start - expandBy ] %>%
      .[ , End := End + expandBy ] %>%
      .[ , Length := End - Start ] %>%
      .[ Length > 0 ] # Dropping events with negative time due to shortening

  events[[ channel ]] <- channelEvents
}

events <- rbindlist(events) %>%
  .[ , EventStart := round(Start * 20000) ] %>%
  .[ , EventEnd := round(End * 20000) ] %>%
  .[ , AnimalID := snakemake@wildcards$animal_id ] %>%
  .[ , CellName := snakemake@wildcards$cell_name ] %>%
  .[]

fwrite(events, snakemake@output$all_events)

finalEvents <- events %>%
  .[ Length >= minLength ] %>%
  .[ , MovementID := make.names(1:.N) ]

nrow(events)
nrow(finalEvents)

fwrite(finalEvents, snakemake@output$events)
