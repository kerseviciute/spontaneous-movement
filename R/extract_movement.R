saveRDS(snakemake, '.extract_movement.R.RDS')
# snakemake <- readRDS('.extract_movement.R.RDS')

library(data.table, verbose = FALSE, quietly = TRUE, warn.conflicts = FALSE)
library(dplyr, verbose = FALSE, quietly = TRUE, warn.conflicts = FALSE)

rms_amplitude <- function(signal) {
  signal^2 %>% mean() %>% sqrt()
}


data <- fread(snakemake@input$data)

tkeoThreshold <- snakemake@params$tkeoThreshold
maxTimeApart <- snakemake@params$maxTimeApart
maxLength <- snakemake@params$maxLength
minLength <- snakemake@params$minLength
calmBeforeEvent <- snakemake@params$calmBeforeEvent
minAmplitude <- snakemake@params$minAmplitude

channels <- data[ , unique(Channel) ]

events <- list()
for (channel in channels) {
  channelId <- data[ Channel == channel, unique(ChannelID) ]
  channelData <- data[ Channel == channel ] %>%
    .[ order(Time) ] %>%
    .[ , list(Time, TKEO, EMG) ] %>%
    .[ , Movement := TKEO > tkeoThreshold ] %>%
    .[ , MovementID := rleid(Movement) ] %>%
    .[ , MovementID := make.names(MovementID) ] %>%
    .[ , NextMovementID := shift(MovementID, type = 'lead') ]

  switchData <- channelData[ MovementID != NextMovementID ][ Movement == TRUE ][ order(Time) ]

  if (nrow(switchData) == 0) { next }
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

  # Calculate amplitudes
  for (event in channelData[ , unique(MovementID) ]) {
    signal <- channelData[ MovementID == event ][ order(Time), EMG ]
    channelData[ MovementID == event, Amplitude := rms_amplitude(signal) ]

    startTime <- channelData[ MovementID == event, min(Time) ]
    signal <- channelData[ Time >= startTime & Time < startTime + 0.05 ][ order(Time), EMG ]
    channelData[ MovementID == event, StartAmplitude := rms_amplitude(signal) ]
  }

  channelEvents <- channelData %>%
    .[ Movement == TRUE, list(
      Start = min(Time),
      End = max(Time),
      Amplitude = unique(Amplitude),
      StartAmplitude = unique(StartAmplitude),
      Channel = channel,
      ChannelID = channelId
    ), by = list(MovementID) ] %>%
    # Add information about the start of the next event
    .[ order(Start), PreviousEndedBefore := abs(shift(End, type = 'lag') - Start) ] %>%
    .[ is.na(PreviousEndedBefore), PreviousEndedBefore := Start ] %>%
    # Add information about the end of the previous event
    .[ order(Start), NextStartedAfter := abs(shift(Start, type = 'lead')) - End ] %>%
    .[ is.na(NextStartedAfter), NextStartedAfter := 10.0 - End ] %>%
    .[ , Length := End - Start ]

  events[[ channel ]] <- channelEvents
}

events <- rbindlist(events) %>%
  .[ , EventStart := round(Start * 20000) ] %>%
  .[ , EventEnd := round(End * 20000) ] %>%
  .[ , AnimalID := snakemake@wildcards$animal_id ] %>%
  .[ , CellName := snakemake@wildcards$cell_name ]

fwrite(events, snakemake@output$all_events)

finalEventsStart <- events %>%
  .[ Length >= minLength ] %>%
  .[ Length <= maxLength ] %>%
  .[ PreviousEndedBefore > calmBeforeEvent ] %>%
  .[ , MovementID := make.names(1:.N) ] %>%
  .[ Amplitude > minAmplitude ]

finalEventsEnd <- events %>%
  .[ Length >= minLength ] %>%
  .[ Length <= maxLength ] %>%
  .[ NextStartedAfter > calmBeforeEvent ] %>%
  .[ , MovementID := make.names(1:.N) ] %>%
  .[ Amplitude > minAmplitude ]

nrow(events)
nrow(finalEventsStart)
nrow(finalEventsEnd)

fwrite(finalEventsStart, snakemake@output$events_start)
fwrite(finalEventsEnd, snakemake@output$events_end)
