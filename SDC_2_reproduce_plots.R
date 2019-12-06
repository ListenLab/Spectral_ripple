# Create the data and figures for the Spectral Ripple paper
# Packages used by this script:
library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr)
library(grid)
library(gridExtra)
library(purrr)

# Remove current objects in the workspace
rm(list = ls())

plots_folder <- "C:\\Name\\Your\\Preferred\\Directory\\Here"

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Un-comment and run only one of the following lines
# processor <- "Cochlear"
processor <- "AB_mid_scala"
# processor <- "AB_slim_J"
# processor <- "MedEl"
#=================================#
# Then, run the rest of the script. 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# set the frequency-channel allocation and other details
# for the CI processors
if (processor == "AB_mid_scala"){
  channel_boundaries.ci <- c(250, 416, 494, 587, 697, 828, 983, 1168, 1387, 
                             1648, 1958,2326, 2763, 3281, 3898, 4630, 8700)
  manufacturer_name <- "Advanced Bionics Mid-Scala"
  
  num_channels <- 16
  
  electrode_spacing = 0.975
  array_length <- 15
}

if (processor == "AB_slim_J"){
  channel_boundaries.ci <- c(250, 416, 494, 587, 697, 828, 983, 1168, 1387, 
                             1648, 1958,2326, 2763, 3281, 3898, 4630, 8700)
  manufacturer_name <- "Advanced Bionics Slim J"
  
  num_channels <- 16
  
  electrode_spacing = 1.3
  array_length <- 20
  
}

if (processor == "Cochlear"){
  channel_boundaries.ci <- c(188,313,438, 563, 688, 813,938,
                             1063, 1188, 1313, 1563, 1813,
                             2063, 2313, 2688, 3063, 3563,
                             4063, 4688, 5313, 6063, 6938, 7938)
  
  num_channels <- 22
  
  electrode_spacing = 0.75
  
  manufacturer_name <- "Cochlear"
}
if (processor == "MedEl"){
  # We only have the center frequencies of the device
  # not frequency bands, like for the other devices. 
  # this is because the Med-El device uses bell-shaped filters
  # intead of FFT processing. 
  # So by calculating these filters for the Med-El device,
  # we are doing two important things:
  # 1) violating the shape of the filters
  # 2) doing some fancy calculations to derive the intermediate 
  #   frequency boundaries between channels,
  #   if they were to exist (they dont; see previous note)
  
  center_freqs <- c(125, 234, 385, 582, 840, 1182, 1631, 2227, 3064, 4065, 5658, 7352)
  
  freq_to_mm <- function(freq=NULL, A = 165.4, a = 2.1, length = 35, k = 0.88){
    mm <- log10((freq/A)+k)*length/a
    return(mm)
  }
  mm_to_freq <- function(position=NULL, A = 165.4, a = 2.1, length = 35, k = 0.88){
    freq <- A*((10^(a*position/length))-k)
    return(freq)
  }
  places <- freq_to_mm(center_freqs)
  intermediate_place_boundaries <- 
    c(
      # border of first channel
      freq_to_mm(100),
      # plce-derived boundaries between center frequencies of channels 2-11
      # center freq plus...
      freq_to_mm(center_freqs)[1:11] + 
        # half the distance between the channel and the next
        0.5*diff(places)
    )
  
  channel_boundaries.ci <- c(mm_to_freq(intermediate_place_boundaries), 8500) %>%
    round(., 1)
  # starts at 100, ends at 8500
  
  # ... and that's how we got the following numbers:
  channel_boundaries.ci <- c(100, 175,  303, 476, 701,
                             998, 1390, 1907, 2614,
                             3531, 4798, 6451, 8500)
  
  manufacturer_name <- "Med-El"
  num_channels <- 12
  electrode_spacing = 2.4
  # 12 electrodes over 26.4mm (26.4/11 = 2.4)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
freq_to_mm <- function(freq=NULL, A = 165.4, a = 2.1, length = 35, k = 0.88){
  mm <- log10((freq/A)+k)*length/a
  return(mm)
}
mm_to_freq <- function(position=NULL, A = 165.4, a = 2.1, length = 35, k = 0.88){
  freq <- A*((10^(a*position/length))-k)
  return(freq)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Turn all the values into a data frame for easy handling. 
channel_boundaries_df <- data.frame(Frequency = channel_boundaries.ci)

# Make data frame of channel frequency lower and upper limits
chan <- data.frame(lower = c(channel_boundaries.ci[1:(length(channel_boundaries.ci)-1)]),
                   upper = c(channel_boundaries.ci[2:(length(channel_boundaries.ci))])) %>%
  mutate(channel = 1:n())

# Get center frequencies
chan$center <- floor(rowMeans(chan[,c("lower","upper")]))

# Code channel bandwidths
chan$oct_diff <- with(chan, (log(upper / lower))/ log(2))
chan$linear_bw <- with(chan, upper - lower)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Calculate maximum ripple density 
# that would fit into the narrowest filter
# before full cycle
0.5 * 1/min(chan$oct_diff)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Shift in insertion depth for most apical electrode
# in an ideal case, it would be this:
distance_from_apex <- freq_to_mm(chan$center[1])
# but you can override it here, to create a spectral shift
distance_from_apex <- 10

# Indicate where the electrodes will fall in the cochlea
electrode_positions <- 
  seq(distance_from_apex, distance_from_apex+((num_channels-1)*electrode_spacing), by = electrode_spacing)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#                         FUNCTIONS                              #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
make_ripple_spec <- function(RPO, phase = 0){
  # Create idealized log-spaced ripple from scratch
  # returns a data frame
  # with the specified number of RPO (the only input argument)
  # along with frequency in Hz
  #-----------------------------------------#
  # Fixed parameters:
  # RPO = 1 # This is an input argument
  mod_depth = 1
  log2_freq_offset = log(100,2) - 1
  log2_freq_offset = log(10,2) - 1
  ripples_per_octave = RPO

  # Frequency range (log spaced on the front end)
  log_freq_start = log(2^5, 2)
  log_freq_end = log(2^15, 2)  
  
  num_frequency_samples <- 4096*2
  
  # The vector of frequencies to sample the spectrum
  Frequencies <- 2^(seq(log_freq_start, log_freq_end, length.out = num_frequency_samples))
  Log_Frequencies <- log(Frequencies, 2)
  
  #=============================================================#
  # there are 6 octaves here (128-256-512-1024-2048-4096, 8192)
  # and there are 4096 samples
  # so there are 4096 / 6 (483) samples per octave
  log_freq_sample_density <- (log_freq_end - log_freq_start) / num_frequency_samples
  log_freq_samplerate = 1/log_freq_sample_density
  # This is another way of calculating the log freq samplerate;
  # it gives the same number as the line above
  log_freq_samplerate2 <- 1/diff(Log_Frequencies)[1]
  # 'log_freq_samplerate' should be the number of frequencies sampled for every octave
  
  #--- Make the Ripple Pattern Here ---#
  ripple <- (1-(mod_depth/2))+
    (sin((((log(Frequencies,2))-log2_freq_offset)*2*pi*ripples_per_octave)+
           ((phase/360)*(2*pi))))*
    (0.5*mod_depth)
  
  # Make a data frame of the rippled amplitude spectrum
  df_ripple <- data.frame(Frequency = Frequencies,
                          Power = ripple,
                          num_frequency_samples = num_frequency_samples,
                          log_freq_samplerate = log_freq_samplerate,
                          RPO = RPO,
                          phase = phase)
  
  # Output the data frame 
  # that contains all the info needed for plotting & analysis
  return(df_ripple)
}

send_through_CI_processor <- function(ripple_spec){
  # function that takes a spectrum 
  # (data frame with Frequency & Power)
  # which is the output of `make_ripple_spec()`
  # and sends it through a simulation of a 
  # Cochlear Implant speech processor
  # and returns the same spectrum (data frame)
  # with original Power renamed as as `Power_ideal`
  # and a new column with simulated output as `Power`
  # to refer to the "electric" spectrum
  spectral_modulation_depth <- 30
  
  ripple_spec %>% 
    # assign channel number
    #(0 for freqs below chan 1)
    mutate(channel = findInterval(Frequency, channel_boundaries.ci)) %>%
    mutate(channel = ifelse(channel %in% 1:num_channels, channel, NA)) %>%
    # Rename the original power column to be Power_ideal
    mutate(Power_ideal = Power) %>%
    # for each channel... 
    group_by(channel) %>%
    # indicate where the electrode would fall in the cochlea
    mutate(electrode_position_mm = electrode_positions[channel]) %>%
    # average all the power that falls within that channel
    mutate(Power = mean(Power)) %>%
    # exit by-channel grouping
    ungroup() %>%
    # instead of having two side-by-side columns 
    # of ideal acoustic power and simulated electric power,
    # stack them so that it's long-format data,
    # with the designation of acoustic/electric indicated on a separate column
    # (this make it a lot easier for plotting in R)
    melt(., measure.vars = c("Power","Power_ideal"), value.name = "Power") %>%
    # make sure Acoustic/Electric isn't a Factor
    mutate(variable = as.character(variable)) %>%
    # Rename designations to be more transparent
    within(., variable[variable=="Power_ideal"] <- "Acoustic") %>%
    within(., variable[variable=="Power"] <- "Electric") %>%
    # `Stimulus` indicated Acoustic / Electric
    rename(Stimulus = variable) %>%
    mutate(Power_dB = (Power * spectral_modulation_depth) + intensity_offset) %>%
    # return the whole data frame. 
    return()
}

# Wrap both functions in one
make_CI_spectrum <- function(RPO, phase) {
  make_ripple_spec(RPO, phase) %>% 
    send_through_CI_processor() %>%
    select(-num_frequency_samples, -log_freq_samplerate) %>%
    return()
}

# Take CI ripple spectrum, 
# simulate activation pattern in cochlea
activate_in_cochlea <- function(CI_ripple_spectrum, 
                                maintain_ideal = FALSE,
                                peak_pick = FALSE){
  # for each channel, find the cochlear position
  # and the level of activation
  # given the shift in mm (distance from apex)
  #---------------------------------------------#
  # maintain_ideal: don't simulate actual electrode position;
  # just treat the electrodes as if they were 
  # tonotopically perfect
  #---------------------------------------------#
  # CI_ripple_spectrum <- make_CI_spectrum(1) %>%
  #   dplyr::filter(!is.na(channel)) %>%
  #   dplyr::filter(Stimulus == "Electric")
  # CI_ripple_spectrum <- make_CI_spectrum(RPO, phase) %>%
  #   group_by(RPO, phase)
  
  electrode_positions <- CI_ripple_spectrum %>%
    pull(electrode_position_mm) %>%
    unique() %>%
    na.omit() %>%
    as.numeric()
  
  ideal_target_frequencies <- chan$center
  
  cochlear_activation_sites <- CI_ripple_spectrum %>%
    dplyr::filter(!is.na(channel)) %>%
    dplyr::filter(Stimulus == "Electric") %>%
    mutate(Frequency = mm_to_freq(electrode_position_mm)) %>%
    select(RPO, phase, Stimulus, channel, 
           electrode_position_mm, Frequency, Power, Power_dB) %>%
    unique() 
  
  if (maintain_ideal == TRUE){
    cochlear_activation_sites <- CI_ripple_spectrum %>%
      dplyr::filter(!is.na(channel)) %>%
      dplyr::filter(Stimulus == "Electric") %>%
      mutate(Frequency = mm_to_freq(electrode_position_mm)) %>%
      select(RPO, phase, Stimulus, channel, electrode_position_mm, Frequency, Power, Power_dB) %>%
      unique() %>%
      mutate(Frequency = ideal_target_frequencies[channel])  %>%
      mutate(electrode_position_mm = freq_to_mm(Frequency))
    
  }
  
  if (peak_pick){
    # only keep the top 8 peaks
    cochlear_activation_sites <- 
      cochlear_activation_sites %>%
      arrange(desc(Power_dB)) %>%
      within(., Power_dB[9:nrow(.)] <- -60) %>%
      arrange(RPO, channel)
  }
  
  return(cochlear_activation_sites)
}

spread_db_per_mm <- function(cf, dB, rolloff_per_mm){
  # simulate_spread of activation
  log_freq_start = log(2^5, 2)
  log_freq_end = log(2^15, 2)  
  
  num_frequency_samples <- 4096*2
  
  # The vector of frequencies to sample the spectrum
  Frequencies <- 2^(seq(log_freq_start, log_freq_end, length.out = num_frequency_samples))
  Log_Frequencies <- log(Frequencies, 2)
  
  aA = 165.4
  a = 2.1
  cochlear_length = 35
  k = 0.88
  
  mm_positions <- freq_to_mm(Frequencies)
  dB_adjustment <- -abs(mm_positions - freq_to_mm(cf))*rolloff_per_mm
  dB_output <- dB + dB_adjustment
  
  output <- data.frame(rolloff = rolloff_per_mm, 
                       cf = cf, 
                       Frequency = Frequencies, 
                       dB_adjustment,
                       Power_dB = dB_output)
  return(output)
}

simulate_spread_one_channel <- function(activation_points, dB_rolloff_per_mm = 3.33){
  # takes a single row from the data frame output 
  # from the 'activate in cochlea' function,
  # which contains center_freq and dB
  
  output <- activation_points %>%
    do(., spread_db_per_mm(cf = .$Frequency,
                           dB = .$Power_dB, 
                           rolloff_per_mm = dB_rolloff_per_mm)) %>%
    mutate(Stimulus = "Electric") %>%
    mutate(electrode_position_mm = freq_to_mm(cf)) %>%
    mutate(cochlear_position = freq_to_mm(Frequency)) 
  
  return(output)
}

simulate_spread_all_channels <- function(cochlear_activation_spots, dB_rolloff_per_mm, mod_depth){
  cochlear_activation_spots %>%
    group_by(RPO, phase, channel, electrode_position_mm) %>%
    do(., simulate_spread_one_channel(., dB_rolloff_per_mm = dB_rolloff_per_mm)) %>%
    select(-cf) %>%
    bind_rows() %>%
    ungroup() %>%
    mutate(Power = Power_dB / mod_depth) %>%
    return()
}

pressure_Pa_to_dB <- function(pressure_Pa){
  reference_level <- 0.00002
  dB <- 20 * log10(pressure_Pa / reference_level)
  return(dB)
}

dB_to_pressure_Pa <- function(dB){
  pressure_Pa <- 0.00002 * 10^(dB/20)
  return(pressure_Pa)
}

`%dB+%` <- function(lhs, rhs){
  # function for adding decibel values 
  # # example:
  # 60 %dB+% 60
  intermediate_sum <- dB_to_pressure_Pa(lhs) + dB_to_pressure_Pa(rhs)
  output <- pressure_Pa_to_dB(intermediate_sum)
  return(output)
}

sum_spread_channels <- function(electric_spectral_channels, peak_dB){
  #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
  # add up the power at each position,
  # making sure to add dB values 
  # using the special %dB+% function
  #----------------------------------#
  # electric_spectral_channels <- activation_each_channel
  apical_frequency_boundary <- electric_spectral_channels$apical_frequency_boundary[1]
  basal_frequency_boundary <- electric_spectral_channels$basal_frequency_boundary[1]
  
  electric_spectral_channels_summed <- electric_spectral_channels %>%
    group_by(RPO, phase, Frequency) %>%
    # Reduce function adds up ALL the dB values at that frequency
    summarise(Power_dB = Reduce(`%dB+%`, Power_dB)) %>%
    ungroup()
  
  temp_peak_dB <- max(electric_spectral_channels_summed$Power_dB)
  dB_adjustment <- peak_dB - temp_peak_dB
  
  electric_spectral_channels_summed$Power_dB <- 
    electric_spectral_channels_summed$Power_dB + dB_adjustment
  
  # Find intensity values at spectral edges
  apical_edge_Power_dB <- 
    electric_spectral_channels_summed$Power_dB[
      which.min(abs(electric_spectral_channels_summed$Frequency - apical_frequency_boundary))]
  
  basal_edge_Power_dB <- 
    electric_spectral_channels_summed$Power_dB[
      which.min(abs(electric_spectral_channels_summed$Frequency - basal_frequency_boundary))]
  
  edge_dB_min <- min(c(apical_edge_Power_dB, basal_edge_Power_dB))
  
  # Copy into new value column
  electric_spectral_channels_summed$Power_dB_cut <- electric_spectral_channels_summed$Power_dB
  
  # Walk down both sides,
  # assign the minimum edge dB value to extend down the sidebands
  electric_spectral_channels_summed$Power_dB_cut[
    (electric_spectral_channels_summed$Frequency < apical_frequency_boundary & 
       electric_spectral_channels_summed$Power_dB_cut < edge_dB_min)] <- edge_dB_min
  
  electric_spectral_channels_summed$Power_dB_cut[
    (electric_spectral_channels_summed$Frequency > basal_frequency_boundary & 
       electric_spectral_channels_summed$Power_dB_cut < edge_dB_min)] <- edge_dB_min
  
  # DC offset
  electric_spectral_channels_summed$Power_dB_cut_AC <- 
    electric_spectral_channels_summed$Power_dB_cut - edge_dB_min
  
  return(electric_spectral_channels_summed)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# RIPPLE PARAMETERS FOR FIGURE 1
RPO <- 1
dB_rolloff_per_mm <- 9.33
phase <- 225
spectral_modulation_depth <- 30
intensity_offset <- 0
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
make_channel_activation_pattern <- function(RPO, phase, maintain_ideal, peak_pick = FALSE){
  # simulate spread of activation
  
  output <- make_CI_spectrum(RPO, phase) %>% 
    group_by(RPO, phase) %>%
    activate_in_cochlea(maintain_ideal = maintain_ideal, peak_pick = peak_pick) %>%
    simulate_spread_all_channels(dB_rolloff_per_mm = dB_rolloff_per_mm, 
                                 mod_depth = spectral_modulation_depth) 
  
  # Here's where we need to do a conditional apical & basal site,
  # in case of peak picking 
  if (peak_pick == FALSE){
    apical_frequency_boundary <- mm_to_freq(min(output$electrode_position_mm))
    basal_frequency_boundary <- mm_to_freq(max(output$electrode_position_mm)) 
  }
  if (peak_pick == TRUE){
    
    picked_electrode_positions <- 
      make_CI_spectrum(RPO, phase) %>% 
      group_by(RPO, phase) %>%
      activate_in_cochlea(maintain_ideal = maintain_ideal, peak_pick = peak_pick) %>%
      dplyr::filter(Power_dB > -60) %>%
      pull(electrode_position_mm)
    
    apical_frequency_boundary <- mm_to_freq(min(picked_electrode_positions))
    basal_frequency_boundary <- mm_to_freq(max(picked_electrode_positions)) 
  }
  
  output$apical_frequency_boundary <- apical_frequency_boundary
  output$basal_frequency_boundary <- basal_frequency_boundary
  
  return(output)
}
# # examples:
# make_channel_activation_pattern(1, 225, TRUE) %>%
#   ggplot()+
#   aes(x = Frequency, y = Power_dB, group = channel)+
#   geom_line()+
#   scale_x_log10(breaks = octave_freq_breaks)

# ggplot(df_activation_pattern_channels_ideal)+
#   aes(x = Frequency, y = Power_dB, group = channel)+
#   geom_line()+
#   geom_vline(xintercept = c(500, 1000, 2000, 4000, 8000))+
#   scale_x_log10(breaks = c(500, 1000, 2000, 4000, 8000))

mean_Power_dB_orig <- 
  make_CI_spectrum(RPO, phase) %>% 
  activate_in_cochlea() %>%
  ungroup() %>%
  pull(Power_dB) %>%
  mean()

make_summed_activation_pattern <- function(RPO, phase, maintain_ideal, peak_pick = FALSE){
  # Function to create summed activation pattern
  # across all channels & cochlear locations
  activation_each_channel <- 
    make_channel_activation_pattern(RPO, phase, 
                                    maintain_ideal = maintain_ideal, 
                                    peak_pick = peak_pick)
  
  peak_dB = max(activation_each_channel$Power_dB)
  
  activation_each_channel %>%
    sum_spread_channels(peak_dB = peak_dB) %>%
    return()
}

# # Spots where the electrodes would land, and their power
df_cochlear_activation_spots <- make_CI_spectrum(RPO, phase = phase) %>%
  activate_in_cochlea() %>%
  mutate(cochlear_position = freq_to_mm(Frequency))

# Where the channels are analyzed in the spectrum,
analysis_freqs <- chan %>%
  mutate(type = "Analysis") %>%
  mutate(cochlear_position = freq_to_mm(center)) %>%
  select(type, channel, Frequency = center, cochlear_position) %>%
  unique()

# Where the channels are delivered in the cochlea
carrier_places <- df_cochlear_activation_spots %>%
  mutate(type = "Carrier") %>%
  mutate(cochlear_position = freq_to_mm(Frequency)) %>%
  select(type, channel, Frequency, cochlear_position) %>%
  unique()

# Data frame with that transformation
df_frequency_shift <- bind_rows(analysis_freqs, carrier_places)
rm(analysis_freqs, carrier_places)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
make_spectra_comparison <- function(RPO, phase, full_freq_sampling = FALSE){
  # wrapper function to generate all the spectra that are needed for the figure
  
  message(paste0("Making spectra comparison for ",RPO,"  RPO, ", phase," phase"))
  
  spectrum_acoustic_ideal <- make_CI_spectrum(RPO, phase) %>% 
    dplyr::filter(Stimulus == "Acoustic") %>%
    mutate(Spectrum = "Acoustic ideal") 
  
  spectrum_electric_ideal <- make_CI_spectrum(RPO, phase) %>% 
    dplyr::filter(Stimulus == "Electric") %>%
    mutate(Spectrum = "Electric discrete")
  
  spectrum_electric_summed <- 
    make_summed_activation_pattern(RPO, phase, maintain_ideal = FALSE) %>%
    mutate(Spectrum = "Summed activation")

  #----------------------------------------------------------------#
  peak_pick_selection <-ifelse(processor== "Cochlear", TRUE, FALSE)
  
  spectrum_electric_summed_picked <- 
    make_summed_activation_pattern(RPO, phase, maintain_ideal = FALSE, peak_pick = peak_pick_selection) %>%
    mutate(Spectrum = "Summed activation picked")
  
  spectrum_electric_summed_ideal <- 
    make_summed_activation_pattern(RPO, phase, maintain_ideal = TRUE) %>%
    mutate(Spectrum = "Summed activation ideal")
  
  df_spectra <- spectrum_acoustic_ideal %>%
    select(RPO, phase, Frequency, Power_dB, Spectrum) %>%
    mutate(Power_dB_cut = Power_dB) %>%
    bind_rows(., 
              {spectrum_electric_ideal %>%
                  select(RPO, phase, Frequency, Power_dB, Spectrum) %>%
                  mutate(Power_dB_cut = Power_dB)}) %>%
    bind_rows(., 
              {spectrum_electric_summed %>%
                  select(RPO, phase, Frequency, Power_dB, Power_dB_cut, Spectrum) }) %>%
    bind_rows(., 
              {spectrum_electric_summed_picked %>%
                  select(RPO, phase, Frequency, Power_dB, Power_dB_cut, Spectrum) }) %>%
    bind_rows(., 
              {spectrum_electric_summed_ideal %>%
                  select(RPO, phase, Frequency, Power_dB, Power_dB_cut, Spectrum)})
  
  # Remove data beyond reasonable frequency limits
  df_spectra[
    (df_spectra$Spectrum == "Electric discrete" & 
       df_spectra$Frequency < min(channel_boundaries.ci)), "Power_dB"] <- NA
  df_spectra[
    (df_spectra$Spectrum == "Electric discrete" & 
       df_spectra$Frequency > max(channel_boundaries.ci)), "Power_dB"] <- NA
  
  # Order the spectra in plotting sequence
  df_spectra$Spectrum <- 
    factor(df_spectra$Spectrum,
           levels = c("Acoustic ideal","Electric discrete",
                      "Summed activation ideal",
                      "Summed activation",
                      "Summed activation picked"))
  
  # Make prettier labels for plotting
  df_spectra$Spectrum_label <- df_spectra$Spectrum
  ordered_spectrum_labels <- c("Acoustic\n(ideal)",
                               "Electric\n(topnotopic, discrete)",
                               "Summed\nelectric activation\n(tonotopic)",
                               "Summed\nelectric activation",
                               "Summed\nelectric activation\n(peak-picked)")
  levels(df_spectra$Spectrum_label) <- ordered_spectrum_labels
  
  return(df_spectra)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Generate the big data frame for the first plot
dB_rolloff_per_mm <- 9.3
RPOs_to_plot <- c(0.5, 1.5, 4)
phases_to_plot <- rep(225, 3)
df_spectra <- map2(.x = RPOs_to_plot, 
                   .y = phases_to_plot, 
                   .f = make_spectra_comparison) %>%
  bind_rows()

# Idealized (discrete) electrical activation
channel_activation_patterns_ideal <- 
  map2(.x = RPOs_to_plot,
       .y = phases_to_plot, 
       .f = make_channel_activation_pattern, 
       maintain_ideal = TRUE) %>%
  bind_rows() %>%
  mutate(
    Spectrum = "Summed activation ideal",
    Spectrum_label = "Summed\nelectric activation\n(tonotopic)")

# Simulated spread of excitation 
channel_activation_patterns <- 
  map2(.x = RPOs_to_plot,
       .y = phases_to_plot, 
       .f = make_channel_activation_pattern, 
       maintain_ideal = FALSE) %>%
  bind_rows() %>%
  mutate(
    Spectrum = "Summed activation",
    Spectrum_label = "Summed\nelectric activation")

# Now with peak-picking
channel_activation_patterns_picked <- 
  map2(.x = RPOs_to_plot,
       .y = phases_to_plot, 
       .f = make_channel_activation_pattern, 
       maintain_ideal = FALSE, peak_pick = TRUE) %>%
  bind_rows() %>%
  mutate(
    Spectrum = "Summed activation picked",
    Spectrum_label = "Summed\nelectric activation\n(peak-picked)")

# Put those data frames together
channel_activation_patterns_all <- 
  bind_rows(channel_activation_patterns,
            channel_activation_patterns_picked,
            channel_activation_patterns_ideal) %>%
  mutate(Spectrum = factor(Spectrum,
                           levels = levels(df_spectra$Spectrum)))

# Clean up intermediate objects
rm(channel_activation_patterns,
    channel_activation_patterns_picked,
    channel_activation_patterns_ideal)

# only maintain peak-picked for Cochlear devices
if (processor != "Cochlear"){
  channel_activation_patterns_all <- channel_activation_patterns_all %>%
    dplyr::filter(Spectrum != "Summed activation picked")
  
  df_spectra <- df_spectra %>%
    dplyr::filter(Spectrum != "Summed activation picked")
}

# Declare an order for presentations
ordered_spectrum_types <- 
  c("Acoustic\n(ideal)",
    "Electric\n(topnotopic, discrete)",
    "Summed\nelectric activation\n(tonotopic)",
    "Summed\nelectric activation\n(on CI array)",
    "Summed\nelectric activation\n(on CI array)\n(peak-picked)")
# levels(df_spectra$Spectrum_label) <- ordered_spectrum_types
channel_activation_patterns_all$Spectrum_label <- channel_activation_patterns_all$Spectrum

# Assign those re-ordered labels
levels(channel_activation_patterns_all$Spectrum_label) <- ordered_spectrum_types
levels(df_spectra$Spectrum_label) <- ordered_spectrum_types

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Plot aesthetics
octave_freq_breaks <- c(125, 250, 500, 1000, 2000, 4000, 8000)

electric_blue <- "#83C2F8"

spectrum_colors <- c(`Acoustic ideal` = "black",
                     `Electric discrete` = electric_blue,
                     `Summed activation` = "#6B0F1A",
                     `Summed activation picked` = "#28536B",
                     `Summed activation ideal` = "gray50")

df_spectra$RPO_label <- paste0(df_spectra$RPO, " RPO")
channel_activation_patterns_all$RPO_label <- paste0(channel_activation_patterns_all$RPO, " RPO")

# Check the size of the data frame object on disc
format(object.size(channel_activation_patterns_all), units = "Mb")

# reduce size of channel_activation_patterns_all
channel_activation_patterns_all_reduced <- 
  channel_activation_patterns_all[
    (channel_activation_patterns_all$Frequency <= 12000 &
       channel_activation_patterns_all$Frequency >= 80),]

sampled_frequencies <- 
  channel_activation_patterns_all$Frequency[
    (channel_activation_patterns_all$Frequency <= 12000 &
       channel_activation_patterns_all$Frequency >= 80)] %>%
  unique() %>%
  sort()

reduced_frequency_indices <- (length(sampled_frequencies) + 1) - 
  (10^seq(0, log10(length(sampled_frequencies)), length.out = 2000)) %>%
  round(., 0) %>%
  unique() %>%
  rev()

sampled_frequencies_reduced <- sampled_frequencies[reduced_frequency_indices]

channel_activation_patterns_all_to_plot <-
  channel_activation_patterns_all[
    channel_activation_patterns_all$Frequency %in% sampled_frequencies_reduced,]

# verify compression of data
format(object.size(channel_activation_patterns_all_to_plot), units = "Mb")

# remove intermediate object
rm(channel_activation_patterns_all)

df_spectra_to_plot <-
  df_spectra[
    df_spectra$Frequency %in% sampled_frequencies_reduced,]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
px_Figure_1_example_spectra <- ggplot(df_spectra_to_plot)+
  aes(x = Frequency, y = Power_dB, color = Spectrum)+
  annotate("rect", xmin = 1, xmax = min(channel_boundaries.ci),
           ymin = -spectral_modulation_depth, ymax = spectral_modulation_depth*2,
           fill = "gray80", alpha = 0.5)+
  annotate("rect", xmin = max(channel_boundaries.ci), xmax = 30000,
           ymin = -spectral_modulation_depth, ymax = spectral_modulation_depth*2,
           fill = "gray80", alpha = 0.5)+
  scale_color_manual(values = spectrum_colors)+
  geom_line(size = 1.2)+
  geom_line(size = 1, aes(y = Power_dB_cut), linetype = "longdash")+
  xlab("Frequency (Hz)")+
  ylab("Relative dB")+
  scale_x_log10(breaks = octave_freq_breaks)+
  coord_cartesian(xlim = c(110, 9500),
                  ylim = c(0, 33))+
  theme_bw()+
  theme(legend.position = "none")+
  theme(panel.grid.minor = element_blank())+
  theme(panel.grid.major = element_blank())+
  theme(strip.text.y = element_text(angle = 0))+
  facet_grid(Spectrum_label ~ RPO_label)+
  ggtitle(processor)
# px_Figure_1_example_spectra

# now add the individual channels underneath
px_Figure_1_example_spectra_w_channels <- 
  ggplot(df_spectra_to_plot)+
  aes(x = Frequency, y = Power_dB, color = Spectrum)+
  annotate("rect", xmin = 1, xmax = min(channel_boundaries.ci),
           ymin = -spectral_modulation_depth, ymax = spectral_modulation_depth*2,
           fill = "gray80", alpha = 0.5)+
  annotate("rect", xmin = max(channel_boundaries.ci), xmax = 30000,
           ymin = -spectral_modulation_depth, ymax = spectral_modulation_depth*2,
           fill = "gray80", alpha = 0.5)+
  scale_color_manual(values = spectrum_colors)+
  # Individual channel activation patterns
  geom_ribbon(data = channel_activation_patterns_all_to_plot,
              aes(ymin = -Inf, ymax = Power_dB,
                  x = Frequency, group = channel),
              fill = "black", alpha = 0.3, color = NA)+
  # Overall summed activation pattern
  geom_line(size = 1.2)+
  # summed activation pattern, 
  # with disebands cut away
  geom_line(size = 0.9, aes(y = Power_dB_cut), linetype = "longdash")+
  xlab("Frequency (Hz)")+
  ylab("Relative dB")+
  scale_x_log10(breaks = octave_freq_breaks)+
  coord_cartesian(xlim = c(110, 9500),
                  ylim = c(0, 33))+
  theme_bw()+
  theme(legend.position = "none")+
  theme(panel.grid.minor = element_line(size = 0.4))+
  theme(panel.grid.major = element_blank())+
  theme(strip.text.y = element_text(angle = 0))+
  facet_grid(Spectrum_label ~ RPO_label)+
  ggtitle(processor)
px_Figure_1_example_spectra_w_channels

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Saving the plot output
plot_1_output_height = ifelse(processor == "Cochlear", 4.25, 3.6)
plot_1_output_width = 9.35

setwd(plots_folder)
ggsave(px_Figure_1_example_spectra_w_channels,
       file = paste0("Figure_1_example_spectra_channels_",processor,".png"),
       height = plot_1_output_height, width = plot_1_output_width, dpi = 300)
ggsave(px_Figure_1_example_spectra_w_channels,
       file = paste0("Figure_1_example_spectra_channels_",processor,".pdf"),
       height = plot_1_output_height, width = plot_1_output_width, device = cairo_pdf)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#                     FIGURE 2                          #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# MANY RPOs for a spectral modulation density analysis
RPOs_to_analyze <- seq(0.2, 7.5, 0.1)

# All the same phase (225 degrees, arbitrary)
phases_to_analyze <- rep(225, length(RPOs_to_analyze))

# Make spectra for each of those RPOs at that phase
df_spectra_across_RPOs <- map2(.x = RPOs_to_analyze, 
                               .y = phases_to_analyze, 
                               .f = make_spectra_comparison,
                               full_freq_sampling = TRUE) %>%
  bind_rows()

# Only maintain peak-picked pattern for Cochlear devices
if (processor != "Cochlear"){
  df_spectra_across_RPOs <- df_spectra_across_RPOs %>%
    dplyr::filter(Spectrum != "Summed activation picked")
}

# Analyze log frequency components in log-spaced ripple
analyze_ripple <- function(df_ripple, HPF = 0.2){
  # df_ripple (input) should have two columns: Power_dB and Frequency
  # Also, this function relies on some other variables that were delcared earlier
  # HPF is a high-pass filter that removes super-low-density modulations
  
  # Extract the power vector (the amplitude spectrum)
  ripple <- df_ripple$Power_dB_cut
  ripple[is.na(ripple)] <- min(ripple, na.rm = TRUE)
  # fft on the log-frequency-sampled ripple
  # make FFT, take absolute value (power spectrum)
  # then take first element up to the element beyond the halfway point
  # (common practice to get power spectrum, 
  # as it's symmetrical)
  ripple_fft <- fft(ripple) %>% abs %>% `[`(., 1:(length(.)/2+1))

  #*******************************************************************# ??????
  # Some parameters used to create the original ripples:
  log_freq_start = log(2^5, 2)
  log_freq_end = log(2^15, 2)  
  
  num_frequency_samples <- 4096*2
  
  # The vector of frequencies to sample the spectrum
  Frequencies <- 2^(seq(log_freq_start, log_freq_end, length.out = num_frequency_samples))
  Log_Frequencies <- log(Frequencies, 2)
  
  #=============================================================#
  # There are 6 octaves here (128-256-512-1024-2048-4096, 8192)
  # and there are 4096 samples
  # so there are 4096 / 6 (483) samples per octave
  log_freq_sample_density <- (log_freq_end - log_freq_start) / num_frequency_samples
  log_freq_samplerate = 1/log_freq_sample_density

  log_freq_Nyquist <- log_freq_samplerate / 2
  num_frequency_samples <- 4096*2
  num_frequency_samples <- nrow(df_ripple)
  
  # Establish log-frequency smapling range
  log_freq_start = log(min(df_ripple$Frequency),2)
  log_freq_end = log(max(df_ripple$Frequency),2)
  log_freq_sample_density <- (log_freq_end - log_freq_start) / num_frequency_samples
  log_freq_samplerate = 1/log_freq_sample_density
  
  freq_bins <- seq(0,log_freq_Nyquist,by=log_freq_samplerate/num_frequency_samples); # frequency vector from 0 to the Nyquist
  #-----------------------------------------------------------------------------#
  # Convert to frequency space rather than bin index
  df_ripple_fft <- data.frame(Log_Mod_Frequency = freq_bins,
                              Power = ripple_fft) %>%
    mutate(bin = 1:n()) 
  
  # Obtain the frequency resolution (separation between frequency samples)
  # by taking the first difference between consecutive values
  df_ripple_log_freq_resolution <- diff(df_ripple_fft$Log_Mod_Frequency)[1]

  # Take only the first half, and exclude first point 
  # (Log_Mod_Frequency of 0, which is DC offset)
  df_ripple_fft <- df_ripple_fft[c(2:(nrow(df_ripple_fft)/2)),]
  
  # Exclude  electric offset
  df_ripple_fft$Power[df_ripple_fft$Log_Mod_Frequency < HPF] <- NA
  #-----------------------------------------------------------------------------#
  # Find peak (this is the estimated peak ripples per octave density that is represented in the spectrum)
  log_fft_peak <- df_ripple_fft[which.max(df_ripple_fft$Power), "Log_Mod_Frequency"]
  
  #-----------------------------------------------------------------------------#
  # plot the estimation of ripples per octave present in the original spectrum
  verify_LFMSpectrum <- ggplot(df_ripple_fft)+
    aes(x = Log_Mod_Frequency, y = Power)+
    geom_point()+
    geom_vline(xintercept = log_fft_peak)+
    geom_line()+
    coord_cartesian(xlim = c(0, 7))+
    scale_x_continuous(name = "Log spectral modulation frequency \n(Ripples per octave)",
                       breaks = 0:7)

  # add some indexical information to the data frame,
  # so that if you bind it with other DFs,
  # it maintains identity
  df_ripple_fft$RPO <- df_ripple$RPO[1]
  df_ripple_fft$peak <- log_fft_peak
  
  # Power of 61440 = 30 dB modulation depth
  # and Power scales linearly with dB,
  # as 30720 corresponds to 15 dB modulation depth
  # 61440 / 30
  # each dB is: 2048
  # which is a quarter of the sampling window
  df_ripple_fft$Power_dB <- df_ripple_fft$Power / (num_frequency_samples/4)
  # we also want to add the FFT-FFT of the idealized spectrum based on the RPO,
  # so we can compare it to the derived FFT-FFT
  # that will be accomplished not in THIS function,
  # but rather by running this function on ideal 
  # or non-ideal (processed) spectra. 
  
  return(df_ripple_fft)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Run the density analysis on all the RPOs
df_spectral_density <- df_spectra_across_RPOs %>%
  group_by(RPO, phase, Spectrum, Spectrum_label) %>%
  # analyze ripple density for each combination of 
  # RPO, phase, and spectrum type
  do(., analyze_ripple(., HPF = 0.2)) %>%
  dplyr::filter(Log_Mod_Frequency <= 8)

# Normalize modulation power within each RPO, phase, spectrum type
df_spectral_density <- df_spectral_density %>%
  group_by(RPO, phase, Spectrum, Spectrum_label) %>%
  mutate(Power_dB_norm = Power_dB / max(Power_dB, na.rm = TRUE))

# Cut away sideband frequency info for the acoustic sitmulus,
# just as was done for the electroc stimulus
min_freq <- channel_boundaries.ci[1]
max_freq <- channel_boundaries.ci[length(channel_boundaries.ci)]

# Add constant value at the sidebands
# equal to half the modulation depth
# (it's like zero-adding)
df_spectra_across_RPOs$Power_dB_cut[
  (df_spectra_across_RPOs$Spectrum == "Acoustic ideal" &
     df_spectra_across_RPOs$Frequency < min_freq)] <- spectral_modulation_depth/2

df_spectra_across_RPOs$Power_dB_cut[
  (df_spectra_across_RPOs$Spectrum == "Acoustic ideal" &
     df_spectra_across_RPOs$Frequency > max_freq)] <- spectral_modulation_depth/2

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Run the density analysis on all the RPOs
df_spectral_density <- df_spectra_across_RPOs %>%
  group_by(RPO, phase, Spectrum, Spectrum_label) %>%
  do(., analyze_ripple(., HPF = 0.2)) %>%
  # only keep densities below 8
  dplyr::filter(Log_Mod_Frequency <= 8)

# Normalize modulation power within each RPO
df_spectral_density <- df_spectral_density %>%
  group_by(RPO, phase, Spectrum, Spectrum_label) %>%
  mutate(Power_dB_norm = Power_dB / max(Power_dB, na.rm = TRUE))

# Produce the sinc function that represents the distortion 
# of modulation depth as a function of RPO
max_ideal_modulation_depth <- 
  df_spectral_density %>%
  dplyr::filter(Spectrum == "Acoustic ideal") %>%
  group_by(RPO) %>%
  summarise(max_depth = max(Power_dB, na.rm = TRUE)) %>%
  mutate(correction_factor = 30 / max_depth)

# Visualize that function
ggplot(max_ideal_modulation_depth)+
  aes(x = RPO, y = max_depth)+
  geom_point()

# ^ That's the intensity adjustment that needs to be made 
# simply on the basis of the desire to produce fair comparisons
df_spectral_density <- left_join(df_spectral_density, max_ideal_modulation_depth)

df_spectral_density$Power_dB_corrected <- 
  df_spectral_density$Power_dB * df_spectral_density$correction_factor

# labels for the simulations (for the plot)
ordered_spectrum_types <- 
  c("Acoustic\n(ideal)",
    "Electric\n(topnotopic, discrete)",
    "Summed\nelectric activation\n(tonotopic)",
    "Summed\nelectric activation\n(on CI array)",
    "Summed\nelectric activation\n(on CI array)\n(peak-picked)")

levels(df_spectral_density$Spectrum_label) <- ordered_spectrum_types

# Maximum modulation depth
# per spectrum type
df_spectral_density %>%
  dplyr::filter(Log_Mod_Frequency > 0.24) %>%
  group_by(Spectrum) %>%
  summarise(max_modulation = max(Power_dB_corrected))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# plot top of fill map as blue,
# fill 0 - 16 as white to black
px_RPO_IO_map <- 
  df_spectral_density %>%
  dplyr::filter(Log_Mod_Frequency > 0.24) %>%
  ggplot(.)+
  aes(x = RPO, y = Log_Mod_Frequency, fill = Power_dB_corrected)+
  geom_tile()+
  scale_fill_gradientn(colours = c("white", "black", "blue"),
                       guide = "legend", breaks = c(0, 6, 18, 30),
                       name = "Modulation\npower (dB)")+
  scale_x_continuous(breaks = 1:7)+
  scale_y_continuous(breaks = 1:7)+
  xlab("Input spectral modulation density (ripples per octave)")+
  ylab("Output\nspectral modulation density\n(ripples per octave)")+
  theme_bw()+
  theme(legend.justification=c(0, 1.6))+
  theme(legend.key.height = unit(1, "line"))+
  coord_fixed(ratio = 1, expand = FALSE)+
  facet_grid(. ~ Spectrum_label)
# px_RPO_IO_map

px_RPO_IO_map2_guide_line <- df_spectral_density %>%
  dplyr::filter(Log_Mod_Frequency > 0.24) %>%
  ggplot(.)+
  aes(x = RPO, y = Log_Mod_Frequency, fill = Power_dB_corrected)+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "black",
                      name = "Modulation\npower (dB)")+
  scale_x_continuous(breaks = 1:7)+
  scale_y_continuous(breaks = 1:7)+
  xlab("Input spectral modulation density (ripples per octave)")+
  ylab(paste0(processor,
              "\n\nOutput\nspectral modulation\ndensity\n(ripples per octave)"))+
  theme_bw()+
  theme(legend.justification=c(0, 1.6))+
  theme(legend.key.height = unit(1, "line"))+
  theme(axis.title.y = element_text(angle = 0))+
  coord_fixed(ratio = 1, expand = FALSE)+
  geom_abline(slope = 1, intercept = 0,
              size = 1.1,
              color = "blue",
              linetype = "dotted",
              alpha = 0.75)+
  facet_grid(. ~ Spectrum_label)
# px_RPO_IO_map2_guide_line

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Same thing, but for Normalized modulation depth
px_RPO_IO_map_norm <-
  df_spectral_density %>%
  dplyr::filter(Log_Mod_Frequency > 0.24) %>%
  ggplot(.)+
  aes(x = RPO, y = Log_Mod_Frequency, fill = Power_dB_norm)+
  geom_tile()+
  geom_abline(slope = 1, intercept = 0,
              size = 1.1,
              color = "blue",
              linetype = "dotted",
              alpha = 0.75)+
  scale_fill_gradient(low = "white", high = "black",
                      name = "Proportion of\nmodulation\npower\nwithin\ninput\nRPO", 
                      breaks = seq(0,1,0.25))+
  scale_x_continuous(breaks = 1:7)+
  scale_y_continuous(breaks = 1:7)+
  xlab("Input spectral modulation density (ripples per octave)")+
  ylab(paste0(processor,
              "\n\nOutput\nspectral modulation\ndensity\n(ripples per octave)"))+
  theme_bw()+
  theme(legend.justification=c(0, 1.5))+
  theme(legend.key.height = unit(1, "line"))+
  theme(axis.title.y = element_text(angle = 0))+
  coord_fixed(ratio = 1, expand = FALSE)+
  facet_grid(. ~ Spectrum_label)
# px_RPO_IO_map_norm

# Arrange multiple plots on one panel
px_Figure_2_ripple_density_IO <- 
  cowplot::plot_grid(px_RPO_IO_map2_guide_line,
                     px_RPO_IO_map_norm,
                     align = "v",axis = "l",
                     labels = "AUTO",
                     label_x = 0.05, hjust = 1, label_size=17,
                     ncol = 1, nrow = 2)
grid.newpage()
grid.draw(px_Figure_2_ripple_density_IO)

# remove sub-plot objects
rm(px_RPO_IO_map, px_RPO_IO_map_norm, px_RPO_IO_map2, px_RPO_IO_map2_guide_line)

# save the figure
width_fig_2 <- ifelse(processor == "Cochlear",9, 7.8)
ggsave(px_Figure_2_ripple_density_IO, 
       file = paste0("Figure_2_Ripple_density_IO_maps_",processor,".png"),
       height = 5.5, width = width_fig_2, dpi = 300)

ggsave(px_Figure_2_ripple_density_IO, 
       file = paste0("Figure_2_Ripple_density_IO_maps_",processor,".pdf"),
       height = 5.5, width = width_fig_2)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#                     FIGURE 3                          #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

df_spectral_density$RPO <- round(df_spectral_density$RPO, 1)
df_spectral_density$RPO_label <- paste(df_spectral_density$RPO, "\nRPO")

RPOs_to_demo <- c(0.5, 1.1, 1.7, 3.0, 4, 5)

df_density_OI_demo <- df_spectral_density %>%
  dplyr::filter(Log_Mod_Frequency > 0.24) %>%
  dplyr::filter(RPO %in% RPOs_to_demo) %>%
  dplyr::filter(Spectrum == "Summed activation") %>%
  # re-normalize
  group_by(RPO, RPO_label) %>%
  mutate(Power_dB_norm = Power_dB / max(Power_dB))

pick_peaks_in_mod_spectrum <- function(dsub){
  # roll along to grab sub-peaks
  # if Power_dB_norm is above 0.5
  # if sign of derivative changed from positive to negative
  # dsub <- df_density_OI_demo %>%
  #   dplyr::filter(RPO==4)
  
  dsub$diffs <- c(NA, diff(dsub$Power_dB_norm))
  dsub$sign <- sign(dsub$diffs)
  
  dsub$deflection <- NA
  for (row_index in 2:(nrow(dsub)-1)){
    dsub$deflection[row_index] <- 
      ifelse((dsub$sign[row_index] == 1 & dsub$sign[row_index+1] == -1), 1, 0)
  }
  # remove unneeded columns
  dsub$diffs <- NULL
  dsub$sign <- NULL
  
  return(dsub)
}

# Pick peaks in the spectral modulation spectrum 
# for each RPO 
df_density_OI_demo_peaks <- df_density_OI_demo %>%
  group_by(RPO, RPO_label) %>%
  do(., pick_peaks_in_mod_spectrum(.)) %>%
  dplyr::filter(Power_dB_norm == 1 | deflection == 1) %>%
  dplyr::filter(Power_dB_norm >= 0.5) %>%
  # assign labels to peaks and sub-peaks
  mutate(peak_type = ifelse(Power_dB_norm == 1, "Peak", "Sub-peak"))

RPO_references <- df_density_OI_demo %>%
  select(RPO, RPO_label) %>%
  unique()

max_power_dB <- max(df_density_OI_demo$Power_dB)

# plot extras
line_size = 1

RPO_annotation <- data.frame(RPO = 0.5, 
                             label = "Ripples\nper octave",
                             arrow_ypos = 0.5)

RPO_annotation$RPO_label <- paste(RPO_annotation$RPO, "\nRPO")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Produce the plot
px_mod_spectrum <- ggplot(df_density_OI_demo)+
  aes(x = Log_Mod_Frequency, y = Power_dB)+
  # add reference lines in the back
  geom_vline(data = RPO_references, aes(xintercept = RPO),
             color = "gray60", size = 2.3, alpha = 0.8)+
  geom_line(size = line_size)+
  geom_point(data = subset(df_density_OI_demo_peaks, peak_type == "Peak"),
             size = 1.9)+
  xlab("Output spectral modulation density\n(ripples per octave)")+
  ylab("Modulation power (dB)")+
  scale_x_continuous(breaks = 1:7)+
  scale_y_continuous(breaks = seq(0,8,2))+
  coord_cartesian(ylim = c(0, max_power_dB*1.1))+
  theme_bw()+
  theme(strip.text.y = element_text(angle = 0))+
  theme(panel.grid.minor = element_blank())+
  facet_grid(RPO_label ~ .)

px_mod_spectrum_norm <- px_mod_spectrum +
  aes(y = Power_dB_norm)+
  # white space above 1.0
  annotate("rect", fill = "white", 
           xmin = -Inf, xmax = Inf,
           ymin = 1, ymax = Inf)+
  geom_hline(yintercept = 1, color = "gray80")+
  geom_point(data = df_density_OI_demo_peaks,
             aes(color = peak_type, shape = peak_type),
             size = 1.9, fill = "white", stroke = 1.3)+
  scale_shape_manual(values = c(16, 21))+
  scale_color_manual(values = c("black", "firebrick"))+
  ylab("\n\nNormalized Spectral modulation power (proportion)")+
  scale_y_continuous(breaks = seq(0,1,0.5))+
  coord_cartesian(ylim = c(0, 1.05))+
  theme(panel.grid.minor = element_blank())+
  theme(legend.position = c(0.7, 0.65),
        legend.title = element_blank(), 
        legend.background = element_rect(color = "black"))
# px_mod_spectrum_norm

# Alter the axes so that the plots can fit nicely together
px_mod_spectrum_w_title <- px_mod_spectrum +
  theme(axis.title.y = element_blank())+
  ggtitle(paste0(processor,"\n\nSpectral modulation power\n(dB)"))+
  theme(title = element_text(size = 10))

px_mod_spectrum_norm_w_title <- px_mod_spectrum_norm +
  theme(axis.title.y = element_blank())+
  ggtitle("\nNormalized\nspectral modulation power\n(proportion)")+
  theme(title = element_text(size = 10))

px_mod_spectrum_norm_w_labels <- px_mod_spectrum_norm_w_title +
  geom_segment(data = RPO_annotation, inherit.aes = FALSE,
               aes(y = arrow_ypos, yend = arrow_ypos),
               x = 5.7, xend = 8.3, arrow = arrow(length = unit(0.4, "line")))+
  geom_label(data = RPO_annotation, inherit.aes = FALSE,
             aes(label = label), size = 2.7,
             x = 5.8, y = 0.5)+
  guides(shape = guide_legend(override.aes = list(size = 1.8)),
         color = guide_legend(override.aes = list(size = 1.8)))+
  theme(legend.position = c(0.7, 0.76))+
  theme(legend.text = element_text(size = 8),
        legend.key.size = unit(0.7, "line"),
        legend.spacing = unit(0.2, "line"))
# px_mod_spectrum_norm3

# Arrange the plots as one image
px_Figure_3 <- cowplot::plot_grid(px_mod_spectrum_w_title, 
                                  px_mod_spectrum_norm_w_labels,
                               ncol = 2)

rm(px_mod_spectrum, 
   px_mod_spectrum_norm_w_title, 
   px_mod_spectrum_norm, 
   px_mod_spectrum_norm_w_title,
   px_mod_spectrum_norm_w_labels)

# Save the plot
ggsave(px_Figure_3, file = paste0("Figure_3_IO_demos_",processor,".png"),
       height = 5.3, width = 5.8, dpi = 200)

ggsave(px_Figure_3, file = paste0("Figure_3_IO_demos_",processor,".pdf"),
       height = 5.3, width = 5.8, device = cairo_pdf)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#                     FIGURE 4                          #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

dB_rolloff_per_mm <- 9.3
RPOs_to_plot <- seq (0.5, 7, 0.5)
# add 180 to each of the numbers in sequence from 0 to 170 in increments of 10
phases_to_plot <- seq(0, 170, 10) %>% c(., . + 180)

# Fully cross all those combinations of parameters
parameters <- expand.grid(RPO = RPOs_to_plot,
                          phase = phases_to_plot)

# Make spectra for ALL of those parameters
# This step takes a LONG time. 
df_spectra_phases <- map2(.x = parameters$RPO, 
                   .y = parameters$phase, 
                   .f = make_spectra_comparison) %>%
  bind_rows()

# That's a pretty large object
format(object.size(df_spectra_phases), units = "Mb")

# only maintain peak-picked pattern for Cochlear devices
if (processor != "Cochlear"){
  df_spectra_phases <- df_spectra_phases %>%
    dplyr::filter(Spectrum != "Summed activation picked")
}

# Set the nice labels
df_spectra_phases$Spectrum_label <- df_spectra_phases$Spectrum
ordered_spectrum_labels <- c("Acoustic\n(ideal)",
                             "Electric\n(topnotopic, discrete)",
                             "Summed\nelectric activation\n(tonotopic)",
                             "Summed\nelectric activation",
                             "Summed\nelectric activation\n(peak-picked)")
levels(df_spectra_phases$Spectrum_label) <- ordered_spectrum_labels

# Code which phase group
df_spectra_phases <- df_spectra_phases %>%
  mutate(phase_opposite = phase + 180) %>%
  mutate(phase_opposite = phase_opposite %% 360)

# Establish phase pairs
df_spectra_phases <- df_spectra_phases %>%
  group_by(phase) %>%
  mutate(phase_group = paste(min(c(phase,phase_opposite)),
                             max(c(phase,phase_opposite)),
                             sep = "-")) %>%
  # for each phase group (pair)... 
  group_by(phase_group) %>%
  # code a 1 or 2 for the lower or higher phase
  # within this pair of phases
  mutate(phase_factor = as.numeric(as.factor(phase)))

# Separate phase group
df_phase_groups <- colsplit(df_spectra_phases$phase_group, "-",c("phase_1","phase_2"))

# Add phase info to the big df
df_spectra_phases <- bind_cols(df_spectra_phases, df_phase_groups)
phase_labels <- c("Original","Flipped")
df_spectra_phases$phase_label <- phase_labels[df_spectra_phases$phase_factor]

# Remove working columns
df_spectra_phases[,c("phase_opposite","phase_factor","phase_1","phase_2")] <- NULL

# Cut down the data frame to only contain the columns we need to analyze & plot
df_spectra_phases <- df_spectra_phases %>%
  select(RPO, phase, phase_group, phase_label, 
         Spectrum, Spectrum_label, 
         Frequency, Power_dB, Power_dB_cut)

# Indicate where the electrodes will fall in the cochlea
electrode_positions <- 
  seq(distance_from_apex, distance_from_apex+((num_channels-1)*electrode_spacing), by = electrode_spacing)

#=======================================================#
df_spectra_phases <- df_spectra_phases %>%
  dplyr::filter(phase_group != "180-360")

# Convert from long to wide format,
# so that values can be directly subtracted
spread_columns <- function(d_subset){
  pg <- unique(d_subset$phase_group)
  rpo <- unique(d_subset$RPO)
  message(paste0("Processing ",rpo," RPO ",pg," phase group"))
  
  d_subset %>%
    unique() %>%
    ungroup() %>%
    dplyr::filter(phase_label == "Original") %>%
    cbind(., {d_subset %>%
        unique() %>%
        ungroup() %>%
        dplyr::filter(phase_label == "Flipped") %>%
        select(Flipped = Power_dB)}) %>%
    return()
}

# Run that splitting processing for ALL the parameters
# This step takes a LONG time
# ... but it updates you as it runs. 
df_spectra_phases_flip <- df_spectra_phases %>%
  group_by(RPO, phase_group, Spectrum, Spectrum_label) %>%
  do(., spread_columns(.)) %>%
  bind_rows()

# Calculate dB difference when inverting the phase
df_spectra_phases_flip$flip_depth <- abs(df_spectra_phases_flip$Flipped - df_spectra_phases_flip$Power_dB)

df_spectra_phases_flip$lower <- with(df_spectra_phases_flip, ifelse(Power_dB < Flipped, Power_dB, Flipped))
df_spectra_phases_flip$upper <- with(df_spectra_phases_flip, ifelse(Power_dB > Flipped, Power_dB, Flipped))

# For values that extend below zero,
# set them to be zero,
# to avoid comparisons between a number and (~ -Inf)
df_spectra_phases_flip$Power_dB_0 <- df_spectra_phases_flip$Power_dB
df_spectra_phases_flip$Flipped_0 <- df_spectra_phases_flip$Flipped

# Correct values to zero
df_spectra_phases_flip$Power_dB_0[df_spectra_phases_flip$Power_dB_0 < 0] <- 0
df_spectra_phases_flip$Flipped_0[df_spectra_phases_flip$Flipped_0 < 0] <- 0

# Calculate dB difference when inverting the phase
df_spectra_phases_flip$flip_depth_0 <- abs(df_spectra_phases_flip$Flipped_0 - df_spectra_phases_flip$Power_dB_0)

# Re-track the greater and lesser value of the two phases
# at every frequency point
df_spectra_phases_flip$lower_0 <- with(df_spectra_phases_flip, ifelse(Power_dB_0 < Flipped_0, Power_dB_0, Flipped_0))
df_spectra_phases_flip$upper_0 <- with(df_spectra_phases_flip, ifelse(Power_dB_0 > Flipped_0, Power_dB_0, Flipped_0))

df_spectra_phases_flip$RPO_label <- paste0(df_spectra_phases_flip$RPO, " RPO")

# Make a smaller subset of that data to plot the demonstration figure (Fig 4)
df_spectra_flip_to_plot <- df_spectra_phases_flip %>%
  dplyr::filter(RPO %in% c(0.5, 1.5, 5.5), 
                # Spectrum == "Electric discrete",
                phase_group == "0-180")%>%
  dplyr::filter(Frequency > 125, 
                Frequency < 10000)

# Separate data frame to plot lines for the dense spectra
# (in a smaller line size)
df_spectra_flip_to_plot_dense <- 
  subset(df_spectra_flip_to_plot, (RPO > 2 & Spectrum == "Acoustic ideal"))

# Plot colors
maroon <- "#7a0019"
gold <- "#ffcc33"
light_blue <- "#65B1F3"

#----------------------------------------------------#
px_Figure_4_w_out_lines <- 
  ggplot(df_spectra_flip_to_plot)+
  aes(x = Frequency)+
  # space between phased ripples
  geom_segment(aes(y = lower_0, yend = upper_0,
                   x = Frequency, xend = Frequency,
                   color = flip_depth_0+0.1))+
  scale_color_gradient(high = "black", low = "white",
                       name = "Maximum\nchange\nin dB")+
  # ripples, but excluding high-density acoustic ripples
  # one ripple
  geom_line(
    data = subset(df_spectra_flip_to_plot, !(RPO > 2 & Spectrum == "Acoustic ideal")),
    aes(y = Power_dB_0), color = light_blue, size = 1.1)+
  geom_line(
    data = subset(df_spectra_flip_to_plot, !(RPO > 2 & Spectrum == "Acoustic ideal")),    
    aes(y = Power_dB_0), color = "white", size = 0.4)+
  # other ripple (flipped)
  geom_line(
    data = subset(df_spectra_flip_to_plot, !(RPO > 2 & Spectrum == "Acoustic ideal")),    
    aes(y = Flipped_0), color = maroon, size = 1.1)+
  geom_line(
    data = subset(df_spectra_flip_to_plot, !(RPO > 2 & Spectrum == "Acoustic ideal")),    
    aes(y = Flipped_0), color = "white", size = 0.4)+
  xlab("Frequency (Hz)")+
  ylab("Relative dB")+
  scale_x_log10(breaks = octave_freq_breaks)+
  coord_cartesian(xlim = c(110, 9500),
                  ylim = c(0, 35))+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  theme(strip.text.y = element_text(angle = 0))+
  theme(axis.text.x = element_text(angle = 90))+
  facet_grid(Spectrum_label ~ RPO_label) +
  # lines showing the dense ripple (anything above 2)
  # one ripple
  geom_line(
    data = df_spectra_flip_to_plot_dense,
    aes(y = Power_dB_0), color = light_blue, size = 0.5)+
  geom_line(
    data = df_spectra_flip_to_plot_dense,    
    aes(y = Power_dB_0), color = "white", size = 0.18)+
  # other ripple (flipped)
  geom_line(
    data = df_spectra_flip_to_plot_dense,    
    aes(y = Flipped_0), color = maroon, size = 0.5)+
  geom_line(
    data = df_spectra_flip_to_plot_dense,    
    aes(y = Flipped_0), color = "white", size = 0.18)
# px_Figure_4_w_out_lines
# NOTE: the dense lines don't show up in the R Studio preview pane,
# but they still export normally when you save the plot. 

# add line segment to highlight analysis areas
analysis_low <- 600
analysis_high <- 5000

analysis_full_low <- 125
analysis_full_high <- 10000

df_spectra_phases_flip$Spectrum_label <- droplevels(df_spectra_phases_flip$Spectrum_label)

index_of_last_spectrum_type <- length(levels(df_spectra_phases_flip$Spectrum_label))

df_analysis_range <- 
  data.frame(Spectrum_label = 
               levels(df_spectra_phases_flip$Spectrum_label)[index_of_last_spectrum_type],
             RPO = 5.5,
             RPO_label = "5.5 RPO",
             analysis_low,
             analysis_high,
             analysis_full_low,
             analysis_full_high,
             fullrange_x = 350,
             constrained_x = 1650)

df_analysis_range_left <- 
  data.frame(Spectrum_label = levels(df_spectra_phases_flip$Spectrum_label)[index_of_last_spectrum_type],
             RPO = 0.5,
             RPO_label = "0.5 RPO",
             analysis_low,
             analysis_high,
             analysis_full_low,
             analysis_full_high,
             fullrange_x = 350,
             constrained_x = 2000)

px_Figure_4 <- px_Figure_4_w_out_lines +
  # range on the right lower panel
  geom_segment(data = df_analysis_range, inherit.aes = FALSE,
               aes(x = analysis_low, xend = analysis_high),
               y = 26, yend = 26, color = "black", size = 2)+
  geom_segment(data = df_analysis_range, inherit.aes = FALSE,
               aes(x = analysis_full_low, xend = analysis_full_high),
               y = 32, yend = 32, color = "black", size = 2)+
  geom_label(data = df_analysis_range, inherit.aes = FALSE,
             aes(x = constrained_x), y = 26,
             size = 2.8, label.padding = unit(0.17, "line"),
             label = "cons.")+
  geom_label(data = df_analysis_range, inherit.aes = FALSE,
             aes(x = fullrange_x), y = 32,
             size = 2.8, label.padding = unit(0.17, "line"),
             label = "full")+
  ggtitle(processor)
px_Figure_4

height_fig_4 <- ifelse(processor == "Cochlear", 4, 3.3)
# Save the figure
ggsave(px_Figure_4, file = paste0("Figure_4_flipped_rippled_",processor,".pdf"), height = height_fig_4, width = 8, device = cairo_pdf)
ggsave(px_Figure_4, file = paste0("Figure_4_flipped_rippled_", processor,".png"), height = height_fig_4, width = 8, dpi = 600)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#                     FIGURE 5                          #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Identify maximum dB difference resulting from phase flip
# for overall freq range and also a constrained freq range
df_spectra_flip_max <- df_spectra_phases_flip %>%
  dplyr::filter(!is.na(Frequency)) %>%
  group_by(Spectrum, Spectrum_label, RPO, phase) %>%
  mutate(max_depth_0 = max(flip_depth_0[
    (Frequency > analysis_full_low & 
       Frequency < analysis_full_high)], na.rm = TRUE),
    max_depth_constr_0 = max(flip_depth_0[
      (Frequency > analysis_low & 
         Frequency < analysis_high)], na.rm = TRUE)) %>%
  mutate(is_max_freq = flip_depth_0 == max_depth_0,
         is_max_freq_constr = flip_depth_0 == max_depth_constr_0)

# For the Electric discrete in df_spectra_flip2,
# replace frequency with its nearest center frequency
find_nearest_df <- function(freq){
  num_channels <- length(channel_boundaries.ci)
  center_freqs <- floor(rowMeans(data.frame(channel_boundaries.ci[1:num_channels], channel_boundaries.ci[2:(num_channels+1)])))
  
  nearest_df <- center_freqs[which.min(abs(center_freqs - freq))]
  return(nearest_df)
}

df_spectra_flip_max$Frequency <- 
  ifelse(df_spectra_flip_max$Spectrum=="Electric discrete",
         find_nearest_df(df_spectra_flip_max$Frequency),
         df_spectra_flip_max$Frequency)

max_freqs <- df_spectra_flip_max %>%
  dplyr::filter(is_max_freq == TRUE) %>%
  select(Spectrum, Spectrum_label, RPO, phase, Frequency, flip_depth_0) %>%
  mutate(flip_depth_0 = round(flip_depth_0*100, 0)/100) %>%
  unique() %>%
  mutate(constraint = "Full\nfrequency analysis range\n(125 - 10000 Hz)")

max_freqs_constr <- df_spectra_flip_max %>%
  dplyr::filter(is_max_freq_constr == TRUE) %>%
  select(Spectrum, Spectrum_label, RPO, phase, Frequency, flip_depth_0) %>%
  unique() %>%
  mutate(constraint = "Constrained\nfrequency analysis range\n(600-5000 Hz)")

df_max_freqs <- bind_rows(max_freqs, max_freqs_constr)

# Make a new version of the spectrum labels,
# just for this plot
df_max_freqs$Spectrum_label2 <- df_max_freqs$Spectrum_label
levels(df_max_freqs$Spectrum_label2) <-
  c("Acoustic\n(ideal)",                         
    "Electric\n(topnotopic, discrete)",
    "Summed\nelectric\nactivation\n(tonotopic)",
    "Summed\nelectric\nactivation",
    "Summed\nelectric\nactivation\n(peak-picked)")

px_Figure_5 <- df_max_freqs %>%
  dplyr::filter(!is.na(Frequency)) %>%
  dplyr::filter(Spectrum != "Acoustic ideal") %>%
  ggplot(.)+
  aes(x = RPO, y = flip_depth_0, color = log10(Frequency))+
  # all the little points
  geom_point(
    data = {df_max_freqs %>%
        dplyr::filter(!is.na(Frequency)) %>%
        dplyr::filter(Spectrum != "Acoustic ideal")},
    aes(x = RPO + 0.15),
    size = 1.6, shape = 16,
    position = position_jitter(width = 0.16))+
  # MEAN points
  geom_point(shape = 21,
             fill = "white", color = "black",
             size = 1.6, stroke = 1,
             stat = "summary", fun.y = "mean")+
  coord_cartesian(ylim = c(0, 33))+
  theme_bw()+
  scale_color_gradient2(low = "blue", high = "red", mid = "gray50", 
                        name = "Frequency\n(Hz)\nat max\ndifference",
                        midpoint = log10(2000), 
                        breaks = log10(c(500, 1000, 2000, 4000, 8000)),
                        labels = c(500, 1000, 2000, 4000, 8000))+
  scale_x_continuous(breaks = 1:7,
                     name = "Ripples per octave")+
  ylab("Maximum\ndifference\nin spectral\npower (dB)\nupon phase\ninversion")+
  theme(strip.text.y = element_text(angle = 0),
        axis.title.y = element_text(angle = 0, vjust = 0.5))+
  guides(color = guide_colorbar(
    draw.ulim = TRUE, draw.llim = TRUE, barheight=unit(13, "line")))+
  facet_grid(Spectrum_label2 ~ constraint)+
  ggtitle(processor)
px_Figure_5

ggsave(px_Figure_5, file = paste0("Figure_5_depth_across_phases_",processor,".pdf"), height = 4, width = 7, device = cairo_pdf)
ggsave(px_Figure_5, file = paste0("Figure_5_depth_across_phases_",processor,".png"), height = 4, width = 7, dpi = 600)
#
# End.