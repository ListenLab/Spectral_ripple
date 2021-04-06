# create spectra data frames for Figures 1 and 2
# Packages used by this script:
library(ggplot2)
library(dplyr)
library(reshape2)
library(cowplot)
library(purrr)
# Remove current objects from workspace
rm(list = ls())

# Set the path for where you want the plots to be saved:
plots_dir <- "C:\\Type\\Your\\File\\Path\\Here"
setwd(plots_dir)

# Remove unnecessary dplyr verbose messages
options(dplyr.summarise.inform=FALSE)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Choose your CI processor
# Un-comment and run only one of the following lines
# processor <- "Cochlear"
# processor <- "AB_mid_scala"
# processor <- "AB_slim_J"
# processor <- "MedEl"
processor <- "Octave"
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Greenwood functions
freq_to_mm <- function(freq=NULL, A = 165.4, a = 2.1, length = 35, k = 0.88){
  mm <- log10((freq/A)+k)*length/a
  return(mm)
}
mm_to_freq <- function(position=NULL, A = 165.4, a = 2.1, length = 35, k = 0.88){
  freq <- A*((10^(a*position/length))-k)
  return(freq)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# set the frequency-channel allocation and other details
# for the CI processors
if (processor == "AB_mid_scala"){
  channel_boundaries.ci <- c(250, 416, 494, 587, 697, 828, 983, 1168, 1387, 
                             1648, 1958,2326, 2763, 3281, 3898, 4630, 8700)
  processor_name <- "Advanced Bionics Mid-Scala"
  num_channels <- 16
  electrode_spacing = 0.975
  array_length <- 15
  # fill this in if you want
  electrode_locations_ideal <- NA
}

if (processor == "AB_slim_J"){
  channel_boundaries.ci <- c(250, 416, 494, 587, 697, 828, 983, 1168, 1387, 
                             1648, 1958,2326, 2763, 3281, 3898, 4630, 8700)
  processor_name <- "Advanced Bionics Slim J"
  num_channels <- 16
  electrode_spacing = 1.3
  array_length <- 20
  # fill this in if you want
  electrode_locations_ideal <- NA
}

if (processor == "Cochlear"){
  channel_boundaries.ci <- c(188,313,438, 563, 688, 813,938,
                             1063, 1188, 1313, 1563, 1813,
                             2063, 2313, 2688, 3063, 3563,
                             4063, 4688, 5313, 6063, 6938, 7938)
  num_channels <- 22
  electrode_spacing = 0.75
  # Note: there is variable spacing between electrodes 
  # on various cochlear arrays
  # fill this in if you want
  electrode_locations_ideal <- NA
  processor_name <- "Cochlear"
}
if (processor == "MedEl"){
  # We only have the center frequencies of the device
  # not frequency bands, like for the other devices. 
  # this is because the Med-El device uses bell-shaped filters
  # instead of FFT processing. 
  # So by calculating these filters for the Med-El device,
  # we are doing two important things:
  # 1) violating the shape of the filters
  # 2) doing some fancy calculations to derive the intermediate 
  #   frequency boundaries between channels,
  #   if they were to exist (they don't; see previous note)
  
  center_freqs <- c(125, 234, 385, 582, 840, 1182, 1631, 2227, 3064, 4065, 5658, 7352)
  
  places <- freq_to_mm(center_freqs)
  intermediate_place_boundaries <- 
    c(
      # border of first channel
      freq_to_mm(100),
      # place-derived boundaries between center frequencies of channels 2-11
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
  processor_name <- "Med-El"
  num_channels <- 12
  electrode_spacing = 2.4
  # 12 electrodes over 26.4mm (26.4/11 = 2.4)
  # fill this in if you want
  electrode_locations_ideal <- NA
}
if (processor == "Octave"){
  channel_boundaries.ci <- 
    c(222.7249, 280.6152, 353.5522, 445.4468, 561.2265,
      707.0994, 890.8874, 1122.4452, 1414.1891, 1781.7624,
      2244.8748, 2828.3585, 3563.5002, 4489.7185,
      5656.6778, 7126.9510)
  
  center_freqs <- 
    c(250.0000, 314.9795, 396.8484, 499.9965, 629.9547,
      793.6914, 999.9861, 1259.9007, 1587.3717,
      1999.9584, 2519.7839, 3174.7214, 3999.8891,
      5039.5328, 6349.3988)
  processor_name <- "perfect 1/3-octave spacing"
  
  num_channels <- 15
  
  electrode_spacing = 1.3
  electrode_locations_ideal <- freq_to_mm(center_freqs)
  array_length <- 20
  
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# INITIATE RIPPLE PARAMETERS
# RPO <- 1
dB_rolloff_per_mm <- 2
phase <- 225 
# this phase is not special, 
# but it allows integer RPOs to align very nicely
# atop clean octave intervals marked on the graph
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Turn into a data frame for easy handling. 
channel_boundaries_df <- data.frame(Frequency = channel_boundaries.ci)

# Make data frame of channel frequency lower and upper limits
chan <- data.frame(lower = c(channel_boundaries.ci[1:(length(channel_boundaries.ci)-1)]),
                   upper = c(channel_boundaries.ci[2:(length(channel_boundaries.ci))])) %>%
  mutate(channel = 1:n())

chan$center <- floor(rowMeans(chan[,c("lower","upper")]))

# override that if we're explicitly declaring center freqs
if(processor %in% c("MedEl", "Octave")){
  chan$center <- center_freqs
}

# Code channel bandwidths
chan$oct_diff <- with(chan, (log(upper / lower))/ log(2))
chan$linear_bw <- with(chan, upper - lower)

# Shift in insertion depth for the most apical electrode
# in an ideal case, it would be this:
distance_from_apex <- freq_to_mm(chan$center[1])
# but you can override it here, to create a spectral shift
distance_from_apex <- 10

# Indicate where the electrodes will fall in the cochlea
electrode_positions <- 
  seq(distance_from_apex, distance_from_apex+((num_channels-1)*electrode_spacing), by = electrode_spacing)

# override ideal electrode positions if you have
# computed the ideal positions
# and declared them in a vector
if(exists("electrode_locations_ideal")){
  if (!is.na(electrode_locations_ideal)){
  electrode_positions <- electrode_locations_ideal
  }
}

spectral_modulation_depth <- 30
# DC shift of spectral power
# (doesn't have any effect for this analysis,
# since there's no compression simulated here)
intensity_offset <- 0
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#               * * *   MAIN  FUNCTIONS  * * *                   #
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
  # phase = 0
  
  # Frequency range (log spaced on the front end)
  log_freq_start = log(2^5, 2)
  log_freq_end = log(2^15, 2)  
  
  num_frequency_samples <- 4096*2
  
  # The vector of frequencies to sample the spectrum
  Frequencies <- 2^(seq(log_freq_start, log_freq_end, length.out = num_frequency_samples))
  Log_Frequencies <- log(Frequencies, 2)
  
  #=============================================================#
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
  
  # Verify with plot (If you're debugging within this function)
  # px_ripple_check <- ggplot(df_ripple)+
  #   aes(x = Frequency, y = Power)+
  #   geom_line()+
  #   scale_x_log10(breaks = c(125, 250, 500, 1000, 2000, 4000, 8000))
  # px_ripple_check
  
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
  
  ripple_spec %>% 
    # assign channel number
    mutate(channel = findInterval(Frequency, channel_boundaries.ci)) %>%
    # nullify channel "zero"
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
    # `Stimulus` indicates Acoustic / Electric
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

activate_in_cochlea <- function(CI_ripple_spectrum, 
                                maintain_ideal = FALSE,
                                peak_pick = FALSE){
  # For each channel, find the cochlear position
  # and the level of activation
  # given the shift in mm (distance from apex)
  #---------------------------------------------#
  # maintain_ideal: don't simulate actual electrode position;
  # just treat the electrodes as if they were 
  # tonotopically matched
  #---------------------------------------------#
  electrode_positions <- CI_ripple_spectrum %>%
    pull(electrode_position_mm) %>%
    unique() %>%
    na.omit() %>%
    as.numeric()
  
  ideal_target_frequencies <- chan$center
  
  # Position of electrodes on the basilar membrane:
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
      # mutate(electrode_position_mm = electrode_position_mm + shift_mm) %>%
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
      # arrange with loudest channels on top
      arrange(desc(Power_dB)) %>%
      # anything below the top 8 channels is set to -60
      within(., Power_dB[9:nrow(.)] <- -60) %>%
      # re-arrange 
      arrange(RPO, channel)
  }
  return(cochlear_activation_sites)
}

spread_db_per_mm <- function(cf, dB, rolloff_per_mm){
  # Simulate_spread of activation
  # spanning from 32 Hz to 32768 Hz
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
  #---------------------------------------------#
  # Takes a single row from the data frame output 
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

simulate_spread_all_channels <- 
  function(cochlear_activation_spots, dB_rolloff_per_mm, mod_depth){
    # Execute the simulate-spread function
    # for each cochlear site,
    # maintaining data frame metadata about
    # spectral density, phase, channel, etc. 
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
  # 10-log instead of 20log
  # because the electrical signals are out of phase
  reference_level <- 0.00002
  dB <- 10 * log10(pressure_Pa / reference_level)
  return(dB)
}

dB_to_pressure_Pa <- function(dB){
  # 10-log instead of 20log
  # because the electrical signals are out of phase
  pressure_Pa <- 0.00002 * 10^(dB/10)
  return(pressure_Pa)
}

`%dB+%` <- function(lhs, rhs){
  # function for adding decibel values 
  intermediate_sum <- dB_to_pressure_Pa(lhs) + dB_to_pressure_Pa(rhs)
  output <- pressure_Pa_to_dB(intermediate_sum)
  return(output)
}

# # example:
# 60 %dB+% 60

sum_spread_channels <- function(electric_spectral_channels, peak_dB){
  #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
  # add up the power at each position,
  # making sure to add dB values 
  # using the special %dB+% function
  #----------------------------------#
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
  
  # This is where we add the CUT values 
  # and then MELT them into long format 
  apical_edge_Power_dB <- 
    electric_spectral_channels_summed$Power_dB[
      which.min(abs(electric_spectral_channels_summed$Frequency - apical_frequency_boundary))]
  
  basal_edge_Power_dB <- 
    electric_spectral_channels_summed$Power_dB[
      which.min(abs(electric_spectral_channels_summed$Frequency - basal_frequency_boundary))]
  
  edge_dB_min <- min(c(apical_edge_Power_dB, basal_edge_Power_dB))
  
  # copy into new value column
  electric_spectral_channels_summed$Power_dB_cut <- electric_spectral_channels_summed$Power_dB
  
  # walk down both the lower and upper edges of the spectrum
  electric_spectral_channels_summed$Power_dB_cut[
    (electric_spectral_channels_summed$Frequency < apical_frequency_boundary & 
       electric_spectral_channels_summed$Power_dB_cut < edge_dB_min)] <- edge_dB_min
  
  electric_spectral_channels_summed$Power_dB_cut[
    (electric_spectral_channels_summed$Frequency > basal_frequency_boundary & 
       electric_spectral_channels_summed$Power_dB_cut < edge_dB_min)] <- edge_dB_min
  
  # "DC component"
  # fine energy level at the edge of the spectrum
  # to create a "zero-pad" signal that can be analyzed via FFT
  electric_spectral_channels_summed$Power_dB_cut_AC <- 
    electric_spectral_channels_summed$Power_dB_cut - edge_dB_min
  
  # # verify
  # ggplot(electric_spectral_channels_summed)+
  #   aes(x = Frequency, y = Power_dB)+
  #   geom_line()+
  #   geom_line(aes(y = Power_dB_cut+3), color = "blue")+
  #   annotate("label", x = 100, y = edge_dB_min+3, label = "0")+
  #   geom_line(aes(y = Power_dB_cut_AC), color = "red")+
  #   scale_x_log10()
  
  return(electric_spectral_channels_summed)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# String all the functions together
make_channel_activation_pattern <- function(RPO, phase, maintain_ideal, peak_pick = FALSE){
  # Wrap all the relevant functions into one
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
# # example:
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

make_summed_activation_pattern <- 
  function(RPO, phase, maintain_ideal, peak_pick = FALSE){
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# # Create the activation patterns
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
make_spectra_comparison <- function(RPO, phase, full_freq_sampling = FALSE){
  # Function to put all the different types of spectra
  # into a single labeled data frame
  # progress message for the user:
  message(paste0("building spectra for ",RPO, " RPO, ",phase," phase"))
  
  spectrum_acoustic_ideal <- make_CI_spectrum(RPO, phase) %>% 
    dplyr::filter(Stimulus == "Acoustic") %>%
    mutate(Spectrum = "Acoustic ideal") 
  
  spectrum_electric_ideal <- make_CI_spectrum(RPO, phase) %>% 
    dplyr::filter(Stimulus == "Electric") %>%
    mutate(Spectrum = "Electric discrete")
  
  spectrum_electric_summed <- 
    make_summed_activation_pattern(RPO, phase, maintain_ideal = FALSE) %>%
    mutate(Spectrum = "Summed activation")
  
  spectrum_electric_summed_picked <- 
    make_summed_activation_pattern(RPO, phase, maintain_ideal = FALSE, peak_pick = TRUE) %>%
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
  
  # Nullify the power level below the min frequency
  df_spectra[
    (df_spectra$Spectrum == "Electric discrete" & 
       df_spectra$Frequency < min(channel_boundaries.ci)), "Power_dB"] <- NA
  # Nullify the power level above the max frequency
  df_spectra[
    (df_spectra$Spectrum == "Electric discrete" & 
       df_spectra$Frequency > max(channel_boundaries.ci)), "Power_dB"] <- NA
  
  # Order the types of spectra (for plotting)
  df_spectra$Spectrum <- 
    factor(df_spectra$Spectrum,
           levels = c("Acoustic ideal","Electric discrete",
                      "Summed activation ideal",
                      "Summed activation",
                      "Summed activation picked"))
  
  # Add plot-friendly labels
  df_spectra$Spectrum_label <- df_spectra$Spectrum
  ordered_spectrum_labels <- c("Acoustic\n(ideal)",
                               "Electric\n(tonotopic, discrete)",
                               "Summed\nelectric activation\n(tonotopic)",
                               "Summed\nelectric activation",
                               "Summed\nelectric activation\n(peak-picked)")
  levels(df_spectra$Spectrum_label) <- ordered_spectrum_labels
  
  return(df_spectra)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# PLOT 1
# parameters to plot:
RPOs_to_plot <- c(0.5, 1.5, 4)
inverted_phase <- 
  ifelse(phase > 180, phase - 180, phase + 180)
phases_to_plot <- c(rep(phase, length(RPOs_to_plot)),
                    rep(inverted_phase, length(RPOs_to_plot)))

# doubling the RPOs to plot since we're plotting two phases
RPOs_to_plot <- rep(RPOs_to_plot, 2)

# Make the acoustic / electric spectra
df_spectra <- map2(.x = RPOs_to_plot, 
                   .y = phases_to_plot, 
                   .f = make_spectra_comparison) %>%
  bind_rows()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Activation pattern for each channel
channel_activation_patterns_ideal <- 
  map2(.x = RPOs_to_plot,
       .y = phases_to_plot, 
       .f = make_channel_activation_pattern, 
       maintain_ideal = TRUE) %>%
  bind_rows() %>%
  mutate(
    Spectrum = "Summed activation ideal",
    Spectrum_label = "Summed\nelectric activation\n(tonotopic)")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# activation patterns positioned within the cochlea
# (by electrode position, not tonotopic)
channel_activation_patterns <- 
  map2(.x = RPOs_to_plot,
       .y = phases_to_plot, 
       .f = make_channel_activation_pattern, 
       maintain_ideal = FALSE) %>%
  bind_rows() %>%
  mutate(
    Spectrum = "Summed activation",
    Spectrum_label = "Summed\nelectric activation")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Remove peak-picked versions if it isn't a Cochlear processor
# (easier to do it here rather than conditional binds above)
if (processor != "Cochlear"){
  df_spectra <- df_spectra[df_spectra$Spectrum != "Summed activation picked",]
  channel_activation_patterns <- channel_activation_patterns[channel_activation_patterns$Spectrum != "Summed activation picked",]
  channel_activation_patterns_ideal <- 
    channel_activation_patterns_ideal[channel_activation_patterns_ideal$Spectrum != "Summed activation picked",]
}
if (processor == "Cochlear"){
  channel_activation_patterns_picked <- 
    map2(.x = RPOs_to_plot,
         .y = phases_to_plot, 
         .f = make_channel_activation_pattern, 
         maintain_ideal = FALSE, peak_pick = TRUE) %>%
    bind_rows() %>%
    mutate(
      Spectrum = "Summed activation picked",
      Spectrum_label = "Summed\nelectric activation\n(peak-picked)")
} else {channel_activation_patterns_picked <- NULL}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Put all the data together into one data frame
channel_activation_patterns_all <- 
  bind_rows(channel_activation_patterns,
            channel_activation_patterns_picked,
            channel_activation_patterns_ideal) %>%
  mutate(Spectrum = factor(Spectrum,
                           levels = levels(df_spectra$Spectrum)))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Only maintain peak-picked data for Cochlear devices
if (processor != "Cochlear"){
  channel_activation_patterns_all <- channel_activation_patterns_all %>%
    dplyr::filter(Spectrum != "Summed activation picked")
  
  df_spectra <- df_spectra %>%
    dplyr::filter(Spectrum != "Summed activation picked")
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ensure proper order for plotting
ordered_spectrum_types <- 
  c("Acoustic\n(ideal)",
    "Electric\n(tonotopic, discrete)",
    "Summed\nelectric activation\n(tonotopic)",
    "Summed\nelectric activation\n(on CI array)",
    "Summed\nelectric activation\n(on CI array)\n(peak-picked)")

channel_activation_patterns_all$Spectrum_label <- channel_activation_patterns_all$Spectrum
levels(channel_activation_patterns_all$Spectrum_label) <- ordered_spectrum_types

df_spectra$RPO_label <- paste0(df_spectra$RPO, " RPO")
channel_activation_patterns_all$RPO_label <- paste0(channel_activation_patterns_all$RPO, " RPO")
levels(df_spectra$Spectrum_label) <- ordered_spectrum_types

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Reduce size of channel_activation_patterns_all
# by only keeping frequencies
# within the plotting range
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

# Produce a vector of indices that are log spaced
reduced_frequency_indices <- (length(sampled_frequencies) + 1) - 
  (10^seq(0, log10(length(sampled_frequencies)), length.out = 1600)) %>%
  round(., 0) %>%
  unique() %>%
  rev()

sampled_frequencies_reduced <- sampled_frequencies[reduced_frequency_indices]

channel_activation_patterns_all_to_plot <-
  channel_activation_patterns_all[
    channel_activation_patterns_all$Frequency %in% sampled_frequencies_reduced,]

# Shrink the data frame by decimating frequency info
df_spectra_to_plot <-
  df_spectra[
    df_spectra$Frequency %in% sampled_frequencies_reduced,]

# Omit peak-picked spectra for non-Cochlear devices
df_spectra_to_plot$Spectrum_label <- droplevels(df_spectra_to_plot$Spectrum_label)
df_spectra_to_plot$Spectrum <- droplevels(df_spectra_to_plot$Spectrum)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# create wide data frame 
# for plotting inter-phase energy difference
# make wide data frame for the acoustic and elec-ideal
df_spec_sub <- df_spectra_to_plot %>%
  dplyr::filter(Spectrum %in% c("Acoustic ideal","Electric discrete")) %>%
  dplyr::filter(Frequency > min(channel_boundaries.ci),
                Frequency < max(channel_boundaries.ci))

df_spec_sub_w <- df_spec_sub %>% 
  dplyr::filter(Spectrum %in% c("Acoustic ideal","Electric discrete")) %>%
  dplyr::filter(phase == 225) %>% 
  mutate(Power_225 = Power_dB) %>%
  bind_cols(.,
            {df_spec_sub %>% dplyr::filter(phase == 45) %>% 
                mutate(Power_45 = Power_dB) %>%
                select(Power_45)})

# define upper and lower edges of the ribbon between phases
df_spec_sub_w$lower <- ifelse(df_spec_sub_w$Power_225 < df_spec_sub_w$Power_45, df_spec_sub_w$Power_225, df_spec_sub_w$Power_45)
df_spec_sub_w$upper <- ifelse(df_spec_sub_w$Power_225 > df_spec_sub_w$Power_45, df_spec_sub_w$Power_225, df_spec_sub_w$Power_45)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Plot details
octave_freq_breaks <- c(125, 250, 500, 1000, 2000, 4000, 8000)

electric_blue <- "#83C2F8"
spectrum_colors <- c(`Acoustic ideal` = "black",
                     `Electric discrete` = electric_blue,
                     `Summed activation` = "#6B0F1A",
                     `Summed activation picked` = "#28536B",
                     `Summed activation ideal` = "gray50")

spectrum_labels <- data.frame(
  Spectrum_label = factor(levels(df_spectra_to_plot$Spectrum_label)),
  Spectrum = factor(levels(df_spectra_to_plot$Spectrum)),
  # which column to put the label in
  RPO_label = paste0(min(RPOs_to_plot)," RPO"),
  # x-axis position of the label
  Frequency = 125,
  # y-axis position of the label
  Power_dB = 27,
  label = LETTERS[1:length(unique(df_spectra_to_plot$Spectrum))])

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# remove frequencies that extend beyond the viewing window
df_spectra_to_plot[
  (df_spectra_to_plot$Spectrum=="Acoustic ideal" & 
     df_spectra_to_plot$Frequency < min(channel_boundaries.ci)),
  c("Power_dB","Power_dB_cut")] <- NA

df_spectra_to_plot[
  (df_spectra_to_plot$Spectrum=="Acoustic ideal" & 
     df_spectra_to_plot$Frequency > max(channel_boundaries.ci)),
  c("Power_dB","Power_dB_cut")] <- NA

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# PLOT
px_Figure_1_example_spectra_channels <- 
  ggplot(subset(df_spectra_to_plot, phase == 225))+
  aes(x = Frequency, y = Power_dB, color = Spectrum)+
  # lower limit of frequency sampling
  annotate("rect", xmin = 1, xmax = min(channel_boundaries.ci),
           ymin = -spectral_modulation_depth, ymax = spectral_modulation_depth*2,
           fill = "gray80", alpha = 0.5)+
  # upper limit of frequency sampling
  annotate("rect", xmin = max(channel_boundaries.ci), xmax = 30000,
           ymin = -spectral_modulation_depth, ymax = spectral_modulation_depth*2,
           fill = "gray80", alpha = 0.5)+
  # individual channel activation patterns
  geom_ribbon(
    # for only one phase
    #(otherwise it's too cluttered)
    data = subset(channel_activation_patterns_all_to_plot, phase == 225),
    aes(ymin = -Inf, ymax = Power_dB,
        x = Frequency, group = channel),
    fill = "black", alpha = 0.3, color = NA)+
  # Overall summed activation pattern
  geom_line(size = 1.2)+
  # summed activation pattern, 
  # with sidebands cut away
  geom_line(size = 0.9, aes(y = Power_dB_cut), linetype = "longdash")+
  # label on the left column
  geom_label(data = spectrum_labels, aes(label = label))+
  geom_text(data = spectrum_labels, aes(label = label), color = "black")+
  xlab("Frequency (Hz)")+
  ylab("Relative dB")+
  scale_x_log10(breaks = octave_freq_breaks)+
  coord_cartesian(xlim = c(110, 9500),
                  ylim = c(0, 33))+
  scale_color_manual(values = spectrum_colors)+
  theme_bw()+
  theme(legend.position = "none")+
  theme(panel.grid.minor = element_line(size = 0.4))+
  theme(panel.grid.major = element_blank())+
  theme(strip.text.y = element_text(angle = 0))+
  facet_grid(Spectrum_label ~ RPO_label)+
  ggtitle(paste0("Example ripple spectral for ",processor_name," device"))
px_Figure_1_example_spectra_channels

# plot energy difference
# between phase inversions
px_Figure_1_example_spectra_channels_phase_inv <- 
  px_Figure_1_example_spectra_channels + 
  # space between phase inversions
  geom_ribbon(data = df_spec_sub_w,
              aes(ymin = lower, ymax = upper),
              color = NA,
              fill = "khaki", alpha = 0.5)+
  # phase-inverted (thin) line
  geom_line(data = 
              {df_spectra_to_plot %>% 
                  dplyr::filter(phase == 225-180,
                                Frequency > min(channel_boundaries.ci),
                                Frequency < max(channel_boundaries.ci))},
            # linetype = "dashed",
            color = "red",
            alpha = 0.3, size = 0.45)+
  # re-plot the black lines on top
  # # main (thick) line
  geom_line(size = 1.5,
            data = subset(df_spec_sub, (phase == 225 & RPO < 4)))+
  # thick line on panel B 4 rpo acoustic
  geom_line(size = 0.8,
            data = subset(df_spec_sub, 
                          (phase == 225 & RPO == 4 & Spectrum == "Acoustic ideal")))+
  # thick line on panel B 4 rpo
  geom_line(size = 1.5,
            data = subset(df_spec_sub, (phase == 225 & RPO == 4 & Spectrum == "Electric discrete")))+
  # thin white line for the electric spectra
  geom_line(size = 0.25, color = "white",alpha = 0.5,
            data = subset(df_spec_sub, (phase == 225 & Spectrum == "Electric discrete")))
px_Figure_1_example_spectra_channels_phase_inv

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# save the plot
plot_output_height = ifelse(processor == "Cochlear", 4.9, 4.14)
plot_output_width = 9.35

ggsave(px_Figure_1_example_spectra_channels,
       file = paste0("Figure_1_example_spectra_channels_" , processor,".png"),
       height = plot_output_height, width = plot_output_width, dpi = 300)
ggsave(px_Figure_1_example_spectra_channels,
       file = paste0("Figure_1_example_spectra_channels_" , processor,".pdf"),
       height = plot_output_height, width = plot_output_width, device = cairo_pdf)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#                                                       #
#               * * *   FIGURE 2   * * *                #
#                                                       #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# MANY RPOs for a spectral modulation density analysis
RPOs_to_analyze <- seq(0.2, 7.5, 0.2)
# all in the same phase
phases_to_analyze <- rep(225, length(RPOs_to_analyze))

# Generate spectra for each of those RPO densities,
# combine into a single big data frame
df_spectra_across_RPOs <- map2(.x = RPOs_to_analyze, 
                               .y = phases_to_analyze, 
                               .f = make_spectra_comparison,
                               full_freq_sampling = TRUE) %>%
  bind_rows()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Analyze log frequency density in log-spaced ripple
analyze_ripple <- function(df_ripple, HPF = 0.2){
  # df_ripple should have two columns: Power_dB and Frequency
  # also, this function needs all the other variables that were declared up there
  
  # Extract the power vector (the amplitude spectrum)
  ripple <- df_ripple$Power_dB_cut
  ripple[is.na(ripple)] <- min(ripple, na.rm = TRUE)
  # fft on the log-frequency-sampled ripple
  # make FFT, take absolute value (power spectrum)
  # then take first element up to the element beyond the halfway point
  # (common practice to get power spectrum, 
  # as it's symmetrical)
  ripple_fft <- fft(ripple) %>% abs %>% `[`(., 1:(length(.)/2+1))
  
  # Some parameters used to create the ripples:
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
  
  log_freq_Nyquist <- log_freq_samplerate / 2
  num_frequency_samples <- 4096*2
  # num_frequency_samples <- df_ripple$num_frequency_samples[1]
  # I wonder if this part isn't correct? 
  # maybe `log_freq_samplerate` isn't really the end of the vector here. 
  # or maybe the `by=` portion is too large?
  num_frequency_samples <- nrow(df_ripple)
  
  # Establish log-frequency smapling range
  log_freq_start = log(min(df_ripple$Frequency),2)
  log_freq_end = log(max(df_ripple$Frequency),2)
  log_freq_sample_density <- (log_freq_end - log_freq_start) / num_frequency_samples
  log_freq_samplerate = 1/log_freq_sample_density
  
  freq_bins <- seq(0,log_freq_Nyquist,by=log_freq_samplerate/num_frequency_samples); # frequency vector from 0 to the Nyquist
  # summary(freq_bins)
  #-----------------------------------------------------------------------------#
  # convert to frequency space rather than bin index
  df_ripple_fft <- data.frame(Log_Mod_Frequency = freq_bins,
                              Power = ripple_fft) %>%
    mutate(bin = 1:n()) 
  
  # Obtain the frequency resolution (separation between frequency samples)
  # by taking the first difference between consecutive values
  df_ripple_log_freq_resolution <- diff(df_ripple_fft$Log_Mod_Frequency)[1]
  
  # Take only the first half of the FFT
  # and exclude 1st point (Log_Mod_Frequency of 0, which is DC offset)
  df_ripple_fft <- df_ripple_fft[c(2:(nrow(df_ripple_fft)/2)),]
  
  # Exclude electric offset
  df_ripple_fft$Power[df_ripple_fft$Log_Mod_Frequency < HPF] <- NA
  #-----------------------------------------------------------------------------#
  # Find peak (this is the estimated peak ripples per octave 
  #   density that is represented in the spectrum)
  log_fft_peak <- df_ripple_fft[which.max(df_ripple_fft$Power), "Log_Mod_Frequency"]
  
  #-----------------------------------------------------------------------------#
  # Plot the estimation of ripples per octave present in the original spectrum
  verify_LFMSpectrum <- ggplot(df_ripple_fft)+
    aes(x = Log_Mod_Frequency, y = Power)+
    geom_point()+
    geom_vline(xintercept = log_fft_peak)+
    geom_line()+
    coord_cartesian(xlim = c(0, 7))+
    scale_x_continuous(name = "Log spectral modulation frequency \n(Ripples per octave)",
                       breaks = 0:7)
  # # un-comment to see the verification:
  # verify_LFMSpectrum
  
  # Add some indexical information to the data frame,
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
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# set steady-state amplitude in the middle of the depth
# at edge frequencies
min_freq <- channel_boundaries.ci[1]
max_freq <- channel_boundaries.ci[length(channel_boundaries.ci)]

# Lower edge
df_spectra_across_RPOs$Power_dB_cut[
  (df_spectra_across_RPOs$Spectrum == "Acoustic ideal" &
     df_spectra_across_RPOs$Frequency < min_freq)] <- spectral_modulation_depth/2

# Upper edge
df_spectra_across_RPOs$Power_dB_cut[
  (df_spectra_across_RPOs$Spectrum == "Acoustic ideal" &
     df_spectra_across_RPOs$Frequency > max_freq)] <- spectral_modulation_depth/2

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Run the analysis on all the RPOs
df_spectral_density <- df_spectra_across_RPOs %>%
  group_by(RPO, phase, Spectrum, Spectrum_label) %>%
  do(., analyze_ripple(., HPF = 0.2)) %>%
  dplyr::filter(Log_Mod_Frequency <= 8)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Normalize modulation power within each RPO
df_spectral_density <- df_spectral_density %>%
  group_by(RPO, phase, Spectrum, Spectrum_label) %>%
  mutate(Power_dB_norm = Power_dB / max(Power_dB, na.rm = TRUE))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Produce the sinc function that represents the distortion 
max_ideal_modulation_depth <- 
  df_spectral_density %>%
  dplyr::filter(Spectrum == "Acoustic ideal") %>%
  group_by(RPO) %>%
  summarise(max_depth = max(Power_dB, na.rm = TRUE)) %>%
  mutate(correction_factor = spectral_modulation_depth / max_depth)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Plot correction factor 
# for maximum possible modulation depth
ggplot(max_ideal_modulation_depth)+
  aes(x = RPO, y = max_depth)+
  geom_point()+
  geom_line()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Encode that into the main data frame
df_spectral_density <- left_join(df_spectral_density, max_ideal_modulation_depth)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Correct possible modulation depth
df_spectral_density$Power_dB_corrected <- 
  df_spectral_density$Power_dB * df_spectral_density$correction_factor

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Labels for the simulations (for the plot)
ordered_spectrum_types <- 
  c("Acoustic\n(ideal)",
    "Electric\n(tonotopic, discrete)",
    "Summed\nelectric activation\n(tonotopic)",
    "Summed\nelectric activation\n(on CI array)",
    "Summed\nelectric activation\n(on CI array)\n(peak-picked)")

levels(df_spectral_density$Spectrum_label) <- ordered_spectrum_types

# Only maintain peak-picked for Cochlear devices
if (processor != "Cochlear"){
  df_spectral_density <- df_spectral_density %>%
    dplyr::filter(Spectrum != "Summed activation picked")
}

# Calculate maxima
df_spectral_density %>%
  dplyr::filter(Log_Mod_Frequency > 0.24) %>%
  group_by(Spectrum) %>%
  summarise(max_modulation = max(Power_dB_corrected))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# PLOT
# Exclude super-low densities (0.24)
# because their power will wreck the scale
px_RPO_IO_map <- 
  df_spectral_density %>%
  dplyr::filter(Log_Mod_Frequency > 0.24) %>%
  ggplot(.)+
  aes(x = RPO, y = Log_Mod_Frequency, fill = Power_dB_corrected)+
  geom_tile()+
  scale_fill_gradientn(colours = c("white", "black", "blue"),
                       guide = "legend", 
                       # plot top of fill map as blue,
                       # fill 0 - 16 as white to black
                       breaks = c(0,3, 6, 18, 30),
                       name = "Modulation\npower (dB)")+
  scale_x_continuous(breaks = 1:7)+
  scale_y_continuous(breaks = 1:7)+
  xlab("Input spectral modulation density (ripples per octave)")+
  ylab(paste0(processor_name,"\n\n",
              "Output\nspectral modulation density\n(ripples per octave)"))+  theme_bw()+
  theme(legend.justification=c(0, 1.6))+
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))+
  theme(legend.key.height = unit(1, "line"))+
  # Ensure perfect diagonal
  coord_fixed(ratio = 1, expand = FALSE)+
  facet_grid(. ~ Spectrum_label)
# px_RPO_IO_map

# Normalized modulation depth:
# Not part of the paper,
# but useful to speculate about mapping
# a spectral modulation across 
# the entire dynamic range.
px_RPO_IO_map_norm <-
  df_spectral_density %>%
  dplyr::filter(Log_Mod_Frequency > 0.24) %>%
  ggplot(.)+
  aes(x = RPO, y = Log_Mod_Frequency, fill = Power_dB_norm)+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "black",
                      name = "Proportion of\nmodulation\npower within\ninput RPO", 
                      breaks = seq(0,1,0.25))+
  scale_x_continuous(breaks = 1:7)+
  scale_y_continuous(breaks = 1:7)+
  xlab("Input spectral modulation density (ripples per octave)")+
  ylab("Normalized Output\nspectral modulation density\n(ripples per octave)")+
  theme_bw()+
  theme(legend.justification=c(0, 4.5))+
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))+
  theme(legend.title = element_text(vjust = 1))+
  theme(legend.key.height = unit(1, "line"))+
  coord_fixed(ratio = 1, expand = FALSE)+
  facet_grid(. ~ Spectrum_label)
# px_RPO_IO_map_norm

# # Assemble the plots together in ONE image
# # with *unified* color-fill legend
# px_Figure_2_ripple_density_IO <- 
#   cowplot::plot_grid(px_RPO_IO_map,
#                        # ggtitle(paste0("Ripple input-output map for ",processor_name," device")),
#                        # NULL,
#                      px_RPO_IO_map_norm,
#                      align = "v",
#                      axis = "b",
#                      labels = "AUTO", 
#                      label_x = 0.05, hjust = 1, label_size=17,
#                      ncol = 1, nrow = 2)
# px_Figure_2_ripple_density_IO


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Set up the plots with customized (scoped)
# color-fill maps

px_RPO_IO_map_ideal <- 
  df_spectral_density %>%
  dplyr::filter(Log_Mod_Frequency > 0.5) %>%
  dplyr::filter((Spectrum %in% c("Acoustic ideal","Electric discrete"))) %>%
  ggplot(.)+
  aes(x = RPO, y = Log_Mod_Frequency, fill = (Power_dB))+
  geom_tile()+
  scale_fill_gradientn(colours = c("white", "black", "blue"),
                       guide = "legend", 
                       name = "Modulation power (dB)")+
  scale_x_continuous(breaks = 1:7)+
  scale_y_continuous(breaks = 1:7)+
  xlab("Input spectral modulation density (ripples per octave)")+
  ylab("Output\nspectral modulation density\n(ripples per octave)")+
  theme_bw()+
  theme(legend.position = "bottom")+
  guides(fill = 
           guide_legend(title.position = "top", 
                        title.hjust = 0.5,
                        label.position = "bottom"))+ 
  coord_fixed(ratio = 1, expand = FALSE, 
              xlim = c(0.5, 7.5), ylim = c(0.5, 7.5))+
  facet_grid(. ~ Spectrum_label)
# px_RPO_IO_map_ideal

px_RPO_IO_map_no_ideal <- 
  df_spectral_density %>%
  dplyr::filter(Log_Mod_Frequency > 0.5) %>%
  dplyr::filter(!(Spectrum %in% c("Acoustic ideal","Electric discrete"))) %>%
  ggplot(.)+
  aes(x = RPO, y = Log_Mod_Frequency, fill = (Power_dB))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "black",
                      guide = "legend", 
                      name = "Modulation power (dB)")+
  scale_x_continuous(breaks = 1:7)+
  scale_y_continuous(breaks = 1:7)+
  xlab("Input spectral modulation density (ripples per octave)")+
  ylab("Output\nspectral modulation density\n(ripples per octave)")+
  theme_bw()+
  theme(legend.position = "bottom")+
  guides(fill = 
           guide_legend(title.position = "top", 
                        title.hjust = 0.5,
                        label.position = "bottom"))+
  coord_fixed(ratio = 1, expand = FALSE, 
              xlim = c(0.5, 7.5), ylim = c(0.5, 7.5))+
  facet_grid(. ~Spectrum_label)
# px_RPO_IO_map_no_ideal



# arrange in sub-plots to combine
px_RPO_IO_map_no_ideal_sub_plot <- 
  px_RPO_IO_map_no_ideal +
  theme(axis.title.y = element_text(color = "white"))
# px_RPO_IO_map_no_ideal_sub_plot


px_RPO_IO_map_ideal_sub_plot <- 
  px_RPO_IO_map_ideal +
  scale_fill_gradient(low = "white", high = "black",
                      guide = "legend", 
                      name = "Modulation power (dB)")+
  facet_grid(. ~ paste0("\n",Spectrum_label, 
                        ifelse(processor == "Cochlear","\n","")))
# px_RPO_IO_map_ideal_sub_plot

# dynamically set the widths 
# because plot-arrangement packages
# can't align plots *and* keep aspect ratio :(
rel_width2 <- 
  ifelse(processor == "Cochlear", 2.65,2)

px_Figure_2_ripple_density_IO <- 
  cowplot::plot_grid(
    px_RPO_IO_map_ideal_sub_plot+
      ggtitle(paste0("Ripple input-output for ",processor_name," device")),
    px_RPO_IO_map_no_ideal_sub_plot+
      ggtitle(""), 
    rel_widths = c(2, rel_width2), 
    nrow = 1, align = "v", axis = "tb"
  )
px_Figure_2_ripple_density_IO

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Save the plot
p2_width <- ifelse(processor == "Cochlear", 11.75, 9.5)
# p2_width <- ifelse(processor == "Cochlear", 13, 10.5)

ggsave(px_Figure_2_ripple_density_IO,
       file = paste0("Figure_2_ripple_density_IO_maps_" , processor,".png"),
       height = 5.5, width = p2_width, dpi = 300)

ggsave(px_Figure_2_ripple_density_IO,
       file = paste0("Figure_2_ripple_density_IO_maps_" , processor,".pdf"),
       height = 5.5, width = p2_width, device = cairo_pdf)

# END

