# Packages used by this script:
library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr)
library(grid)
library(gridExtra)

# Remove current objects
rm(list = ls())

# Set the path for where you want the plots to be saved:
plots_dir <- "L:\\Manuscripts\\Spectral_Ripple\\Figure_scripts\\CI_ripple_scripts"
setwd(plots_dir)

#=====================================#
# Un-commeent and run only one of the following lines
processor <- "Cochlear"
# processor <- "AB"
# processor <- "MedEl"
#=====================================#
# set the frequency-channel allocation
if (processor == "AB"){
  channel_boundaries.ci <- c(250, 416, 494, 587, 697, 828, 983, 1168, 1387, 
                             1648, 1958,2326, 2763, 3281, 3898, 4630, 8700)
  manufacturer_name <- "Advanced Bionics"
}
if (processor == "Cochlear"){
  channel_boundaries.ci <- c(188,313,438, 563, 688, 813,938,
                             1063, 1188, 1313, 1563, 1813,
                             2063, 2313, 2688, 3063, 3563,
                             4063, 4688, 5313, 6063, 6938, 7938)
  manufacturer_name <- "Cochlear"
}
if (processor == "MedEl"){
  # We only have the center frequencies of the device
  # not frequency bands, liek for the other devices. 
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
  
}

# Turn into a data frame for easy handling. 
channel_boundaries_df <- data.frame(Frequency = channel_boundaries.ci)

# Make data frame of channel frequency lower and upper limits
chan <- data.frame(lower = c(channel_boundaries.ci[1:(length(channel_boundaries.ci)-1)]),
                   upper = c(channel_boundaries.ci[2:(length(channel_boundaries.ci))])) %>%
  mutate(channel = 1:n())

# Code channel bandwidths
chan$oct_diff <- with(chan, (log(upper / lower))/ log(2))
chan$linear_bw <- with(chan, upper - lower)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# FUNCTIONS

make_ripple_spec <- function(RPO){
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
  phase = 0
  
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
  ripple=(1-(mod_depth/2))+
    (sin((((log(Frequencies,2))-log2_freq_offset)*2*pi*ripples_per_octave)+
           ((phase/360)*(2*pi))))*
    (0.5*mod_depth)
  
  # Make a data frame of the rippled amplitude spectrum
  df_ripple <- data.frame(Frequency = Frequencies,
                          Power = ripple,
                          num_frequency_samples = num_frequency_samples,
                          log_freq_samplerate = log_freq_samplerate,
                          RPO = RPO)
  
  # Verify with plot (If you're debugging within this function)
  px_ripple_check <- ggplot(df_ripple)+
    aes(x = Frequency, y = Power)+
    geom_line()+
    scale_x_log10(breaks = c(125, 250, 500, 1000, 2000, 4000, 8000))
  # commented line that woudl otherwise plot it out upon execution
  # px_ripple_check
  
  # Output the data frame 
  # that contains all the info needed for plotting & analysis
  return(df_ripple)
}


send_through_CI_processor <- function(ripple_spec){
  # function that takes a spectrum 
  # (data frame with Frequency & Power)
  # and sends it through a simulation of a 
  # Cochlear Implant speech processor
  # and returns the same spectrum (data frame)
  # with original Power renamed as as `Power_ideal`
  # and a new column with simulated output as `Power`
  
  ripple_spec %>% 
    # assign channel number
    #(0 for freqs below chan 1)
    mutate(channel = findInterval(Frequency, channel_boundaries.ci)) %>%
    # Rename the original power column to be Power_ideal
    mutate(Power_ideal = Power) %>%
    # for each channel... 
    group_by(channel) %>%
    # average all the power that falls within that channel
    mutate(Power = mean(Power)) %>%
    # exit by-channel grouping
    ungroup %>%
    # instead of having two side-by-side columns 
    # of ideal acoustic power and simulated electric power,
    # stack them so that it's long-frmat data,
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
    # return the whole data frame. 
    return()
}

# Combine both functions in one
make_CI_spectrum <- function(RPO) {
  RPO %>% 
    make_ripple_spec() %>% 
    send_through_CI_processor() %>%
    return()
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# These are the RPO densities that we'll plot in the tall figure
RPOs_to_plot <- c(0.3, 0.5, 0.7, 0.9,
                  1.1, 1.3, 1.5, 1.7, 
                  1.9,2.1, 2.3,2.5, 3, 
                  3.5, 4, 4.5, 5, 5.5, 6)

# For each of those RPOs, make the CI spectrum
# (which comes accompanied by the ideal acoustic spectrum)
df_all_ripples <- lapply(RPOs_to_plot, make_CI_spectrum) %>%
  bind_rows()

#-------------------------------#
octave_freqs <- c(125, 250, 500, 1000, 2000, 4000, 8000)

electric_color <- "#20B4E6"
acoustic_color <- "#525252"

# Make nicely formatted labels for the plot
df_all_ripples$RPO_label <- paste0(df_all_ripples$RPO, " ripples\nper octave")
df_all_ripples$RPO_label_short <- paste0(df_all_ripples$RPO, "\nRPO")
df_all_ripples$RPO_label_short_oneline <- paste0(df_all_ripples$RPO, " RPO")

# Plot a small subset of those ripples 
# for demonstration of Acoustic & Electric spectra
px_Figure_1_electric_spectra <- ggplot(subset(df_all_ripples, RPO %in% c(0.5, 1.5, 4)))+
  aes(x = Frequency, y = Power, color=Stimulus)+
  # Acoustic spectrum
  geom_ribbon(data = subset(df_all_ripples, (Stimulus=="Acoustic") &
                              RPO %in% c(0.5, 1.5, 4)),
              aes(ymin = -0.5, ymax = Power), color = "gray40",
              fill = "#DEDEDE", show.legend = FALSE)+
  # Vertical lines separating the channels
  geom_vline(data = channel_boundaries_df, 
             aes(xintercept = Frequency),
             color = "gray50")+
  geom_line(data = {df_all_ripples %>%
      dplyr::filter(RPO %in% c(0.5, 1.5, 4),
                    Stimulus=="Acoustic")},
      size = 1)+
  # just electric
  geom_line(data = {df_all_ripples %>%
      dplyr::filter(RPO %in% c(0.5, 1.5, 4),
                    Stimulus=="Electric",
                    Frequency >= min(channel_boundaries.ci) &
                      Frequency <= max(channel_boundaries.ci))},
      size = 1.8)+
  # FANCY WHITE INLAY
  geom_line(data = {df_all_ripples %>%
      dplyr::filter(RPO %in% c(0.5, 1.5, 4),
                    Stimulus=="Electric",
                    Frequency >= min(channel_boundaries.ci) &
                      Frequency <= max(channel_boundaries.ci))},
      color = "white", size = 0.7)+
  # log axis
  scale_x_log10(breaks = octave_freqs,
                name = "Frequency (Hz)")+
  coord_cartesian(xlim = c(175, 8005), ylim = c(0, 1))+
  scale_color_manual(values = c(acoustic_color, electric_color),
                     labels = c("Acoustic          ",
                                "Electric"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.direction = "horizontal",
        legend.position = "bottom",
        strip.text.y = element_text(angle = 0),
        legend.title = element_blank())+
  facet_grid(RPO_label ~ .)+
  ggtitle(paste0(manufacturer_name," processor channels"))
px_Figure_1_electric_spectra

# # Save it
# ggsave(px_Figure_1_electric_spectra, 
#        file = paste0("Figure_1_",processor,"_electric_spectra_dr_strip.jpg"),
#        height = 3, width = 6.63, dpi = 600)
# 
# ggsave(px_Figure_1_electric_spectra, 
#        file = paste0("Figure_1_",processor,"_electric_spectra_dr_strip.eps"),
#        height = 3, width = 6.63, device = cairo_ps)
# Note: On the manuscript plot, we inserted an intermediate white line 
# in the legend for Electric to match the style in the plot. 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

px_Fig_2_electric_spectra_tall <- ggplot(df_all_ripples)+
  aes(x = Frequency, y = Power, color=Stimulus, size = Stimulus)+
  # Acoustic spectrum
  geom_ribbon(data = subset(df_all_ripples, Stimulus=="Acoustic"),
              aes(ymin = -0.5, ymax = Power), color = "gray40",
              fill = "#DEDEDE", show.legend = FALSE)+
  geom_line(data = subset(df_all_ripples, Stimulus=="Acoustic"),
            size = 0.6)+
  # Vertical lines separating the channels
  geom_vline(data = channel_boundaries_df, aes(xintercept = Frequency),
             color = "gray50")+
  # only draw electric stim within the stimulated frequency range
  geom_line(data = subset(df_all_ripples, 
                          (Stimulus=="Electric" & 
                             Frequency >= min(channel_boundaries.ci) &
                             Frequency <= max(channel_boundaries.ci))),
            size = 1.7)+
  # fancy white inlay
  geom_line(data = subset(df_all_ripples, 
                          (Stimulus=="Electric" & 
                             Frequency >= min(channel_boundaries.ci) &
                             Frequency <= max(channel_boundaries.ci))),
            color = "white", size = 0.7)+
  scale_x_log10(breaks = octave_freqs,
                name = "Frequency (Hz)")+
  scale_y_continuous(name="Normalized Spectral Power",
                     breaks = c(0,0.5,1))+
  scale_size_manual(values = c(0.5, 0.9))+
  coord_cartesian(xlim = c(175, 8005), ylim = c(0, 1))+
  scale_color_manual(values = c(acoustic_color, electric_color),
                     labels = c("Acoustic          ",
                                "Electric"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        # legend.direction = "horizontal",
        legend.position = "none",
        legend.title = element_blank())+
  theme(strip.text.y = element_text(size = 11, angle = 0),
        axis.text.x = element_text(size = 15),
        axis.title = element_text(size = 16))+
  facet_grid(RPO_label_short_oneline ~ .)+
  ggtitle(paste0(manufacturer_name," processor channels"))
px_Fig_2_electric_spectra_tall

# # Save it
# ggsave(px_Fig_2_electric_spectra_tall, 
#        file = paste0("Figure_2_",processor,"_electric_spectra_tall.jpg"), 
#        height = 8, width = 7, dpi = 600)
# 
# ggsave(px_Fig_2_electric_spectra_tall,
#        file = paste0("Figure_2_",processor,"_electric_spectra_tall.eps"),
#        height = 8, width = 7, device = cairo_ps)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Analyze log frequency components in a spectrum
analyze_ripple <- function(df_ripple){
  # input df_ripple should have two columns: Power and Frequency
  # also, this function needs all the other variables that were delcared up there
  
  # Extract the power vector (the amplitude spectrum)
  ripple <- df_ripple$Power
  
  # fft on the log-frequency-sampled ripple
  # make FFT, take absolute value (power spectrum)
  # then take first element up to the element beyond the halfway point
  # (common practice to get power spectrum, 
  # as it's symmetrical)
  ripple_fft <- fft(ripple) %>% abs %>% `[`(., 1:(length(.)/2+1))

  #*******************************************************************# ??????
  log_freq_Nyquist <- df_ripple$log_freq_samplerate[1] / 2
  num_frequency_samples <- df_ripple$num_frequency_samples[1]
  num_frequency_samples <- nrow(df_ripple)
  
  # Establish log-frequency smapling range
  log_freq_start = log(min(df_ripple$Frequency),2)
  log_freq_end = log(max(df_ripple$Frequency),2)
  log_freq_sample_density <- (log_freq_end - log_freq_start) / num_frequency_samples
  log_freq_samplerate = 1/log_freq_sample_density
  
  # Frequency vector from 0 to the Nyquist
  freq_bins <- seq(0,log_freq_Nyquist,by=log_freq_samplerate/num_frequency_samples) 
  #-----------------------------------------------------------------------------#
  # Convert to frequency space rather than bin index
  df_ripple_fft <- data.frame(Log_Mod_Frequency = freq_bins,
                              Power = ripple_fft) %>%
    mutate(bin = 1:n()) 
  
  # Obtain the frequency resolution (separation between frequency samples)
  # by taking the first difference between consecutive values
  df_ripple_log_freq_resolution <- diff(df_ripple_fft$Log_Mod_Frequency)[1]
  
  # Take only the first half, 
  # and exclude first point (which is 0 Hz DC offset)
  df_ripple_fft <- df_ripple_fft[c(2:(nrow(df_ripple_fft)/2)),]
  
  #-----------------------------------------------------------------------------#
  # Find peak (this is the estimated peak ripples per octave density that is represented in the spectrum)
  log_fft_peak <- df_ripple_fft[which.max(df_ripple_fft$Power), "Log_Mod_Frequency"]
  
  #-----------------------------------------------------------------------------#
  # Plot the estimation of ripples per octave present in the original spectrum
  verify_LFMSpectrum <- ggplot(df_ripple_fft)+
    aes(x = Log_Mod_Frequency, y = Power)+
    geom_point()+
    geom_vline(xintercept = log_fft_peak)+
    geom_line()+
    # limit to 7 RPO
    # (the real data go wayyyy higher than that,
    # but don't contain any energy)
    coord_cartesian(xlim = c(0, 7))+
    scale_x_continuous(name = "Log spectral modulation frequency \n(Ripples per octave)",
                       breaks = 1:7)
  # That plots verifies that an idealized rippled spectrum returns a single peak
  # exactly at the target spectral modulation frequency. 
  
  # add some indexical information to the data frame,
  # so that if you bind it with other DFs,
  # it maintains a unique identity
  df_ripple_fft$RPO <- df_ripple$RPO[1]
  df_ripple_fft$peak <- log_fft_peak
  
  return(df_ripple_fft)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# single pipeline to produce RPO analyses based on processing type
analyze_ripple_ci <- function(RPO=NULL){
  make_ripple_spec(RPO) %>%
    # Assign channels to frequency bins
    mutate(channel = findInterval(Frequency, channel_boundaries.ci)) %>%
    # Create column original power 
    mutate(Power_ideal = Power) %>%
    # For each channel... 
    group_by(channel) %>%
    # Summarize power 
    mutate(Power = mean(Power)) %>%
    ungroup() %>% 
    # Analyze the spectrum for log-frequency spectral modulations
    analyze_ripple() %>%
    # Create an additional column for normalized power
    # scaled from 0 to 1
    mutate(Power_norm = Power / max(Power)) %>%
    # Mark peaks
    mark_peaks_in_LFMSpectrum() %>%
    return()
}

analyze_ripple_ideal <- function(RPO) {
  # for analyzing ideal acoustic ripples
  # with no CI processing
  # (just like the previous function,
  # but without the CI processing)
  make_ripple_spec(RPO) %>% 
    mutate(channel = findInterval(Frequency, channel_boundaries.ci)) %>%
    mutate(Power_ideal = Power) %>%
    # group_by(channel) %>%
    # mutate(Power = mean(Power)) %>%
    # ungroup %>%
    analyze_ripple() %>%
    mutate(Processing = "Ideal") %>%
    mark_peaks_in_LFMSpectrum() %>%
    return()
}

mark_peaks_in_LFMSpectrum <- function(LFMSpectrum){
  # Mark absolute peak
  LFMSpectrum$is_peak <- ifelse(LFMSpectrum$Power == max(LFMSpectrum$Power), 1,0)
  
  # Mark frequencies that are candidate sub-peaks 
  # (at least half the power of the peak)
  #
  # 1) first check if points are at least half the max power
  LFMSpectrum$is_subpeak_candidate <- ifelse(LFMSpectrum$Power >= 0.5*max(LFMSpectrum$Power), 1,0)
  #
  # 2) take change in signed derivative of peaks, 
  # mark when they coincide with subpeak candidacy == 1
  LFMSpectrum$Power_deriv_sign <- c(NA, ifelse(diff(LFMSpectrum$Power) < 0, -1, 1))
  #
  # 3) measure if the direction of change switched sign
  LFMSpectrum$Power_deriv_sign_diff <- c(diff(LFMSpectrum$Power_deriv_sign), NA)
  
  # Which indices mark the changed direction of the signed derivative?
  # subtract 1 because the diff score lags by one
  subpeak_indices <- which(with(LFMSpectrum, Power_deriv_sign_diff == -2 & is_subpeak_candidate == 1))
  
  # Establish column variable for subpeaks
  LFMSpectrum$is_subpeak <- 0
  LFMSpectrum$is_subpeak[subpeak_indices] <- 1
  
  # For things marked as the absolute peak,
  # don't also mark them as sub-peaks. 
  LFMSpectrum$is_subpeak[LFMSpectrum$is_peak==1] <- 0
  
  # verify with plot (always!)
  px_spectral_mod_spectrum_peak_picks <- ggplot(LFMSpectrum[1:100,])+
    aes(x = Log_Mod_Frequency, y = Power_norm)+
    geom_point(aes(color = as.factor(Power_deriv_sign)))+
    geom_point(data = subset(LFMSpectrum, is_subpeak==1),
               color = "red", size = 4)+
    geom_line()
  
  # Clean up intermediate columns
  LFMSpectrum[,c("is_subpeak_candidate","Power_deriv_sign","Power_deriv_sign_diff")] <- list(NULL)
  
  return(LFMSpectrum)
}

#**************************************************#
# EXAMPLE ANALYSIS
RPO_vec <- round(seq(0.1,7,0.1),1)

# analyze ripple log-freq modulation spectra
# for each RPO listed above. 
df_rip_analyzed <- lapply(RPO_vec, analyze_ripple_ci) %>% bind_rows()

# Add pretty labels for plotting
df_rip_analyzed$RPO_label <- paste0(df_rip_analyzed$RPO,"\nRPO")
df_rip_analyzed$RPO_label_oneline <- paste0(df_rip_analyzed$RPO," RPO")


# make DF of labels for plotting,
# specifically for adding target markers and labels
RPO_labels <- data.frame(RPO_location = RPO_vec,
                         RPO = RPO_vec) %>%
  mutate(RPO_label = paste0(RPO,"\nRPO"),
         RPO_label_oneline = paste0(RPO," RPO"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# subset of RPOs to plot for the next figure
RPO_to_plot_fig_3 <- c(0.5, 1.1, 1.7, 3, 4, 5)

df_sub <- subset(df_rip_analyzed, RPO %in% RPO_to_plot_fig_3)

RPO_labels_sub <- subset(RPO_labels, RPO %in% RPO_to_plot_fig_3)


px_Figure_3a_mod_power_raw <- ggplot(df_sub)+
  aes(x = Log_Mod_Frequency, y = Power)+
  geom_vline(data = RPO_labels_sub, aes(xintercept = RPO_location),
             color = "gray80", linetype = "solid", size = 1.5)+
  geom_line(size = 1)+
  theme_bw()+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank())+
  ylab("Raw spectral modulation power")+
  xlab("Log modulation frequency\n(Ripple per octave)")+
  coord_cartesian(xlim = c(-0.03, 7.4), expand = FALSE)+
  scale_x_continuous(breaks = 1:7)+
  # scale_y_continuous(limits = c(0, upper_limit_power),
  #                    breaks = c(0, 500, 1000))+
  theme(strip.text.y = element_text(size = 11, angle = 0),
        axis.text.x = element_text(size = 15),
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        panel.spacing = unit(0.7, "line"))+
  facet_grid(RPO_label ~ .)+
  ggtitle(paste0(manufacturer_name," processor channels"))
px_Figure_3a_mod_power_raw


# DF of peaks with proper aesthetic mappings
peak_df_sub <- df_sub %>%
  dplyr::filter(is_peak == 1| is_subpeak == 1) %>%
  melt(id.vars = c("RPO","RPO_label","RPO_label_oneline", 
                   "Log_Mod_Frequency","Power","Power_norm", 
                   "bin","peak")) %>%
  dplyr::filter(value == 1)

px_Figure_3b_mod_power_norm <- px_Figure_3a_mod_power_raw +
  aes(y = Power_norm)+
  # re-name and reset y axis
  ylab("Normalized spectral modulation power")+
  scale_y_continuous(breaks = c(0, 0.5, 1))+
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray60")+
  # Limits of plot axes
  coord_cartesian(xlim = c(-0.03, 7.4), ylim = c(-0.03, 1.22), expand = FALSE)+
  # add points for subpeaks
  geom_point(data=peak_df_sub,
             aes(fill = variable, color = variable),
             size = 1.7, shape = 21,
             alpha = 0.9, stroke = 1.5)+
  scale_fill_manual(values = c(`is_peak` = "black",
                               `is_subpeak` = "white"),
                    labels = c(`is_peak` = "Max peak",
                               `is_subpeak` = "local peak"),
                    name = "Peaks")+
  scale_color_manual(values = c(`is_peak` = "black",
                                `is_subpeak` = "darkred"),
                     labels = c(`is_peak` = "Max peak",
                                `is_subpeak` = "local peak"),
                     name = "Peaks")+
  theme(legend.position = c(0.7, 0.88))+
  theme(legend.title = element_blank())+
  theme(legend.background = element_rect(color = "black"))+
  ggtitle("")
px_Figure_3b_mod_power_norm  


# make a blank plot to insert between them, to add some space
px_spacer <- ggplot(data.frame(x = 1, y = 1))+
  geom_blank()+
  theme_void()

# Combine those plots together
px_Figure_3_mod_power <- arrangeGrob(px_Figure_3a_mod_power_raw+ggtitle(""), 
                                     px_spacer,
                                     px_Figure_3b_mod_power_norm,
                                     nrow = 1,
                                     widths = c(5, 0.5, 5))

# Save it
# ggsave(px_Figure_3_mod_power, 
#        file = paste0("Figure_3_",processor,"_mod_power.jpg"), 
#        height = 4.5, width = 7, dpi = 600)
# ggsave(px_Figure_3_mod_power,
#        file = paste0("Figure_3_",processor,"_mod_power.eps"),
#        height = 4.5, width = 7, device = cairo_ps)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Identify the log-frequency spectral modulation peaks
# For all the input ripple densities sequenced from 0.1 to 7
peak_df_all <- df_rip_analyzed %>%
  dplyr::filter(is_peak == 1| is_subpeak == 1) %>%
  # Convert from wide to long data frame,
  # maintaining indexical info for all these named variables
  melt(id.vars = c("RPO","RPO_label","RPO_label_oneline", 
                   "Log_Mod_Frequency","Power","Power_norm", 
                   "bin","peak"),
       variable.name = "peak_type") %>%
  # keep the real peaks
  dplyr::filter(value == 1) 


# Make data frame of labels for plotting,
# specifically for adding target markers and labels
RPO_labels <- data.frame(RPO_location = RPO_vec,
                         RPO = RPO_vec) %>%
  mutate(RPO_label = paste0(RPO,"\nRPO"),
         RPO_label_oneline = paste0(RPO," RPO"))


px_Figure_4_mod_peaks_bifurcation <- ggplot(peak_df_all)+
  aes(x = RPO, y = Log_Mod_Frequency, color = peak_type, fill = peak_type, size = Power)+
  # ensure that the plot is square,
  # so that the diagonal is always 45 degrees
  coord_fixed(ratio = 1, xlim = c(-0.04, 7.4), ylim = c(-0.04, 7.4))+
  xlab("Ripple-per-octave density of input spectrum")+
  ylab("Ripple-per-octave power\nof output spectrum")+
  theme_bw()+
  geom_point(shape = 21,
             alpha = 1, stroke = 1.2)+
  scale_fill_manual(values = c(`is_peak` = "black",
                               `is_subpeak` = "white"),
                    labels = c(`is_peak` = "Max peak",
                               `is_subpeak` = "local peak"),
                    name = "Peaks")+
  scale_color_manual(values = c(`is_peak` = "black",
                                `is_subpeak` = "darkred"),
                     labels = c(`is_peak` = "Max peak",
                                `is_subpeak` = "local peak"),
                     name = "Peaks")+
  theme(legend.title = element_blank())+
  theme(legend.position = c(0.38, 0.71),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14))+
  guides(size = guide_legend(override.aes = list(shape = 1)),
         color = guide_legend(override.aes = list(size = 3.5)))+
  scale_size_continuous(range = c(0.5, 6),
                        breaks = c(max(peak_df_all$Power)*0.99,
                                   max(peak_df_all$Power)/2,
                                   max(peak_df_all$Power)/4),
                        labels = c("Max possible modulation power",
                                   "3/4 possible modulation power",
                                   "1/2 possible modulation power"))+
  guides(size = guide_legend(override.aes = list(shape = 1)),
         color = guide_legend(override.aes = list(size = 3.5)))+
  ggtitle(paste0(manufacturer_name," processor channels"))
px_Figure_4_mod_peaks_bifurcation


# Save it
# ggsave(px_Figure_4_mod_peaks_bifurcation, 
#        file = paste0("Figure_4_",processor,"_mod_peaks_bifurcation.jpg"),
#        height = 4.2, width = 4.7, dpi = 600)
# ggsave(px_Figure_4_mod_peaks_bifurcation, 
#        file = paste0("Figure_4_",processor,"_mod_peaks_bifurcation.eps"),
#        height = 4.2, width = 4.7, device = cairo_ps)


make_ripple_phase_spec <- function(RPO, phase){
  # like function make_ripple_spec
  # but with an input argument for phase
  # The earlier function is probably redundant,
  # But these different functions made it easy to send
  # into lapply or Map functionals. 
  #
  # create log-spaced ripple from scratch
  # returns a data frame
  # with the specified number of RPO (the only input argument)
  # along with frequency in Hz
  #-----------------------------------------#
  mod_depth = 1
  log2_freq_offset = log(10,2) - 1
  ripples_per_octave = RPO
  # frequency range (log spaced on the front end)
  log_freq_start = log(2^5, 2)
  log_freq_end = log(2^15, 2)  
  
  num_frequency_samples <- 4096*2
  
  # The vector of frequencies
  Frequencies <- 2^(seq(log_freq_start, log_freq_end, length.out = num_frequency_samples))
  Log_Frequencies <- log(Frequencies, 2)
  
  #=============================================================#
  # there are 6 octaves here (128-256-512-1024-2048-4096, 8192)
  # and there are 4096 samples
  # so there are 4096 / 6 samples per octave
  log_freq_sample_density <- (log_freq_end - log_freq_start) / num_frequency_samples
  log_freq_samplerate = 1/log_freq_sample_density
  # this is another way of calculating the log freq samplerate;
  # it gives the same number as the line above
  log_freq_samplerate2 <- 1/diff(Log_Frequencies)[1]
  # 'log_freq_samplerate' should be the number of frequencies sampled for every octave
  
  #--- MAKE THE RIPPLE AMPLITUDE SPECTRUM ---#
  ripple=(1-(mod_depth/2))+
    (sin((((log(Frequencies,2))-log2_freq_offset)*2*pi*ripples_per_octave)+
           ((phase/360)*(2*pi))))*
    (0.5*mod_depth)
  
  # Produce the output data frame
  df_ripple <- data.frame(Frequency = Frequencies,
                          Power = ripple,
                          num_frequency_samples = num_frequency_samples,
                          log_freq_samplerate = log_freq_samplerate,
                          RPO = RPO,
                          phase = phase)
  
  # verify with plot
  px_ripple_check <- ggplot(df_ripple)+
    aes(x = Frequency, y = Power)+
    geom_line()+
    scale_x_log10(breaks = c(125, 250, 500, 1000, 2000, 4000, 8000))
  # px_ripple_check
  
  return(df_ripple)
}

# combine both functions in one
make_CI_spectrum_w_phase <- function(RPO, phase) {
  make_ripple_phase_spec(RPO, phase) %>% 
    send_through_CI_processor() %>%
    return()
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# make a large data frame with inverted-phase ripples
# for each ripple density,
# measure the maximum difference in power within each channel
# across all phases

# Parameters
RPOs_to_cross <- seq(0.5, 8, 0.5)
phases_to_cross <- seq(0, 350, 10)

# All the combinations of parameters
params_to_cross <- expand.grid(RPO = RPOs_to_cross,
                               phase = phases_to_cross)

# Send those parameters into the function,
# and make a big data frame out of it. 
df_all_phases <- Map(f = make_CI_spectrum_w_phase,
                     RPO = params_to_cross$RPO,
                     phase = params_to_cross$phase) %>%
  bind_rows() %>%
  dplyr::filter(Frequency <= 10000,
                Frequency >= 100)

# identities of phases that could be inverted 180 degrees
# phase_groups <- paste(0:179, 180:359, sep = "-")
phase_groups_A <- 0:179
phase_groups_B <- 180:359
phase_groups <- paste(phase_groups_A, phase_groups_B, sep = "-")

df_all_phases <- df_all_phases %>%
  mutate(phase_opposite = phase + 180) %>%
  # circle shift numbers higher than 360
  mutate(phase_opposite = phase_opposite %% 360)


# Establish phase group
df_all_phases <- df_all_phases %>%
  group_by(phase) %>%
  # indicate phase group (e.g. 0-180, 10-190, etc.)
  # and convert 180-0 and 190-10 into 0-180 and 10-190
  # by taking the minimum first, then the maximum phase
  # within each group
  # (this will avoid duplicates for later)
  mutate(phase_group = paste(min(c(phase,phase_opposite)),
                             max(c(phase,phase_opposite)),
                             sep = "-")) %>%
  group_by(phase_group) %>%
  mutate(phase_factor = as.numeric(as.factor(phase)))

# Separate phase group
df_phase_groups <- colsplit(df_all_phases$phase_group, "-",c("phase_1","phase_2"))

# Add phase info to the ig df
df_all_phases <- bind_cols(df_all_phases, df_phase_groups)
phase_labels <- c("Original","Flipped")
df_all_phases$phase_label <- phase_labels[df_all_phases$phase_factor]

# add plot-friendly labels 
df_all_phases$RPO_label <- paste0(df_all_phases$RPO, " ripples\nper octave")
df_all_phases$RPO_label_short <- paste0(df_all_phases$RPO, "\nRPO")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# begin with a test case of only phase flipping
RPOs_to_plot <- c(0.5, 1.5, 4)
phases_to_plot <- c(0, 180)

octave_freqs <- c(125, 250, 500, 1000, 2000, 4000, 8000)

electric_color <- "#20B4E6"
acoustic_color <- "#525252"

df_inverted_phase <- df_all_phases %>%
  dplyr::filter(RPO %in% RPOs_to_plot,
                phase %in% phases_to_plot,
                Frequency >= min(channel_boundaries.ci),
                Frequency <= max(channel_boundaries.ci),
                Stimulus=="Electric")

# Put the power values for 0 and 180 side-by-side
# by converting from long to wide data frame
# (gather)
df_inverted_phase_wide <- df_inverted_phase %>%
  ungroup() %>%
  select(RPO, RPO_label, phase, Power, Frequency, channel) %>%
  spread(data = ., key = phase, Power)

# For each frequency bin and each channel,
# Mark minimum power across both phases
df_inverted_phase_wide$lower <- with(df_inverted_phase_wide, ifelse(`0` < `180`, `0`, `180`))
# ... and maximum power
df_inverted_phase_wide$upper <- with(df_inverted_phase_wide, ifelse(`0` > `180`, `0`, `180`))

# For each frequency bin and each phase,
# marke the difference in power between flipped phases
df_inverted_phase_wide <- df_inverted_phase_wide %>%
  # find the difference between whichever phase is higher/lower
  mutate(phase_flip_depth = upper - lower) %>%
  group_by(RPO) %>%
  # For each RPO, mark a column that tells whether that frequency
  # is the one with maximum change in power resulting from phase inversion
  mutate(is_max_depth = phase_flip_depth == max(phase_flip_depth)) %>%
  dplyr::filter(Frequency > 150, Frequency < 8000)

# colors
phase_flip_color <- "#F5B700"
phase_flip_color <- "#D8C18B"

px_Fig_5_demo_phase_inversion <- ggplot(df_inverted_phase)+
  aes(x = Frequency, y = Power, color = as.factor(phase))+
  geom_vline(data = channel_boundaries_df, aes(xintercept = Frequency),
             color = "#E5E5E5")+
  # space between phase flips
  geom_ribbon(data = df_inverted_phase_wide, 
              inherit.aes = FALSE,
              aes(x = Frequency, ymin = lower, ymax = upper,
                  fill = is_max_depth))+
  scale_fill_manual(values= c(`TRUE` = "#797979", `FALSE` = "gray80"))+
  # electric spectrum
  geom_line(size = 1.8)+
  # FANCY WHITE INLAY
  geom_line(color = "white", size = 0.7, aes(group = phase))+
  scale_x_log10(breaks = octave_freqs,
                name = "Frequency (Hz)")+
  coord_cartesian(xlim = c(175, 8005), ylim = c(0, 1))+
  scale_color_manual(values = c(phase_flip_color, electric_color))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  facet_grid(RPO_label ~ . )+
  ggtitle(paste0(manufacturer_name," processor channels"))
px_Fig_5_demo_phase_inversion


# aesthetic changes so the line-within-a-line looks continuous. 

px_Fig_5_demo_phase_inversion_2 <- ggplot(df_inverted_phase)+
  aes(x = Frequency, y = Power, color = as.factor(phase))+
  geom_vline(data = channel_boundaries_df, aes(xintercept = Frequency),
             color = "#E5E5E5")+
  # space between phase flips
  geom_ribbon(data = df_inverted_phase_wide, 
              inherit.aes = FALSE,
              aes(x = Frequency, ymin = lower, ymax = upper,
                  fill = is_max_depth))+
  scale_fill_manual(values= c(`TRUE` = "#797979", `FALSE` = "gray80"))+
  # electric spectrum
  geom_line(data = subset(df_inverted_phase, phase == 0), size = 1.8)+
  # FANCY WHITE INLAY
  geom_line(data = subset(df_inverted_phase, phase == 0), color = "white", size = 0.7, aes(group = phase))+
  # electric spectrum
  geom_line(data = subset(df_inverted_phase, phase == 180), size = 1.8)+
  # FANCY WHITE INLAY
  geom_line(data = subset(df_inverted_phase, phase == 180), color = "white", size = 0.7, aes(group = phase))+
  scale_x_log10(breaks = octave_freqs,
                name = "Frequency (Hz)")+
  scale_y_continuous(breaks = seq(0,1,0.5))+
  ylab("Normalized spectral power")+
  coord_cartesian(xlim = c(175, 8005), ylim = c(0, 1))+
  scale_color_manual(values = c(phase_flip_color, electric_color))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  theme(strip.text.y = element_text(angle = 0))+
  facet_grid(RPO_label ~ . )+
  ggtitle(paste0(manufacturer_name," processor channels"))
px_Fig_5_demo_phase_inversion_2

# # Save it
# ggsave(px_Fig_5_demo_phase_inversion_2, 
#        file = paste0("Figure_5_",processor,"_demo_phase_inversion.jpg"),
#        height = 3.4, width = 5.5, dpi = 300)
# ggsave(px_Fig_5_demo_phase_inversion_2, 
#        file = paste0("Figure_5_",processor,"_demo_phase_inversion.eps"),
#        height = 3.4, width = 5.5, device = cairo_ps)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Quantify the difference in spectral power
# when phase is inverted

channel_power_sum <- df_all_phases %>%
  dplyr::filter(Stimulus == "Electric") %>%
  dplyr::filter(channel %in% 1:22) %>%
  ungroup() %>%
  group_by(Stimulus, RPO, channel, phase_group) %>%
  summarise(Power_min = min(Power),
            Power_max = max(Power)) %>%
  mutate(Power_diff = Power_max - Power_min) 

# Make separate data frames that specifically indicate
# max (for a given phase configuration), 
# median (across all channel & phases)
# and overall max (scanning across all phases)

channel_power_sum_max <- channel_power_sum %>%
  ungroup() %>%
  group_by(RPO, phase_group) %>%
  summarise(max_depth_channel = channel[Power_diff == max(Power_diff)],
            Power_diff_max = max(Power_diff))

channel_power_sum_median <- channel_power_sum %>%
  ungroup() %>%
  group_by(RPO) %>%
  summarise(Power_diff_median = median(Power_diff))

channel_power_sum_max_overall <- channel_power_sum_max %>%
  ungroup() %>%
  group_by(RPO) %>%
  summarise(phase_group_max = phase_group[Power_diff_max == max(Power_diff_max)],
            max_depth_channel = max_depth_channel[Power_diff_max == max(Power_diff_max)],
            Power_diff_median = median(Power_diff_max),
            Power_diff_max = max(Power_diff_max))


# Plot it out
px_Fig_6_modulation_depths <- ggplot(channel_power_sum) +
  aes(x = RPO, y = Power_diff, group = RPO) +
  # don't show outliers in this layer;
  # they will be shown in the next layer
  geom_boxplot(outlier.colour = NA)+
  geom_point(position = position_jitter(width = 0.12, height = 0),
             alpha = 0.2, shape = 16, 
             aes(color = "single channel, single phase"))+
  # Overall max
  geom_point(data = channel_power_sum_max_overall,
             inherit.aes = FALSE,
             aes(x = RPO, y = Power_diff_max, color = "Across all phases"),
             size = 4, shape = 17)+
  # Median (with white background coded as stroke)
  geom_point(data = channel_power_sum_median,
             inherit.aes = FALSE,
             aes(x = RPO, y = Power_diff_median, fill = "Median depth"),
             size = 3.5, color = "white", shape = 21, stroke = 1.3)+
  scale_x_continuous(breaks = 0:8)+
  xlab("Spectral density (RPO)")+
  ylab("Maximum difference\nin spectral power\nupon phase inversion")+
  coord_cartesian(ylim = c(-0.05, 1.5), xlim = c(-0.05, 8.5), expand = FALSE)+
  scale_y_continuous(breaks = seq(0,1, 0.25))+
  theme_bw() +
  # label the channel with overall max 
  geom_label(data = channel_power_sum_max_overall,
             inherit.aes = FALSE,
             aes(x = RPO, y = Power_diff_max + 0.15, label = max_depth_channel),
             color = "firebrick", size = 3)+
  scale_color_manual(values = c(`Across all phases` = "firebrick",
                                `single channel, single phase` = "black"),
                     name = "Maximum achievable\nmodulation depth change")+
  scale_fill_manual(values = c(`Median depth` = "blue"), name = "")+
  # Customize the legend to make it more helpful
  guides(color = guide_legend(override.aes = list(shape = c(17, 16),
                                                  size = c(4, 1.5)),
                              order = 1),
         fill = guide_legend(name = NULL, title.position = "right"))+
  # Position the legend just so
  theme(legend.position = c(0.75, 0.7),
        legend.spacing.y = unit(0.4, "line"),
        legend.key.height = unit(0.6, "line"))+
  annotate("label",
           label = "Channel with maximum change\nin modulation depth",
           x = 1.8, y = 1.36, color = "firebrick",
           size = 2.9)+
  ggtitle(paste0(manufacturer_name," processor channels"))
px_Fig_6_modulation_depths

# # Save the plot
# ggsave(px_Fig_6_modulation_depths, 
#        file = paste0("Figure_6_",processor,"_modulation_depths.jpg"),
#        height = 3.15, width = 6.1, dpi = 300)
# ggsave(px_Fig_6_modulation_depths, 
#        file = paste0("Figure_6_",processor,"_modulation_depths.eps"),
#        height = 3.15, width = 6.1, device = cairo_ps)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Plot bandwidth of all the channels
all_linear_BWs <- unique(sort(chan$linear_bw))
smallest_bw <- all_linear_BWs[1]
largest_bw <- all_linear_BWs[length(all_linear_BWs)]

px_channel_linear <- ggplot(chan)+
  aes(x = channel, y = linear_bw)+
  geom_point()+
  geom_label(aes(label = channel))+
  ylab("Bandwidth (Hz)")+
  scale_y_continuous()+
  coord_cartesian(ylim = c(smallest_bw*0.75, largest_bw*1.1))+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid.minor = element_blank())+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle(paste0(manufacturer_name," processor channels"))
px_channel_linear


px_channel_octaves <- ggplot(chan)+
  aes(x = channel, y = oct_diff)+
  geom_point()+
  geom_label(aes(label = channel))+
  ylab("Bandwidth (octaves)")+
  coord_cartesian(ylim = c(-0.02, 0.9))+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle("")
px_channel_octaves

px_Figure_7_channel_BWs <- arrangeGrob(px_channel_linear, px_channel_octaves, nrow = 2)

# # Save the plot
# ggsave(px_Figure_7_channel_BWs, 
#        file = paste0("Figure_7_",processor,"_channel_BWs.jpg"),
#        height = 4.8, width = 5.1, dpi = 300)
# 
# ggsave(px_Figure_7_channel_BWs, 
#        file = paste0("Figure_7_",processor,"_channel_BWs.eps"),
#        height = 4.8, width = 5.1, device = cairo_ps)

# Tell me what the minimum octave space that can fit in the narrowest channel
min_oct_space <- min(chan$oct_diff)
1/(2*min_oct_space)

# AB min single-channel octave space is 0.2475606
# AB max single-channel spectral period = 2.019707 octave

# Cochlear min single-channel octave space is 0.1443321
# Cochlear max single-channel spectral period = 3.464233 octave

# Med-El min single-channel octave space is 0.39794
# Med-El max single-channel spectral period = 1.256471 octave
