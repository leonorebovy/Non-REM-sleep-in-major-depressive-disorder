## Script to analyse all datasets together
## Project: Non-REM sleep in major depressive disorder
## Author: Leonore Bovy

## Load packages ----
library(tidyverse)
library(plyr)
library(cowplot)

## Source functions ----
source("funcs/theme.R")
source("funcs/rainclouds.R")
source("dataset_a.R")
source("dataset_b.R")
source("dataset_c.R")

## Add NAs for dataset C; does not have finger tapping data
dataset_c$CONS_change <- rep(NA, nrow(dataset_c) )

## Get colnumnnames to match on
colnames_a <- colnames(dataset_a)
colnames_b <- colnames(dataset_b)
colnames_c <- colnames(dataset_c)

common_colnames <- intersect(colnames_a,  colnames_b)
common_colnames <- intersect(common_colnames,  colnames_c)

## Filter only the columns that are common
dataset_a <- dataset_a %>% dplyr::select(common_colnames)
dataset_b <- dataset_b %>% dplyr::select(common_colnames)
dataset_c <- dataset_c %>% dplyr::select(common_colnames)

dataset_a$dataset <- "A"
dataset_b$dataset <- "B"
dataset_c$dataset <- "C"

dataset_a_comb <- dataset_a
dataset_b_comb <- dataset_b
dataset_c_comb <- dataset_c

## Rename to combine based on length of medication
dataset_a_comb$group<- revalue(dataset_a_comb$group, c("Patients" = "Medicated_long"))
dataset_b_comb$group<- revalue(dataset_b_comb$group, c("control" = "Controls",
                                             "patient_unmed" = "Unmedicated",
                                             "patient_7d" = "Medicated_short"))
dataset_c_comb$group<- revalue(dataset_c_comb$group, c("control" = "Controls",
                                             "patient_28" = "Medicated_long",
                                             "patient_7" = "Medicated_short"))

dataset_comb <- rbind(dataset_a_comb, dataset_b_comb, dataset_c_comb)


## Rename to distinguish between datasets
dataset_a$group<- revalue(dataset_a$group, c("Controls" = "Controls (A)",
                                             "Patients" = "Medicated (A)"))
dataset_b$group<- revalue(dataset_b$group, c("control" = "Controls (B)",
                                             "patient_unmed" = "Unmedicated (B)",
                                             "patient_7d" = "Medicated 7d (B)"))
dataset_c$group<- revalue(dataset_c$group, c("control" = "Controls (C)",
                                             "patient_28" = "Medicated 28d (C)",
                                             "patient_7" = "Medicated 7d (C)"))

dataset_all <- rbind(dataset_a, dataset_b, dataset_c)


## Plots ----

## spindles
plot_rain_spindens_all <- ggplot(dataset_all,aes(x=group,y=spin_dens,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, seed = 18),size=2.5, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=spin_dens),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Spindle density [/epoch]')+ xlab('')+ ylim(0,3.5) +
  raincloud_theme +theme(axis.text.x = element_text(size = 15)) + 
  guides(fill=FALSE,colour=FALSE)+ 
  scale_color_manual(values=color_dataset_all) + 
  scale_fill_manual(values=color_dataset_all) + 
  ggtitle("Group difference in spindle density") +
  geom_signif(annotation=formatC("**", digits=1), y_position=3.3,xmin=1.2, xmax=2.2, 
              textsize = 5,colour="BLACK") +
  geom_signif(annotation=formatC("**", digits=1), y_position=3.3,xmin=6.2, xmax=8.2, 
              textsize = 5,colour="BLACK")

plot_rain_spinamp_all <- ggplot(dataset_all,aes(x=group,y=spin_amp,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, seed = 18),size=2.5, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=spin_amp),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Spindle amplitude [µV]')+ xlab('')+ ylim(0,50) +
  raincloud_theme + theme(axis.text.x = element_text(size = 15)) +
  guides(fill=FALSE,colour=FALSE)+ 
  scale_color_manual(values=color_dataset_all) + 
  scale_fill_manual(values=color_dataset_all) + 
  ggtitle("Group difference in spindle amplitude") +
  geom_signif(annotation=formatC("**", digits=1), y_position=45,xmin=3.2, xmax=4.2, 
              textsize = 5,colour="BLACK")


## SOs
plot_rain_soamp_all <- ggplot(dataset_all,aes(x=group,y=so_amp,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, seed = 18),size=2.5, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=so_amp),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('SO amplitude [µV]')+ xlab('')+ ylim(75,300) +
  raincloud_theme +theme(axis.text.x = element_text(size = 15)) + 
  guides(fill=FALSE,colour=FALSE)+ 
  scale_color_manual(values=color_dataset_all) + 
  scale_fill_manual(values=color_dataset_all) + 
  ggtitle("Group difference in SO amplitude") +
  geom_signif(annotation=formatC("**", digits=1), y_position=270,xmin=1.2, xmax=2.2, 
              textsize = 5,colour="BLACK")+
  geom_signif(annotation=formatC("**", digits=1), y_position=270,xmin=6.2, xmax=7.2, 
              textsize = 5,colour="BLACK")+
  geom_signif(annotation=formatC("**", digits=1), y_position=290,xmin=6.2, xmax=8.2, 
              textsize = 5,colour="BLACK")

plot_rain_sodur_all <- ggplot(dataset_all,aes(x=group,y=so_dur,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, seed = 18),size=2.5, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=so_dur),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('SO duration [s]')+ xlab('')+ ylim(1.05,1.6) +
  raincloud_theme + theme(axis.text.x = element_text(size = 15)) +
  guides(fill=FALSE,colour=FALSE)+ 
  scale_color_manual(values=color_dataset_all) + 
  scale_fill_manual(values=color_dataset_all) + 
  ggtitle("Group difference in SO duration") +
  geom_signif(annotation=formatC("**", digits=1), y_position=1.52,xmin=1.2, xmax=2.2, 
              textsize = 5,colour="BLACK")  +
  geom_signif(annotation=formatC("*", digits=1), y_position=1.52, xmin=4.2, xmax=5.2, 
              textsize = 5,colour="BLACK") 


## Coupling
plot_rain_amount_coup_all <- ggplot(dataset_all,aes(x=group,y=coupling_occur,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=2.5, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=coupling_occur),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Amount SO-fast spindles')+ xlab('')+ ylim(0,300) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE) +
  scale_color_manual(values=color_dataset_all) + 
  scale_fill_manual(values=color_dataset_all) + 
  ggtitle("Group difference in amount of SO-fast spindle couplings") 

plot_rain_mean_delay_all  <- ggplot(dataset_all,aes(x=group,y=mean_delay,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=2.5, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=mean_delay),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Mean delay [s]')+ xlab('')+ ylim(0.35,0.75) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE) +
  scale_color_manual(values=color_dataset_all) + 
  scale_fill_manual(values=color_dataset_all) + 
  ggtitle("Group difference in mean delay") 

plot_rain_sd_delay_all <- ggplot(dataset_all,aes(x=group,y=sd_delay,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=2.5, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=sd_delay),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Delay dispersion (SD)')+ xlab('')+ ylim(0.15,0.4) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE) +
  scale_color_manual(values=color_dataset_all) + 
  scale_fill_manual(values=color_dataset_all) + 
  ggtitle("Group difference in delay dispersion (SD)") +
  geom_signif(annotation=formatC("*", digits=1), y_position=0.375,xmin=1.2, xmax=2.2,
              textsize = 6, colour="BLACK")+
  geom_signif(annotation=formatC("*", digits=1), y_position=0.375,xmin=3.2, xmax=5.2,
              textsize = 6, colour="BLACK")+
  geom_signif(annotation=formatC("*", digits=1), y_position=0.37,xmin=6.2, xmax=7.2,
              textsize = 6, colour="BLACK") +
  geom_signif(annotation=formatC("*", digits=1), y_position=0.395,xmin=6.2, xmax=8.2,
              textsize = 6, colour="BLACK")


spin_plot_dataset_all   <- ggarrange(plot_rain_spindens_all, plot_rain_spinamp_all,nrow = 2, labels = c("A", "B") )
so_plot_dataset_all     <- ggarrange(plot_rain_soamp_all,plot_rain_sodur_all,nrow = 2, labels = c("A", "B") )
coupling_plot_dataset_all     <- ggarrange(plot_rain_amount_coup_all, plot_rain_mean_delay_all,plot_rain_sd_delay_all,nrow = 3, labels = c("A", "B", "C") )

# ggsave("figs/spin_plot_dataset_all_v5.png", spin_plot_dataset_all, device = "png", dpi = 300, width = 15, height = 10)
# ggsave("figs/so_plot_dataset_all_v5.png", so_plot_dataset_all, device = "png", dpi = 300, width = 15, height = 10)
# ggsave("figs/coupling_plot_dataset_all_v5.png", coupling_plot_dataset_all, device = "png", dpi = 300, width = 15, height = 11)


## Freq plots 
## Dataset A
freq_spec_dataset_a <- dataset_a_long_pvalues_perm %>%
  ggplot(aes(x = freq, y = mean_densDB, colour = group)) + 
  theme_classic()+ my_theme_freq + ylim(0,25)+
  stat_summary(fun.y="mean", geom = 'smooth', se = FALSE)+              
  scale_color_manual(values=color_dataset_a, labels = c("Controls (A)", "Medicated (A)")) +  
  ylab("Adjusted power [dB]") +
  theme(plot.title  = element_text(hjust = 0.5)) +                                
  theme(      
    axis.text.y  = element_text(size = 22, color="#000000"),
    axis.title.y =  element_text(size = 22, color="#000000"),
    axis.title.x = element_blank(),
    axis.line.x  = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.length  = unit(0.25,"cm"),
    axis.ticks.y  = element_line(size = 1.5),
    legend.position = 'none') + 
  NULL

perm_plot_dataset_a <- ggplot(dataset_a_long_pvalues_perm, aes(x=freq, y=abs(log10(perm_pvalue)), colour="#000000")) +
  theme_classic()+ my_theme_freq +
  stat_summary(fun.y="mean", geom="bar", colour="#939393", fill = "#939393")+ 
  xlab("Frequency [Hz]") + ylab("P-values: Medicated vs. controls") +
  theme(       
    axis.text.y = element_text(size = 18, color="#000000"),
    axis.text.x = element_text(size = 18, color="#000000"), 
    axis.title.y = element_text(size = 22, color="#000000"),
    axis.title.x = element_text(size = 18, color="#000000"), 
    axis.ticks  = element_line(size = 1.5),
    axis.ticks.length  = unit(0.25,"cm"),
    legend.position='none')+ 
  geom_hline(yintercept = abs(log10(0.05)))+
  scale_y_continuous(limits = c(0, 4),labels=c("0",  "10e-1","10e-2", "10e-3","10e-4"))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) 
NULL


## Dataset B
freq_spec_dataset_b_three <- dataset_b_long_pvalues_perm_three %>%
  ggplot(aes(x = freq, y = mean_densDB, colour = group)) + 
  theme_classic()+ my_theme_freq + ylim(0,25)+
  stat_summary(fun.y="mean", geom = 'smooth', se = FALSE)+              
  scale_color_manual(values=color_dataset_b_all, labels = c("Controls - 7days (B)", "Medicated - 7days (B)", "Medicated - 28days (B)")) +  
  ylab(" ") +
  theme(plot.title  = element_text(hjust = 0.5)) +                                
  theme(      
    axis.text.y  = element_text(size = 22, color="#000000"),
    axis.title.x = element_blank(),
    axis.title.y =  element_text(size = 20, color="#000000"),
    axis.line.x  = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.length  = unit(0.25,"cm"),
    axis.ticks.y  = element_line(size = 1.5),
    legend.position = 'none') + 
  NULL

perm_plot_dataset_b_conpat <- ggplot(dataset_b_long_pvalues_perm, aes(x=freq, y=abs(log10(perm_pvalue)), colour="#000000")) +
  theme_classic()+ my_theme_freq +
  stat_summary(fun.y="mean", geom="bar", colour="#939393", fill = "#939393")+ 
  ylab("Unmedicated vs.\n controls") +
  theme( 
    axis.title.x = element_blank(),
    axis.line.x  = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.ticks.length  = unit(0.25,"cm"),
    axis.ticks.y  = element_line(size = 1.5),    axis.text.y = element_text(size = 18, color="#000000"),
    axis.title.y = element_text(size = 16, color="#000000"),
    legend.position='none')+ 
  geom_hline(yintercept = abs(log10(0.05)))+
  scale_y_continuous(limits = c(0, 3),labels=c("0",  "10e-1","10e-2", "10e-3"))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) 
NULL

perm_plot_dataset_b_pats <- ggplot(dataset_b_long_pvalues_perm_pats_ttest, aes(x=freq, y=abs(log10(perm_pvalue_pats_ttest)), colour="#000000")) +
  theme_classic()+ my_theme_freq +
  stat_summary(fun.y="mean", geom="bar", colour="#939393", fill = "#939393")+ 
  xlab("Frequency [Hz]") + ylab("Unmedicated vs.\n medicated") +
  theme(   
    axis.title.x = element_blank(),
    axis.line.x  = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.ticks.length  = unit(0.25,"cm"),
    axis.ticks.y  = element_line(size = 1.5),
    axis.text.y = element_text(size = 18, color="#000000"),
    axis.title.y = element_text(size = 16, color="#000000"),
    legend.position='none')+ 
  geom_hline(yintercept = abs(log10(0.05)))+
  scale_y_continuous(limits = c(0, 4),labels=c("0",  "10e-1","10e-2", "10e-3", "10e-4"))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) 
NULL

perm_plot_dataset_b_conpat_med <- ggplot(dataset_b_long_pvalues_perm_conpat_med, aes(x=freq, y=abs(log10(perm_pvalue_conpat_med)), colour="#000000")) +
  theme_classic()+ my_theme_freq +
  stat_summary(fun.y="mean", geom="bar", colour="#939393", fill = "#939393")+ 
  xlab("Frequency [Hz]") + ylab("Control vs \n medicated") +
  theme(       
    axis.line  = element_line(size = 10),
    axis.ticks   = element_line(size = 2),
    axis.text.y = element_text(size = 18, color="#000000"),
    axis.text.x = element_text(size = 18, color="#000000"), 
    axis.title.y = element_text(size = 16, color="#000000"),
    axis.title.x = element_text(size = 18, color="#000000"), 
    axis.ticks.length  = unit(0.25,"cm"),
      legend.position='none')+ 
  geom_hline(yintercept = abs(log10(0.05)))+
  scale_y_continuous(limits = c(0, 3),labels=c("0",  "10e-1","10e-2", "10e-3"))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) 



## Dataset C
freq_spec_dataset_c_three <- dataset_c_long_pvalues_perm_three %>%
  ggplot(aes(x = freq, y = mean_densDB, colour = group)) + 
  theme_classic()+my_theme_freq + ylim(0,25)+
  stat_summary(fun.y="mean", geom = 'smooth', se = FALSE)+              
  scale_color_manual(values=color_dataset_c_all, labels = c("Controls - 7days (C)", "Medicated - 7days (C)", "Medicated - 28days (C)")) +  
  ylab("") +
  theme(plot.title  = element_text(hjust = 0.5)) +                                
  theme(      
    axis.text.y  = element_text(size = 20, color="#000000"),
    axis.title.x = element_blank(),
    axis.line.x  = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.length  = unit(0.25,"cm"),
    axis.ticks.y  = element_line(size = 1.5),
    legend.position = 'none') + 
  NULL

perm_plot_dataset_c_conpat <- ggplot(dataset_c_long_pvalues_perm, aes(x=freq, y=abs(log10(perm_pvalue)), colour="#000000")) +
  theme_classic()+ my_theme_freq +
  stat_summary(fun.y="mean", geom="bar", colour="#939393", fill = "#939393")+ 
  xlab("Frequency [Hz]") + ylab("Medicated (7d) vs.\n controls") +
  theme(     
    axis.title.x = element_blank(),
    axis.line.x  = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.ticks.length  = unit(0.25,"cm"),
    axis.ticks.y  = element_line(size = 1.5),
    axis.text.y = element_text(size = 18, color="#000000"),
    axis.title.y = element_text(size = 16, color="#000000"),
    legend.position='none')+ 
  geom_hline(yintercept = abs(log10(0.05)))+
  scale_y_continuous(limits = c(0, 3),labels=c("0",  "10e-1","10e-2", "10e-3"))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) 
NULL

perm_plot_dataset_c_pats <- ggplot(dataset_c_long_pvalues_perm_pats, aes(x=freq, y=abs(log10(perm_pvalue_pats_coin_paired)), colour="#000000")) +
  theme_classic()+ my_theme_freq +
  stat_summary(fun.y="mean", geom="bar", colour="#939393", fill = "#939393")+ 
  xlab("Frequency [Hz]")  + ylab("Medicated (7d) vs.\n medicated (28d)") +
  theme(
    axis.title.x = element_blank(),
    axis.line.x  = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.ticks.length  = unit(0.25,"cm"),
    axis.ticks.y  = element_line(size = 1.5),
    axis.text.y = element_text(size = 18, color="#000000"),
    axis.title.y = element_text(size = 16, color="#000000"),
    legend.position='none')+ 
  geom_hline(yintercept = abs(log10(0.05)))+
  scale_y_continuous(limits = c(0, 3),labels=c("0",  "10e-1","10e-2", "10e-3"))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) 
NULL

perm_plot_dataset_c_conpat_med <- ggplot(dataset_c_long_pvalues_perm_conpat_med, aes(x=freq, y=abs(log10(perm_pvalue_conpat_med)), colour="#000000")) +
  theme_classic()+ my_theme_freq +
  stat_summary(fun.y="mean", geom="bar", colour="#939393", fill = "#939393")+ 
  xlab("Frequency [Hz]") + ylab("Medicated (28d) vs.\n controls") +
  theme(       
    axis.line  = element_line(size = 10),
    axis.ticks   = element_line(size = 2),
    axis.text.y = element_text(size = 18, color="#000000"),
    axis.text.x = element_text(size = 18, color="#000000"), 
    axis.title.y = element_text(size = 16, color="#000000"),
    axis.title.x = element_text(size = 18, color="#000000"), 
    axis.ticks.length  = unit(0.25,"cm"),
    legend.position='none')+ 
  geom_hline(yintercept = abs(log10(0.05)))+
  scale_y_continuous(limits = c(0, 3),labels=c("0",  "10e-1","10e-2", "10e-3"))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) 


freq_plot_dataset_a <- ggarrange(freq_spec_dataset_a, perm_plot_dataset_a, heights = c(1, 1), ncol = 1, nrow = 2)
freq_plot_three_dataset_b  <- ggarrange(freq_spec_dataset_b_three, perm_plot_dataset_b_conpat, 
                                        perm_plot_dataset_b_pats, perm_plot_dataset_b_conpat_med, 
                                        heights = c(3, 0.8, 0.8, 1.2), ncol = 1, nrow = 4)
freq_plot_three_dataset_c  <- ggarrange(freq_spec_dataset_c_three, perm_plot_dataset_c_conpat, 
                                        perm_plot_dataset_c_pats, perm_plot_dataset_c_conpat_med, 
                                        heights = c(3, 0.8, 0.8, 1.2), ncol = 1, nrow = 4)

a <- ggarrange(freq_plot_dataset_a, freq_plot_three_dataset_b,freq_plot_three_dataset_c, ncol = 3, 
               labels =c("A", "B","C"), font.label = list(size = 22))
# ggsave("figs/freq_together_v2.png", a, device = "png", dpi = 300, width = 30, height = 15)
# ggsave("figs/freq_together_v2.pdf", a, device = "pdf", dpi = 300, width = 30, height = 15)


## Sleep spindles and age ----
lm_spindens_age_all      <- lm(spin_dens ~age, dataset_all)
lm_spincount_age_all     <- lm(spin_count ~age, dataset_all)
lm_so_amp_age_all        <- lm(so_amp ~age, dataset_all)
lm_couplingoccur_age_all <- lm(coupling_occur ~age, dataset_all)
bf_lm_spindens_age_all = regressionBF(spin_dens ~ age, data = dataset_all)

## Remove outliers from consolidation 
dataset_all <- rm_outl(dataset_all, "CONS_change")
dataset_comb <- rm_outl(dataset_comb, "CONS_change")

## Split age into young and old
dataset_all$age_split  <- as.numeric(dataset_all$age > median(dataset_all$age))
dataset_comb$age_split <- as.numeric(dataset_comb$age > median(dataset_comb$age))

dataset_comb_conpat_med <- dataset_comb %>%  filter(group == "Controls" | group == "Medicated_long") 
dataset_comb_conpat_med$age_split <- as.numeric(dataset_comb_conpat_med$age > median(dataset_comb_conpat_med$age))

## Depression severity ----
## Severity over time in DATASET B and C
dataset_all_bc      <- dataset_all %>% filter(dataset == "B" | dataset == "C")
dataset_all_bc_hamd <- dataset_all_bc %>% filter(group == "Unmedicated (B)" | group == "Medicated 7d (C)")

dataset_all_bc_hamd_long <- dataset_all_bc_hamd %>% gather(hamd, score, hamd_0, hamd_week1, hamd_week4)
dataset_all_bc_hamd_long %>% dplyr::select(hamd, score, group) %>% group_by(hamd, group) %>%  get_summary_stats(score, type = "mean_se")

## All three datasets
dataset_all_abc_hamd <- dataset_all %>% filter(group == "Medicated (A)" | group == "Unmedicated (B)" | group == "Medicated 7d (C)")
dataset_all_abc_hamd_long <- dataset_all_abc_hamd %>% gather(hamd, score, hamd_0, hamd_week1, hamd_week4)
dataset_all_abc_hamd_long %>% dplyr::select(hamd, score, group) %>% group_by(hamd, group) %>%  get_summary_stats(score, type = "mean_se")

## Take all medicated patients and combine HAMD scores
dataset_comb_med  <- dataset_all %>%  filter(group == "Medicated (A)" | group == "Medicated 7d (B)" | group == "Medicated 7d (C)") 
dataset_all_a_med <- dataset_all %>%  filter(group == "Medicated (A)")
dataset_a_hamd    <- dataset_all_a_med$hamd_0
dataset_all_b_med <- dataset_all %>%  filter(group == "Medicated 7d (B)")
dataset_b_hamd    <- dataset_all_b_med$hamd_week1
dataset_all_c_med <- dataset_all %>%  filter(group == "Medicated 7d (C)")
dataset_c_hamd    <- dataset_all_c_med$hamd_week1

dataset_comb_med$hamds <- c(dataset_a_hamd, dataset_b_hamd, dataset_c_hamd)

## Depression severity and spindle parameters
lm_hamds_spindens_all  <- summary(lm(spin_dens ~ hamds, dataset_comb_med))
lm_hamds_spincount_all <- summary(lm(spin_count ~ hamds, dataset_comb_med))
lm_hamds_cons_all      <- summary(lm(CONS_change ~ hamds, dataset_comb_med))

dataset_comb_med_bf <- dataset_comb_med %>% drop_na(hamds)
bf_lm_spindens_hamd_comb = regressionBF(spin_dens ~ hamds, data = dataset_comb_med_bf)
bf_lm_spincount_hamd_comb = regressionBF(spin_count ~ hamds, data = dataset_comb_med_bf)

## Depression severity and age
lm_hamds_age_all            <- summary(lm(hamds ~ age, dataset_comb_med))
lm_hamds_spin_dens_age_all  <- summary(lm(spin_dens ~ age*hamds, dataset_comb_med))
lm_hamds_spincount_all      <- summary(lm(CONS_change ~ no_episodes, dataset_comb_med))
lm_hamds_cons_all           <- summary(lm(CONS_change ~ hamds, dataset_comb_med))



## Finger tapping ----
## Source again to get clean datasets
source("dataset_a.R")
source("dataset_b.R")

dataset_a <- dataset_a %>% dplyr::rename(BASE = START)

## Check for and remove outliers
dataset_a_rm <- dataset_a
which_outl(dataset_a, "BASE")        #1 control
which_outl(dataset_a, "TRAIN_change")  #1 control
which_outl(dataset_a, "CONS_change") #1 patient
dataset_a_rm <- rm_outl(dataset_a_rm, "BASE")        #1
dataset_a_rm <- rm_outl(dataset_a_rm, "TRAIN_change")  #1
dataset_a_rm <- rm_outl(dataset_a_rm, "CONS_change") #1

dataset_a_long <- dataset_a %>% gather(memory, score, BASE, TRAIN_change, CONS_change)
dataset_a_long %>% dplyr::select(memory, score, group) %>% group_by(group, memory) %>%  get_summary_stats(score, type = "mean_se")
dataset_a_rm_long <- dataset_a_rm %>% gather(memory, score, BASE, TRAIN_change, CONS_change)
dataset_a_rm_long %>% dplyr::select(memory, score, group) %>% group_by(group, memory) %>%  get_summary_stats(score, type = "mean_se")

## Check for and remove outliers
dataset_b_rm <- dataset_b
which_outl(dataset_b, "BASE")        #1 control
which_outl(dataset_b, "TRAIN_change")  #2 controls, 1 unmed, 1 med
which_outl(dataset_b, "CONS_change") #1 med
dataset_b_rm <- rm_outl(dataset_b_rm, "BASE")        #1
dataset_b_rm <- rm_outl(dataset_b_rm, "TRAIN_change")  #1
dataset_b_rm <- rm_outl(dataset_b_rm, "CONS_change") #1
which_outl(dataset_b_rm, "CONS_change") #1 med

dataset_b_long <- dataset_b %>% gather(memory, score, BASE, TRAIN_change, CONS_change)
dataset_b_long %>% dplyr::select(memory, score, group) %>% group_by(group, memory) %>%  get_summary_stats(score, type = "mean_se")
dataset_b_rm_long <- dataset_b_rm %>% gather(memory, score, BASE, TRAIN_change, CONS_change)
dataset_b_rm_long %>% dplyr::select(memory, score, group) %>% group_by(group, memory) %>%  get_summary_stats(score, type = "mean_se")

## Group difference in consolidation
dataset_b_rm_conmed <- dataset_b_rm %>% dplyr::filter(group == "control" | group == "patient_7d")
t.test(CONS_change  ~ group, data = dataset_b_rm_conmed) 

dataset_b_rm_conmed_na   <- dataset_b_rm_conmed %>% drop_na(CONS_change)
ttestBF(formula = CONS_change  ~ group, data = dataset_b_rm_conmed_na)

## remove these lines if you want with outliers:!!!
dataset_b<-dataset_b_rm
dataset_a<-dataset_a_rm
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

colnames_a <- colnames(dataset_a)
colnames_b <- colnames(dataset_b)

common_colnames <- intersect(colnames_a,  colnames_b)

dataset_a <- dataset_a %>% dplyr::select(common_colnames)
dataset_b <- dataset_b %>% dplyr::select(common_colnames)

dataset_a$dataset <- "A"
dataset_b$dataset <- "B"

## To distinguish
dataset_a$group<- revalue(dataset_a$group, c("Controls" = "Controls (A)",
                                             "Patients" = "Medicated (A)"))
dataset_b$group<- revalue(dataset_b$group, c("control" = "Controls (B)",
                                             "patient_unmed" = "Unmedicated (B)",
                                             "patient_7d" = "Medicated 7d (B)"))

## Combine
dataset_all <- rbind(dataset_a, dataset_b)
dataset_all_rm_long <- dataset_all %>% gather(memory, score, BASE, TRAIN_change, CONS_change)
dataset_all_rm_long %>% dplyr::select(memory, score, group) %>% group_by(group, memory) %>%  get_summary_stats(score, type = "mean_se")


## Plots 
rainplot_baseline_dataset_all <- ggplot(dataset_all,aes(x=group,y=BASE,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=BASE),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Baseline correctly tapped')+ xlab('')+ my_theme + raincloud_theme + 
  guides(fill=FALSE,colour=FALSE)+ 
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_fill_manual(values = color_dataset_beh) +  scale_colour_manual(values = color_dataset_beh) +
  geom_signif(annotation=formatC("***", digits=1), y_position=15,xmin=1.2, xmax=2.2, 
              textsize = 6, colour="BLACK")

rainplot_training_dataset_all <- ggplot(dataset_all,aes(x=group,y=TRAIN_change*100,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=TRAIN_change*100),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Training \npercent change [%]')+ xlab('')+ 
  my_theme +  raincloud_theme + geom_hline(yintercept=0, linetype="dashed")+
  guides(fill=FALSE,colour=FALSE)+  
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_fill_manual(values = color_dataset_beh) + 
  scale_colour_manual(values = color_dataset_beh) 

rainplot_consolidation_dataset_all  <- ggplot(dataset_all,aes(x=group,y=CONS_change*100,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+ 
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=CONS_change*100),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Consolidation \npercent change [%]')+ xlab('')+ my_theme +   raincloud_theme + geom_hline(yintercept=0, linetype="dashed")+
  guides(fill=FALSE,colour=FALSE)+ 
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_fill_manual(values = color_dataset_beh) +  scale_colour_manual(values = color_dataset_beh) +
  geom_signif(annotation=formatC("***", digits=1), y_position=100,xmin=1.2, xmax=2.2, 
              textsize = 6, colour="BLACK")

beh_plot_dataset_all     <- ggarrange(rainplot_baseline_dataset_all, rainplot_training_dataset_all, 
                                      rainplot_consolidation_dataset_all, ncol = 1, nrow = 3, labels = c("A", "B","C"))

# ggsave("figs/beh_plot_dataset_all.pdf", beh_plot_dataset_all, device = "pdf", dpi = 300, width = 10, height = 10)
# ggsave("figs/beh_plot_dataset_all.png", beh_plot_dataset_all, device = "png", dpi = 300, width = 10, height = 10)


