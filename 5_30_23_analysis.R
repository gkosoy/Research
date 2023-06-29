library(tidyverse)
library(ggpmisc)

#virus control work
####begin####
order_5_16 <- readxl::read_xlsx("5-30-23/5-16-23_chip_order.xlsx")
chip_map_5_16 <- readxl::read_xlsx("5-30-23/5-16-23_vir_chipmap.xlsx")

df_5_16_vir_75 <- read.delim(file="5-30-23/5-16-23_vir_75ms.spots.txt", sep = "\t", 
                              header = TRUE, skip = 15) %>% 
  mutate(exposure = 75) %>%
  right_join(chip_map_5_16) %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  filter(!type == 11)
df_5_16_vir_200 <- read.delim(file="5-30-23/5-16-23_vir_200ms.spots.txt", sep = "\t", 
                             header = TRUE, skip = 15) %>% 
  mutate(exposure = 200) %>%
  right_join(chip_map_5_16) %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  filter(!type == 11)

vir_con_5_16 <- rbind(df_5_16_vir_75, df_5_16_vir_200)
vir_5_16 <- order_5_16  %>% right_join(vir_con_5_16)

vir_sum_func <- function(df) {
  df_info <- df %>% 
    group_by(virus, Main_Column, protein, concentration) %>%
    summarize(mean_thickness = mean(thickness_median_total),
              std_thickness = sd(thickness_median_total),
              count = n(),
              sem_thickness = std_thickness/(sqrt(count)),
              .groups = 'drop') %>%
    group_by(virus, Main_Column) %>%
    mutate(FITC_change = mean_thickness - (mean_thickness %>% 
                                             keep(protein == "FITC Lo"))) %>%
    group_by(protein, concentration) %>%
    mutate(ctrl_change = FITC_change - mean(FITC_change[virus == "con" & Main_Column == 2]))
  
  return(df_info)
}

vir_final_5_16 <- vir_sum_func(vir_5_16)

plot <- ggplot(data = vir_final_5_16 %>%
         filter(protein %in% c("perth","FITC Lo", "C6", "D8")), 
         aes(x=virus, y=ctrl_change))+
  facet_wrap(~protein)+
  geom_point(aes(color = as.factor(concentration)), size = 7) +
  theme_bw()+xlab("virus")+ylab("Thickness")+
  theme(text = element_text(size = 35))+
  labs(col = "concentration")+
  scale_x_discrete(name ="virus", 
                   labels=c("H1N1","Control","H3N2"))
tiff("5-30-23/vir_con_byprot5_30.tiff", units="in", width=18, height=11, res=300)
print(plot)
dev.off()

plot <- ggplot(vir_final_5_16, 
               aes(x=virus, y=ctrl_change, 
                   color = as.factor(protein)))+
  geom_boxplot(size = 5, position = "dodge")+
  theme(text = element_text(size = 25),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 6))
tiff("5-30-23/all_prot_6_1_virus.tiff", units="in", width=18, height=11, res=300)
print(plot)
dev.off()

####end####

#victoria serum
###begin####
order_5_16 <- readxl::read_xlsx("5-30-23/5-16-23_chip_order.xlsx")
chip_map_5_16 <- readxl::read_xlsx("5-30-23/5-16-23_vir_chipmap.xlsx")

df_5_16_vicser_15 <- read.delim(file="5-30-23/5_16_23_vic_ser_15ms.spots_2.txt", sep = "\t", 
                             header = TRUE, skip = 15) %>% 
  mutate(exposure = 15) %>%
  right_join(chip_map_5_16) %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  filter(!type == 11)
df_5_16_vicser_20 <- read.delim(file="5-30-23/5_16_23_vic_ser_20ms.spots_2.txt", sep = "\t", 
                                header = TRUE, skip = 15) %>% 
  mutate(exposure = 20) %>%
  right_join(chip_map_5_16) %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  filter(!type == 11)
df_5_16_vicser_40 <- read.delim(file="5-30-23/5_16_23_vic_ser_40ms.spots_2.txt", sep = "\t", 
                                header = TRUE, skip = 15) %>% 
  mutate(exposure = 40) %>%
  right_join(chip_map_5_16) %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  filter(!type == 11)

full_vic_ser <- rbind(df_5_16_vicser_15, df_5_16_vicser_20, 
                      df_5_16_vicser_40)
ser_vir_5_16 <- order_5_16  %>% right_join(full_vic_ser) %>%
  filter(!exposure %in% c(75, 200))
ser_vir_5_16_cut <- ser_vir_5_16 %>%
  select(thickness_median_total, protein, uL, concentration,
         Main_Column)

sum_ser_func <- function(df) {
  df_info <- df %>% 
    group_by(serum, dilution, Main_Column, protein, concentration) %>%
    summarize(mean_thickness = mean(thickness_median_total),
              std_thickness = sd(thickness_median_total),
              count = n(),
              sem_thickness = std_thickness/(sqrt(count)),
              .groups = 'drop') %>%
    group_by(serum, Main_Column, dilution) %>%
    mutate(FITC_change = mean_thickness - (mean_thickness %>% 
                                             keep(protein == "FITC Lo"))) %>%
    group_by(protein, concentration) %>%
    mutate(ctrl_change = FITC_change - mean(FITC_change[serum == "con" & Main_Column == 3]))
  
  return(df_info)
}
final_vic_ser <- sum_ser_func(ser_vir_5_16)
final_vic_ser$dilution <- as.numeric(final_vic_ser$dilution)

for_aved_vic <- ser_vir_5_16 %>% filter(!protein == "C4") %>%
  filter(virus == "vic")
final_ave_func <- function(df) {
  df_info <- df %>% 
    group_by(serum, dilution, protein, concentration) %>%
    summarize(mean_thickness = mean(thickness_median_total),
              std_thickness = sd(thickness_median_total),
              count = n(),
              sem_thickness = std_thickness/(sqrt(count)),
              .groups = 'drop') %>%
    group_by(serum, dilution) %>%
    mutate(FITC_change = mean_thickness - (mean_thickness %>% 
                                             keep(protein == "FITC Lo")),
           FITC_sd = sqrt((std_thickness^2) + ((std_thickness %>% keep(protein == "FITC Lo"))^2))) %>%
    group_by(protein, concentration) %>%
    mutate(ctrl_change = FITC_change - mean(FITC_change[serum == "con"]),
           ctrl_std = sqrt((FITC_sd^2) + ((FITC_sd[serum == "con"])^2)))
  return(df_info)
}
vic_ser_aved <- final_ave_func(for_aved_vic)
vic_ser_aved$dilution <- as.numeric(vic_ser_aved$dilution)

FITC_local_func <- function(df) {
  df_info <- df %>% 
    group_by(serum, dilution, Main_Column, protein, concentration, Sub_Column) %>%
    summarize(mean_thickness = mean(thickness_median_total),
              std_thickness = sd(thickness_median_total),
              count = n(),
              sem_thickness = std_thickness/(sqrt(count)),
              .groups = 'drop') 
  #  group_by(serum, Main_Column, dilution, Sub_Column) %>%
   # mutate(FITC_change = mean_thickness - (mean_thickness %>% 
    #                                         keep(protein == "FITC"))) %>%
    #group_by(protein, concentration, Sub_Column) %>%
    #mutate(ctrl_change = FITC_change - mean(FITC_change[serum == "con" & Main_Column == 3]))
  
  return(df_info)
}
local_fix_df <- FITC_local_func(ser_vir_5_16)

for_working_plot <- final_vic_ser %>%
  filter(ctrl_change > -7) %>%
  mutate(ctrl_change_2 = case_when(concentration == 250 ~ ctrl_change+2,
                                   concentration == 350 ~ ctrl_change+5,
                                   concentration == 450 ~ ctrl_change+2,
                                   TRUE ~ ctrl_change))
my_eq <- y~((Bmax*x)/(KD+x))
plot <- ggplot(final_vic_ser %>% filter(protein == "anti-IgG"), 
       aes(x=as.numeric(dilution), y=ctrl_change))+
  geom_point(aes(color = as.factor(concentration)), size = 8)+
  facet_wrap(~as.factor(concentration)) +
  theme_bw()+
  xlab("dilution")+ylab("thickness")+
  theme(text=element_text(size = 35),
        legend.title=element_text(size=5))+
  geom_smooth(method="nls", formula=y~(Bmax*x/(Kd + x)), 
              method.args = list(start = list(Bmax=25,Kd=0.005)), se = FALSE, show.legend = TRUE)+
  stat_fit_tidy(method = "nls", 
                method.args = list(formula = my_eq,
                                   start = list(Bmax=25,KD=0.005)),
                label.x = "right",
                label.y = "top",
                aes(label = paste("B[max]~`=`~", signif(after_stat(Bmax_estimate), digits = 2),
                                  #  "%+-%", signif(after_stat(Bmax_se), digits = 2),
                                  "~~~~K[D]~`=`~", signif(after_stat(KD_estimate), digits = 2),
                                  #  "%+-%", signif(after_stat(KD_se), digits = 2),
                                  sep = "")),
                parse = TRUE, size = 10)
tiff("5-30-23/vic_ser_6_1_IgG.tiff", units="in", width=18, height=11, res=300)
print(plot)
dev.off()

plot <- ggplot(vic_ser_aved %>% filter(protein == "perth")%>%
                 filter(serum %in% c("vic","con")), 
               aes(x=dilution, y=ctrl_change, 
                   color=as.factor(concentration))) + 
  geom_point(size = 8)+
  facet_wrap(~as.factor(concentration))+
  geom_errorbar(aes(ymin=(ctrl_change)-ctrl_std, 
                    ymax=(ctrl_change)+ctrl_std),
                position=position_dodge(0.05)) +
  theme_bw()+
  xlab("dilution")+ylab("thickness")+labs(col="concentration")+
  theme(text=element_text(size = 35),
        legend.title=element_text(size=5))+
  geom_smooth(method="nls", formula=y~(Bmax*x/(Kd + x)), 
              method.args = list(start = list(Bmax=25,Kd=0.005)), se = FALSE, show.legend = TRUE)+
  stat_fit_tidy(method = "nls", 
                method.args = list(formula = my_eq,
                                   start = list(Bmax=25,KD=0.005)),
                label.x = "left",
                label.y = 30,
                aes(label = paste("B[max]~`=`~", signif(after_stat(Bmax_estimate), digits = 2),
                                  #  "%+-%", signif(after_stat(Bmax_se), digits = 2),
                                  "~~~~K[D]~`=`~", signif(after_stat(KD_estimate), digits = 2),
                                  #"%+-%", signif(after_stat(KD_se), digits = 2),
                                  sep = "")),
                parse = TRUE, size = 10)
tiff("5-30-23/vic_ser_6_8_2.tiff", units="in", width=18, height=11, res=300)
print(plot)
dev.off()
####end####  

plot <- ggplot(final_cal_ser, 
       aes(x=protein, y=ctrl_change, 
           fill = as.factor(dilution)))+
  geom_bar(stat="identity", position = "dodge")+
  facet_wrap(~concentration)+
  theme(text = element_text(size = 25),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 6))
tiff("5-30-23/all_prot_6_1_cal.tiff", units="in", width=18, height=11, res=300)
print(plot)
dev.off()

#cali serum
####begin####
order_5_16 <- readxl::read_xlsx("5-30-23/5-16-23_chip_order.xlsx")
chip_map_5_16 <- readxl::read_xlsx("5-30-23/5-16-23_vir_chipmap.xlsx")

df_5_16_calser_15 <- read.delim(file="5-30-23/6_1_23_cal_ser_15ms.spots.txt", sep = "\t", 
                             header = TRUE, skip = 15) %>% 
  mutate(exposure = 15) %>%
  right_join(chip_map_5_16) %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  filter(!type == 11)
df_5_16_calser_20 <- read.delim(file="5-30-23/6_1_23_cal_ser_20ms.spots.txt", sep = "\t", 
                                header = TRUE, skip = 15) %>% 
  mutate(exposure = 20) %>%
  right_join(chip_map_5_16) %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  filter(!type == 11)
df_5_16_calser_30 <- read.delim(file="5-30-23/6_1_23_cal_ser_30ms.spots.txt", sep = "\t", 
                                header = TRUE, skip = 15) %>% 
  mutate(exposure = 30) %>%
  right_join(chip_map_5_16) %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  filter(!type == 11)
df_5_16_calser_40 <- read.delim(file="5-30-23/6_1_23_cal_ser_40ms.spots.txt", sep = "\t", 
                                header = TRUE, skip = 15) %>% 
  mutate(exposure = 40) %>%
  right_join(chip_map_5_16) %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  filter(!type == 11)
df_5_16_calser_70 <- read.delim(file="5-30-23/6_1_23_cal_ser_70ms.spots.txt", sep = "\t", 
                                header = TRUE, skip = 15) %>% 
  mutate(exposure = 70) %>%
  right_join(chip_map_5_16) %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  filter(!type == 11)

full_cal_ser <- rbind(df_5_16_calser_15, df_5_16_calser_20, 
                      df_5_16_calser_30, df_5_16_calser_40,
                      df_5_16_calser_70)
ser_cal_5_16 <- order_5_16  %>% right_join(full_cal_ser) %>%
  filter(!exposure %in% c(75, 200))
ser_cal_5_16_cut <- ser_cal_5_16 %>%
  select(thickness_median_total, protein, uL, concentration,
         Main_Column)

sum_ser_func <- function(df) {
  df_info <- df %>% 
    group_by(serum, dilution, Main_Column, protein, concentration) %>%
    summarize(mean_thickness = mean(thickness_median_total),
              std_thickness = sd(thickness_median_total),
              count = n(),
              sem_thickness = std_thickness/(sqrt(count)),
              .groups = 'drop') %>%
    group_by(serum, Main_Column, dilution) %>%
    mutate(FITC_change = mean_thickness - (mean_thickness %>% 
                                             keep(protein == "FITC Lo")),
           FITC_sd = sqrt((std_thickness^2) + ((std_thickness %>% keep(protein == "FITC Lo"))^2))) %>%
    group_by(protein, concentration) %>%
    mutate(ctrl_change = FITC_change - mean(FITC_change[serum == "con" & Main_Column == 3]),
           ctrl_std = sqrt((FITC_sd^2) + ((FITC_sd[serum == "con" & Main_Column == 3])^2)))
  return(df_info)
}
final_cal_ser <- sum_ser_func(ser_cal_5_16)
final_cal_ser$dilution <- as.numeric(final_cal_ser$dilution)

cal_only_ser <- ser_cal_5_16 %>% filter(!serum == "vic")
final_ave_func <- function(df) {
  df_info <- df %>% 
    group_by(serum, dilution, protein, concentration) %>%
    summarize(mean_thickness = mean(thickness_median_total),
              std_thickness = sd(thickness_median_total),
              count = n(),
              sem_thickness = std_thickness/(sqrt(count)),
              .groups = 'drop') %>%
    group_by(serum, dilution) %>%
    mutate(FITC_change = mean_thickness - (mean_thickness %>% 
                                             keep(protein == "FITC Lo")),
           FITC_sd = sqrt((std_thickness^2) + ((std_thickness %>% keep(protein == "FITC Lo"))^2))) %>%
    group_by(protein, concentration) %>%
    mutate(ctrl_change = FITC_change - mean(FITC_change[serum == "con"]),
           ctrl_std = sqrt((FITC_sd^2) + ((FITC_sd[serum == "con"])^2)))
  return(df_info)
}
cal_vir_aved <- final_ave_func(cal_only_ser) 
cal_vir_aved$dilution <- as.numeric(cal_vir_aved$dilution)

my_eq <- y~((Bmax*x)/(KD+x))
plot <- ggplot(final_cal_ser %>% filter(protein == "C6")%>%
                 filter(serum == "cali"), 
               aes(x=as.numeric(dilution), y=ctrl_change+10))+
  geom_point(aes(color = as.factor(concentration)), size = 8)+
  facet_wrap(~as.factor(concentration)) +
  theme_bw()+
  xlab("dilution")+ylab("thickness")+
  theme(text=element_text(size = 35),
        legend.title=element_text(size=5))+
 geom_smooth(method="nls", formula=y~(Bmax*x/(Kd + x)), 
            method.args = list(start = list(Bmax=25,Kd=0.005)), se = FALSE, show.legend = TRUE)+
 stat_fit_tidy(method = "nls", 
             method.args = list(formula = my_eq,
                               start = list(Bmax=25,KD=0.005)),
           label.x = "right",
          label.y = "top",
        aes(label = paste("B[max]~`=`~", signif(after_stat(Bmax_estimate), digits = 2),
                         #  "%+-%", signif(after_stat(Bmax_se), digits = 2),
                        "~~~~K[D]~`=`~", signif(after_stat(KD_estimate), digits = 2),
                        #"%+-%", signif(after_stat(KD_se), digits = 2),
                       sep = "")),
    parse = TRUE, size = 10)
tiff("5-30-23/cal_ser_C6_6_1.tiff", units="in", width=18, height=11, res=300)
print(plot)
dev.off()

plot <- ggplot(cal_vir_aved %>% filter(protein == "C6")%>%
         filter(serum %in% c("cali","con")), 
       aes(x=dilution, y=ctrl_change, 
                         color=as.factor(concentration))) + 
  geom_point(size = 8)+
  facet_wrap(~concentration)+
  geom_errorbar(aes(ymin=(ctrl_change)-ctrl_std, 
                    ymax=(ctrl_change)+ctrl_std),
                position=position_dodge(0.05)) +
  theme_bw()+
  xlab("dilution")+ylab("thickness")+labs(col="concentration")+
  theme(text=element_text(size = 35),
        legend.title=element_text(size=5))+
  geom_smooth(method="nls", formula=y~(Bmax*x/(Kd + x)), 
              method.args = list(start = list(Bmax=25,Kd=0.005)), se = FALSE, show.legend = TRUE)+
  stat_fit_tidy(method = "nls", 
                method.args = list(formula = my_eq,
                                   start = list(Bmax=25,KD=0.005)),
                label.x = "right",
                label.y = "top",
                aes(label = paste("B[max]~`=`~", signif(after_stat(Bmax_estimate), digits = 2),
                                  #  "%+-%", signif(after_stat(Bmax_se), digits = 2),
                                  "~~~~K[D]~`=`~", signif(after_stat(KD_estimate), digits = 2),
                                  #"%+-%", signif(after_stat(KD_se), digits = 2),
                                  sep = "")),
                parse = TRUE, size = 10)
tiff("5-30-23/cal_ser_C6_6_9.tiff", units="in", width=18, height=11, res=300)
print(plot)
dev.off()
####end####