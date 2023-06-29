library(tidyverse)
library(readxl)

#insert data
####begin####
order_5_16_prot <- readxl::read_xlsx("5-30-23_prot/5_16_23_prot_order.xlsx")
chip_map_5_16_prot <- readxl::read_xlsx("5-30-23_prot/5_16_23_prot_chipmap.xlsx")

df_calser_prot_15 <- read.delim(file="5-30-23_prot/6_2_23_cal_ser_15ms.spots.txt", sep = "\t", 
                                header = TRUE, skip = 15) %>% 
  mutate(exposure = 15) %>%
  right_join(chip_map_5_16_prot) %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  filter(!type == 11)
df_calser_prot_20 <- read.delim(file="5-30-23_prot/6_2_23_cal_ser_20ms.spots_3.txt", sep = "\t", 
                                header = TRUE, skip = 15) %>% 
  mutate(exposure = 20) %>%
  right_join(chip_map_5_16_prot) %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  filter(!type == 11)
df_calser_prot_30 <- read.delim(file="5-30-23_prot/6_2_23_cal_ser_30ms.spots.txt", sep = "\t", 
                                header = TRUE, skip = 15) %>% 
  mutate(exposure = 30) %>%
  right_join(chip_map_5_16_prot) %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  filter(!type == 11)
df_calser_prot_100 <- read.delim(file="5-30-23_prot/6_2_23_cal_ser_100ms.spots.txt", sep = "\t", 
                                header = TRUE, skip = 15) %>% 
  mutate(exposure = 100) %>%
  right_join(chip_map_5_16_prot) %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  filter(!type == 11)
df_vicser_prot_15 <- read.delim(file="5-30-23_prot/6_2_23_vic_ser_15ms.spots.txt", sep = "\t", 
                                 header = TRUE, skip = 15) %>% 
  mutate(exposure = 15) %>%
  right_join(chip_map_5_16_prot) %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  filter(!type == 11)
df_vicser_prot_30 <- read.delim(file="5-30-23_prot/6_2_23_vic_ser_30ms.spots.txt", sep = "\t", 
                                header = TRUE, skip = 15) %>% 
  mutate(exposure = 30) %>%
  right_join(chip_map_5_16_prot) %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  filter(!type == 11)
df_vicser_prot_100 <- read.delim(file="5-30-23_prot/6_2_23_vic_ser_100ms.spots.txt", sep = "\t", 
                                header = TRUE, skip = 15) %>% 
  mutate(exposure = 100) %>%
  right_join(chip_map_5_16_prot) %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  filter(!type == 11)
df_vicser_prot_300 <- read.delim(file="5-30-23_prot/6_2_23_vic_ser_300ms.spots.txt", sep = "\t", 
                                header = TRUE, skip = 15) %>% 
  mutate(exposure = 300) %>%
  right_join(chip_map_5_16_prot) %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  filter(!type == 11)

full_prot_ser <- rbind(df_calser_prot_15, df_calser_prot_20,
                       df_calser_prot_30, df_calser_prot_100,
                       df_vicser_prot_15, df_vicser_prot_30,
                       df_vicser_prot_100,df_vicser_prot_300)
ser_prot_5_16 <- order_5_16_prot  %>% right_join(full_prot_ser)
ser_prot_5_16_cut <- ser_prot_5_16 %>%
  select(thickness_median_total, protein, uL, concentration,
         Main_Column)

sum_ser_prot_func <- function(df) {
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
    mutate(ctrl_change = FITC_change - mean(FITC_change[serum == "con" & Main_Column == 1]))
  
  return(df_info)
}
final_prot_ser <- sum_ser_prot_func(ser_prot_5_16)
final_prot_ser$dilution <- as.numeric(final_prot_ser$dilution)

####end####

my_eq <- y~((Bmax*x)/(KD+x))
the_dat <- final_prot_ser %>% filter(serum %in% c("cal","con")) %>%
  filter(protein == "Cal") 
plot <- ggplot(final_prot_ser %>% filter(serum %in% c("cal","con")) %>%
                 filter(protein == "Cal"),
               aes(x=as.numeric(dilution), y=ctrl_change))+
  geom_point(aes(color = protein), size = 8)+
  facet_wrap(~as.factor(concentration)) +
  theme_bw()+ ylim(-2,45)+
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
tiff("5-30-23_prot/cal_prot_ser.tiff", units="in", width=18, height=11, res=300)
print(plot)
dev.off()
