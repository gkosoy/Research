library(tidyverse)
#creating LOD dataframe from ziva data - add days data
#made month 3 till 98 days for one data point
####begin####
days_df <- readxl::read_xlsx("asthma_data_final_with_COVID_infoxGKx2.xlsx")
days_df <- read_csv("asthma_data_final_with_COVID_infoxGKx2.csv")
days_df <- days_df %>% 
  mutate(time_draw = case_when(day_post_second_dose > 22 & day_post_second_dose <= 49 ~ "Early",
                               day_post_second_dose > 40 & day_post_second_dose <= 74 ~ "Mid",
                               day_post_second_dose > 75 & day_post_second_dose <= 99 ~ "Late",
                               TRUE ~ NA_character_)) %>%
  mutate(Protein.Name = case_when(Protein.Name == 'Influenza A California H1N1 07/2009' ~ 'Influenza A California H1N1 07-2009',
                                  Protein.Name == 'Influenza A California H1N1 04/2009' ~ 'Influenza A California H1N1 04-2009',
                                  TRUE ~ Protein.Name)) %>%
  select(c(1, 7:16))


aves_func_A1 <- function(df){
  df_info <- df %>%
    group_by(Well.Name, Protein.Name) %>%
    summarize(mean_thickness_1 = mean(Raw.Thickness, na.rm = TRUE),
              std_thickness = sd(Raw.Thickness, na.rm = TRUE),
              count = n(),
              .groups = 'drop') %>%
    group_by(Protein.Name) %>%
    mutate(corrected_neg = mean_thickness_1 - mean_thickness_1[Well.Name == "A1"],
           corrected_neg_sd = sqrt((std_thickness^2) + (std_thickness[Well.Name == "A1"]^2))) %>%
    group_by(Well.Name) %>%
    mutate(corr_cytc_thick = corrected_neg - corrected_neg[Protein.Name == "cytC"],
           cytc_sd = sqrt((corrected_neg_sd^2)+(corrected_neg_sd[Protein.Name == "cytC"]^2))) %>% 
    mutate(corr_fitc_thick = corrected_neg - corrected_neg[Protein.Name == ""],
           fitc_sd = sqrt((corrected_neg_sd^2)+(corrected_neg_sd[Protein.Name == ""]^2))) %>% 
    group_by(Protein.Name) %>%
    mutate(LOD = mean(mean_thickness_1[Well.Name == "A1"])+(3*std_thickness[Well.Name == "A1"])) %>%
    mutate(LOD_cytc = mean(corr_cytc_thick[Well.Name == "A1"])+(3*cytc_sd[Well.Name == "A1"]),
           LOD_fitc = mean(corr_fitc_thick[Well.Name == "A1"])+(3*fitc_sd[Well.Name == "A1"])) %>%
    as.data.frame()
  return(df_info)
}
aves_func_G2 <- function(df){
  df_info <- df %>%
    group_by(Well.Name, Protein.Name) %>%
    summarize(mean_thickness_1 = mean(Raw.Thickness, na.rm = TRUE),
              std_thickness = sd(Raw.Thickness, na.rm = TRUE),
              count = n(),
              .groups = 'drop') %>%
    group_by(Protein.Name) %>%
    mutate(corrected_neg = mean_thickness_1 - mean_thickness_1[Well.Name == "G2"],
           corrected_neg_sd = sqrt((std_thickness^2) + (std_thickness[Well.Name == "G2"]^2))) %>%
    group_by(Well.Name) %>%
    mutate(corr_cytc_thick = corrected_neg - corrected_neg[Protein.Name == "cytC"],
           cytc_sd = sqrt((corrected_neg_sd^2)+(corrected_neg_sd[Protein.Name == "cytC"]^2))) %>% 
    mutate(corr_fitc_thick = corrected_neg - corrected_neg[Protein.Name == ""],
           fitc_sd = sqrt((corrected_neg_sd^2)+(corrected_neg_sd[Protein.Name == ""]^2))) %>% 
    group_by(Protein.Name) %>%
    mutate(LOD = mean(mean_thickness_1[Well.Name == "G2"])+(3*std_thickness[Well.Name == "G2"])) %>%
    mutate(LOD_cytc = mean(corr_cytc_thick[Well.Name == "G2"])+(3*cytc_sd[Well.Name == "G2"]),
           LOD_fitc = mean(corr_fitc_thick[Well.Name == "G2"])+(3*fitc_sd[Well.Name == "G2"])) %>%
    as.data.frame()
  return(df_info)
}
aves_func_B1 <- function(df){
  df_info <- df %>%
    group_by(Well.Name, Protein.Name) %>%
    summarize(mean_thickness_1 = mean(Raw.Thickness, na.rm = TRUE),
              std_thickness = sd(Raw.Thickness, na.rm = TRUE),
              count = n(),
              .groups = 'drop') %>%
    group_by(Protein.Name) %>%
    mutate(corrected_neg = mean_thickness_1 - mean_thickness_1[Well.Name == "B1"],
           corrected_neg_sd = sqrt((std_thickness^2) + (std_thickness[Well.Name == "B1"]^2))) %>%
    group_by(Well.Name) %>%
    mutate(corr_cytc_thick = corrected_neg - corrected_neg[Protein.Name == "cytC"],
           cytc_sd = sqrt((corrected_neg_sd^2)+(corrected_neg_sd[Protein.Name == "cytC"]^2))) %>% 
    mutate(corr_fitc_thick = corrected_neg - corrected_neg[Protein.Name == ""],
           fitc_sd = sqrt((corrected_neg_sd^2)+(corrected_neg_sd[Protein.Name == ""]^2))) %>% 
    group_by(Protein.Name) %>%
    mutate(LOD = mean(mean_thickness_1[Well.Name == "B1"])+(3*std_thickness[Well.Name == "B1"])) %>%
    mutate(LOD_cytc = mean(corr_cytc_thick[Well.Name == "B1"])+(3*cytc_sd[Well.Name == "B1"]),
           LOD_fitc = mean(corr_fitc_thick[Well.Name == "B1"])+(3*fitc_sd[Well.Name == "B1"])) %>%
    as.data.frame()
  return(df_info)
}

df_1_LOD <- aves_func_A1(ziva_data_1) %>% right_join(sample_1) %>% mutate(date = "3_11")
df_2_LOD <- aves_func_A1(ziva_data_2) %>% right_join(sample_2) %>% mutate(date = "4_13")
df_3_LOD <- aves_func_G2(ziva_data_3) %>% right_join(sample_3) %>% mutate(date = "4_29")
df_4_LOD <- aves_func_B1(ziva_data_4) %>% right_join(sample_4) %>% mutate(date = "5_11")

all_df_LOD <- bind_rows(df_1_LOD, df_2_LOD, df_3_LOD,
                         df_4_LOD)

df_1_LOD <- combine_ziva_sample_info_all(df = df_1_LOD, sample = sample_1, remove_well = NULL, remove_sample = c("1A", "15F", "Blank_1", "Blank_2", "Blank_3", "Blank_4","Blank_5", "29B"))
df_2_LOD <- combine_ziva_sample_info_all(df = df_2_LOD, sample = sample_2, remove_well = NULL, remove_sample = NULL)
df_3_LOD <- combine_ziva_sample_info_all(df = df_3_LOD, sample = sample_3, remove_well = "A1", remove_sample = NULL)
df_4_LOD <- combine_ziva_sample_info_all(df = df_4_LOD, sample = sample_4, remove_well = c("A2","B2","H1"), remove_sample = NULL)

final_asthma_LOD <- bind_rows(df_1_LOD,df_2_LOD,df_3_LOD,df_4_LOD)

final_asthma_LOD <- final_asthma_LOD %>%
  mutate(Protein.Name = case_when(Protein.Name == 'Influenza A California H1N1 07/2009' ~ 'Influenza A California H1N1 07-2009',
                                  Protein.Name == 'Influenza A California H1N1 04/2009' ~ 'Influenza A California H1N1 04-2009',
                                  TRUE ~ Protein.Name)) 

final_asthma_LOD <- final_asthma_LOD %>% 
  mutate(dupi = case_when(biologic == "D" ~ "Dupi",
                          visit == "con" ~ NA_character_,
                          sample == 'blank' ~ NA_character_,
                          sample == 'Blank' ~ NA_character_,
                          TRUE ~ "else")) %>%
  mutate(LOD_corr = corr_cytc_thick - LOD_cytc) %>%
  mutate(thick_cytc_LOD = case_when(LOD_corr > 0 ~ corr_cytc_thick,
                                    LOD_corr < 0 ~ NA_real_))

final_asthma_LOD <- final_asthma_LOD %>%
  mutate(sample = case_when(sample == "Blank" ~ "Blank_1",
                            TRUE ~ sample))

####end####
days_df <- days_df %>% mutate(ID = case_when(sample == "Pos_1" ~ "pos",
                                             sample == "Pos_2" ~ "pos",
                                             TRUE ~ ID))
days_df$ID <- as.character(days_df$ID)
final_asthma_LOD <- final_asthma_LOD %>% mutate(run = case_when(date == "3_11" ~ "1027",
                                                                date == "4_13" ~ "1047",
                                                                date == "4_29" ~ "1052",
                                                                date == "5_11" ~ "1053",
                                                                TRUE ~ NA_character_)) %>%
  mutate(biologic = case_when(date == "3_11" & visit == "con" ~ NA_character_,
                              TRUE ~ biologic))
final_try_1 <- final_asthma_LOD %>% right_join(days_df)


#difference in average day of sample taken 
#day_draw for healthy controls starts at day 0(from first dose day)
####begin####
looking_days <- final_total_2 %>% filter(visit == "B") %>%
  filter(Protein.Name == "SARS-CoV-2 RBD") %>% 
  mutate(day_draw = day_draw - 28) %>%
  unite(day_col, c(day_draw, day_post_second_dose), na.rm = TRUE)
summary(as.integer(looking_days$day_post_second_dose))

#filter ID9 because it's astra-zeneca
plot <- ggplot(looking_days %>% filter(!is.na(asthma)) %>%
         filter(!(ID == "9" & asthma == "yes")) %>% filter(!run %in% c(1008,1013,1017)), 
       aes(x = asthma, y = as.integer(day_col)))+
  geom_boxplot()+
  theme_classic() +
  stat_summary(fun = "mean", color = "black", size = 3)+
  geom_point(size=7, aes(color = time_draw))+ 
  ggtitle(paste("day post second dose - visit B")) +
  theme_classic()+ xlab("asthma") + ylab("Day") +
  theme(text = element_text(size = 40)) 
tiff("full_df_plots_1_4_22/day_post_2nd-B-only_445-447_mean.tiff", units="in", width=13, height=13, res=300)
print(plot)
dev.off()
####end####

#comparing neg and pos controls healthy to asthma
####begin####
#looking at negative and positive controls for proteins
neg_df <- final_total %>% filter(ID == "h0" | ID == "neg" | sample == "neg ctrl00-0") %>%
  filter(Protein.Name %in% c("SARS-CoV-2 RBD", "SARS-CoV-2 Spike (S1+S2)",
                             "Influenza B Massachusetts", "SARS-CoV-2 S1 D614G mutation",
                             "SARS-CoV-2 S1","Influenza B Malaysia","Influenza B Phuket",
                             "SARS-CoV RBD","Influenza A California H1N1 07-2009"))
pos_df <- final_total_2 %>% filter(ID == "h61" | ID == "h62" | ID == "pos") %>%
  filter(Protein.Name %in% c("SARS-CoV-2 RBD", "SARS-CoV-2 Spike (S1+S2)",
                             "Influenza B Massachusetts", "SARS-CoV-2 S1 D614G mutation",
                             "SARS-CoV-2 S1","Influenza B Malaysia","Influenza B Phuket",
                             "SARS-CoV RBD","Influenza A California H1N1 07-2009"))
pos_df_new <- pos_df %>% 
  filter(sample %in% c("pos ctrl61-0", "Pos_1")) %>%
  mutate(asthma_run = case_when(run == "1027" ~ "yes",
                                run == "1047" ~ "yes",
                                run == "1052" ~ "yes",
                                run == "1053" ~ "yes",
                                TRUE ~ "no")) %>%
  filter(!Protein.Name %in% c("Influenza A California H1N1 07-2009",
                              "Influenza B Massachusetts","SARS-CoV-2 Spike (S1+S2)")) %>%
  group_by(asthma_run, Protein.Name)# %>%
  summarise(ave = mean(thick_cytc_LOD),
            std = sd(thick_cytc_LOD))
           

plot <- ggplot(pos_df_new, aes(x = run, y = thick_cytc_LOD, fill = run))+
  geom_bar(stat = "identity", width = 0.5) +
  facet_wrap(~ Protein.Name)+
  theme_bw()+
  ylab("Thickness")+
  xlab("positive control") +
  theme(text = element_text(size = 30)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
tiff("4-14-23_work/pos_con_by_prot.tiff", units="in", width=24, height=15, res=300)
print(plot)
dev.off()
####end####

pos_prot_list <- c("SARS-CoV-2 RBD", "SARS-CoV-2 Spike (S1+S2)",
                   "Influenza B Massachusetts", "SARS-CoV-2 S1 D614G mutation",
                   "SARS-CoV-2 S1","Influenza B Malaysia","Influenza B Phuket",
                   "SARS-CoV RBD","Influenza A California H1N1 07-2009")
protein <- "SARS-CoV-2 RBD"
protein <- "SARS-CoV-2 S1 D614G mutation"
proteins <- c("SARS-CoV-2 RBD", "SARS-CoV-2 S1", 
              "SARS-CoV RBD", "SARS-CoV-2 S1 D614G mutation")
the_day <- "B"
#for (protein in pos_prot_list){
plot <- ggplot(final_total_2 %>% filter(asthma == "yes") %>% 
                 filter(!biologic == "none") %>%
                 filter(Protein.Name %in% proteins) %>%
                 filter(visit == the_day), aes(x = dupi, y = thick_cytc_LOD))+
  geom_boxplot() + 
  stat_compare_means(size = 8, method = "wilcox.test") +
  geom_point(size=7, aes(color = Protein.Name))+ 
  ggtitle(paste("covid proteins - visit ", the_day)) +
  theme_classic()+ xlab("") + ylab("thickness change") +
  theme(text = element_text(size = 30)) +
  ylim(0, 45)
tiff(paste("full_df_plots_1_4_22/", the_day,protein, "none.rm_time_dupi_cytc_LOD.tiff"), units="in", width=13, height=13, res=300)
print(plot)
dev.off()
#}

plot <- ggplot(final_total_2 %>% filter(asthma == "yes") %>% 
                # filter(!biologic == "none") %>%
                 filter(Protein.Name == "SARS-CoV-2 RBD") %>%
                 filter(visit == the_day), aes(x = dupi, y = thick_cytc_LOD))+
  geom_boxplot() + 
  stat_compare_means(size = 13, method = "wilcox.test") +
  geom_point(size=8, aes(color = time_draw))+ 
  ggtitle("SARS-CoV-2 RBD with all asthmatics") +
  theme_classic()+ xlab("") + ylab("thickness change") +
  theme(text = element_text(size = 45))+
  ylim(0,45) +
  theme(plot.margin=unit(c(2,2,0,2.5), "cm"),
        axis.title.y=element_text(vjust=4),
        plot.title = element_text(size = 40, vjust = 4))
tiff(paste("time_work/", the_day,"SARS-CoV-2 RBD", "time_dupi_cytc_LOD.tiff"), units="in", width=15, height=13, res=300)
print(plot)
dev.off()
