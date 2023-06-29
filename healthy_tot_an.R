
#creating healthy controls dataframe
####begin####
healthy_df <- read.csv("final_LOD_no_cov_no452_cytfix.csv")
healthy_df_2 <- read.csv("final_LOD_no_cov_no452_cytfix_pos.csv")
healthy_df <- healthy_df %>%
  mutate(LOD_corr = corr_cytc_thick - LOD_cytc) %>%
  mutate(thick_cytc_LOD = case_when(LOD_corr > 0 ~ corr_cytc_thick,
                                           LOD_corr < 0 ~ NA_real_)) %>%
  mutate(visit = case_when(time_draw == "s_05" |
                           time_draw == "s_06" | 
                             time_draw == "s_07" | 
                             time_draw == "s_08" ~ "B",
                           TRUE ~ NA_character_)) %>%
  mutate(asthma = "no") %>%
  rename(Protein.Name = protein_name) %>%
  mutate(id = as.character(id)) %>%
  mutate(ID = paste0("h", id)) %>% 
  mutate(time_draw = case_when(day_draw > 22 & day_draw <= 49 ~ "Early",
                               day_draw > 40 & day_draw <= 74 ~ "Mid",
                               day_draw > 75 & day_draw <= 99 ~ "Late",
                               TRUE ~ NA_character_)) 

healthy_df_2 <- healthy_df_2 %>%
  mutate(LOD_corr = corr_cytc_thick - LOD_cytc) %>%
  mutate(thick_cytc_LOD = case_when(LOD_corr > 0 ~ corr_cytc_thick,
                                    LOD_corr < 0 ~ NA_real_)) %>%
  mutate(visit = case_when(time_draw == "s_05" |
                             time_draw == "s_06" | 
                             time_draw == "s_07" | 
                             time_draw == "s_08" ~ "B",
                           TRUE ~ NA_character_)) %>%
  mutate(asthma = "no") %>%
  rename(Protein.Name = protein_name) %>%
  mutate(id = as.character(id)) %>%
  mutate(ID = paste0("h", id)) %>% 
  mutate(time_draw = case_when(day_draw > 22 & day_draw <= 49 ~ "Early",
                               day_draw > 50 & day_draw <= 74 ~ "Mid",
                               day_draw > 75 & day_draw <= 99 ~ "Late",
                               TRUE ~ NA_character_)) 

text_234 <- final_total %>% filter(ID %in% c("h2","h3","h4","VN02","VN03","VN04"))

####end####
final_df_c_2 <- final_df_c %>% mutate(thick_cytc_LOD = FITC_change)
final_total <- healthy_df %>% full_join(final_asthma_LOD)
final_total_try <- healthy_df %>% full_join(final_try_1)
final_try_1 <- final_try_1 %>% mutate(run = 1047)
healthy_df_2$run <- as.character(healthy_df_2$run)
final_total_2 <- healthy_df_2 %>% full_join(final_try_1)


protein <- "SARS-CoV-2 RBD"
protein <- "Influenza A Wisconsin H3N2"
the_day <- "B"
pros <- unique(final_total_try$Protein.Name)
for (protein in pos_prot_list){
plot <- ggplot(final_total_2 %>% filter(Protein.Name == protein) %>%
                 filter(visit == the_day) %>% 
                 filter(!run %in% c(1008, 1013, 1017)) %>%
                 filter(!ID %in% c("VN02","VN03","VN04","VN17","VN27")) %>%
                 filter(!is.na(asthma)), 
               aes(x = asthma, y = thick_cytc_LOD))+
  geom_boxplot() + 
  stat_compare_means(size = 8) +
  geom_point(size=7, aes(color = ID))+ 
  ggtitle(paste(protein, " visit ", the_day)) +
  theme_classic()+ xlab("asthma") + ylab("thickness change") +
  theme(text = element_text(size = 30)) 
tiff(paste("full_df_plots_1_4_22/",the_day,protein,"ID 445-447_cytc_LOD.tiff"), units="in", width=13, height=13, res=300)
print(plot)
dev.off()
}

plot <- ggplot(final_df_c,
               aes(as.factor(ID), Protein.Name)) + 
  geom_tile(aes(fill= FITC_change)) +
  scale_fill_continuous(low = "yellow", high = "red", name = "thickness change") +
 # facet_wrap(. ~ ID, ncol = 6) +
  xlab("") + ylab("") +
  theme(text = element_text(size = 20)) 
tiff("HM_cytc_12-27.tiff", units="in", width=14, height=9, res=300)
print(plot)
dev.off() 

plot <- ggplot(final_total_2 %>%
                 filter(visit == the_day) %>% 
                 filter(!ID %in% c("VN02","VN03","VN04","VN17","VN27")) %>%
                 filter(!is.na(asthma)) %>%
                 filter(!is.na(dupi)) %>%
                 filter(!biologic == "none") %>%
                 filter(!run %in% c(1008,1013,1017)), 
               aes(x = Protein.Name, y = thick_cytc_LOD, fill = dupi))+
  geom_boxplot() + 
 # stat_compare_means(size = 4) +
  ggtitle(paste(" visit ", the_day)) +
  theme_classic()+ xlab("protein") + ylab("thickness change") +
  theme(text = element_text(size = 30)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     hide.ns = TRUE, size = 10,
                     label.y = 32)
tiff("full_df_plots_1_4_22/no_bio_remove_dupi_B_cytc_LOD.tiff", units="in", width=25, height=9, res=300)
print(plot)
dev.off() 


