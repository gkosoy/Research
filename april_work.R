library(tidyverse)
library(ggpmisc)
library(ggpubr)
#dupi work
####begin####
data_for_plot <- final_total_2 %>% filter(!is.na(asthma)) %>% 
  filter(!(ID %in% c(1,9,12))) %>% filter(visit == "B") 

proteins <- c("SARS-CoV-2 S1", "Influenza A California H1N1 07-2009",
              "SARS-CoV RBD", "Influenza B Massachusetts", 
              "CoV-229E Spike (S1+S2)", "Influenza A Texas H3N2",
              "SARS-CoV-2 Spike (S1+S2)", "SARS-CoV-2 S1 D614G mutation",
              "SARS-CoV-2 RBD")
proteins <- c("SARS-CoV-2 S1 D614G mutation")
data_for_plot$time_draw <- factor(data_for_plot$time_draw,      # Reordering group factor levels
                         levels = c("Early", "Mid", "Late"))

for (prot in proteins){
plot <- ggplot(data_for_plot %>% filter(Protein.Name == prot) %>%
                 filter(!is.na(dupi)),
               aes(x = dupi, y = thick_cytc_LOD, color = dupi)) +
  geom_boxplot() + geom_point(size = 6) + theme_bw()+
  ggtitle(prot)+
  stat_compare_means(label = "p.signif", hide.ns = FALSE,size = 15,
                     method = "t.test", label.x = 1.5,
                     symnum.args = list(cutpoints = c(0, 0.001, 0.05, 0.1, 1), 
                                        symbols = c("****", "**", "*", "ns"))) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 45)) 
tiff(paste("4-14-23_work/",prot,"dupi_vs_else_B_box.tiff"), units="in", width=20, height=15, res=300)
print(plot)
dev.off()
}




#for scatterplot time course
protein <- "SARS-CoV-2 S1"
for (prot in proteins){
for_plot <- final_total_2 %>% mutate(day_draw = day_draw-30)
for_plot <- unite(for_plot, col='day', c('day_draw', 'day_post_second_dose'), 
                  sep='', na.rm = TRUE)
for_plot <- for_plot %>% filter(Protein.Name == prot)%>% 
  filter(!(is.na(time_draw))) %>%
  filter(!(is.na(dupi))) #%>%
#filter(visit == "B") %>% filter(!(as.numeric(day) > 120))
for_plot_dupi <- for_plot %>% filter(!asthma == 'no') %>%
  filter(!(ID %in% c("9", "12", "1"))) %>% filter(visit == "B")
plot <- ggplot(for_plot_dupi)+
  geom_point(aes(as.numeric(day), thick_cytc_LOD, color = dupi), 
             size = 12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle(prot)+xlab("day")+ylab("thickness")+
  theme_bw()+theme(text = element_text(size = 45))+
  stat_poly_line(aes(as.numeric(day), thick_cytc_LOD, group = dupi,
                     color = dupi),
                 size = 5, level = 0.9)
tiff(paste("4-14-23_work/",prot,"dupi_day_thick_B_0.9.tiff"), units="in", width=15, height=15, res=300)
print(plot)
dev.off()
}

plot <- ggplot(data_for_plot %>% filter(Protein.Name %in% proteins) %>%
                 filter(!is.na(dupi)) %>% filter(!biologic == "none"),
               aes(x = Protein.Name, y = thick_cytc_LOD, fill = dupi)) +
  geom_boxplot() + 
  xlab("Protein Name")+ylab("Thickness")+
  theme_bw()+
  stat_compare_means(aes(group = dupi), method = "t.test")+ 
                  #   label = "p.signif", hide.ns = FALSE,size = 15,
                   #  method = "wilcox.test", label.x = 1.5)+
                   #  symnum.args = list(cutpoints = c(0, 0.001, 0.05, 0.1, 1), 
                    #                    symbols = c("***", "**", "*", "ns"))) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 45),
        axis.text.x = element_text(angle = 90, size = 30, vjust=0.5, hjust=1)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15))
tiff(paste("4-14-23_work/","ttest_p_all_dupi_vs_biologic_B.tiff"), units="in", width=20, height=15, res=300)
print(plot)
dev.off()

####end####

ID <- c("h1","h2","h3","h4","h5","h6","h7","h8","h9","h10",
                  "h11","h12","h13","h14","h15","h16","h17","h18",
                  "h19","h20","h21","h22","h23","h24","h25","h26","h27",
                  "h28","h29","h30")
vaccine_h <- c("pfizer","pfizer","pfizer","pfizer","pfizer","pfizer",
                   "pfizer","pfizer","pfizer","AstraZeneca","pfizer",
                   "pfizer","pfizer","pfizer","pfizer","pfizer","pfizer",
                   "pfizer", "moderna","moderna","pfizer","", "pfizer",
                   "pfizer","pfizer","moderna","moderna","","pfizer","pfizer")
vaccine_healthy_info <- data.frame(ID, vaccine_h)
plot_df <- final_total_2 %>% left_join(vaccine_healthy_info) 
write.csv(plot_df, "CAVA_May_2023_2.csv")
#asthma vs healthy
####begin####
data_for_plot <- final_total_2 %>% filter(!is.na(asthma)) %>% 
  filter(!(ID %in% c(1,9,12))) %>% filter(visit == "B") 

proteins <- c("SARS-CoV-2 S1", 
              "SARS-CoV RBD", "Influenza B Malaysia",
              "SARS-CoV-2 S1 D614G mutation","SARS-CoV-2 RBD")
proteins <- c("SARS-CoV-2 RBD")
data_for_plot$time_draw <- factor(data_for_plot$time_draw,      # Reordering group factor levels
                                  levels = c("Early", "Mid", "Late"))

for (prot in proteins){
  plot <- ggplot(data_for_plot %>% filter(Protein.Name == prot) %>%
                   filter(!is.na(time_draw)),
                 aes(x = asthma, y = thick_cytc_LOD, color = asthma)) +
    geom_boxplot() +
    geom_point(size = 6) +
    theme_bw()+
    ggtitle(prot)+
    facet_wrap(.~time_draw)+
    stat_compare_means(label = "p.signif", hide.ns = FALSE, 
                       size = 15,
                       method = "t.test", label.x = 1.5,
                       symnum.args = list(cutpoints = c(0, 0.001, 0.05, 0.1, 1), 
                                          symbols = c("****", "**", "*", "ns"))) +
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = 45)) 
  tiff(paste("4-14-23_work/",prot,"asthma_vs_not_B_box.tiff"), units="in", width=20, height=15, res=300)
  print(plot)
  dev.off()
}

for (prot in proteins){
for_plot <- final_total_2 %>% mutate(day_draw = day_draw-30)
for_plot <- unite(for_plot, col='day', c('day_draw', 'day_post_second_dose'), 
                  sep='', na.rm = TRUE)
for_plot <- for_plot %>% filter(Protein.Name == prot)%>% 
  filter(!(is.na(time_draw))) %>%
  filter(!(is.na(asthma))) 
for_plot_asthma <- for_plot %>%
  filter(!(ID %in% c("9", "12", "1"))) %>% filter(visit == "B")
plot <- ggplot(for_plot_asthma)+
  geom_point(aes(as.numeric(day), thick_cytc_LOD, color = asthma), 
             size = 12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle(prot)+xlab("day")+ylab("thickness")+
  theme_bw()+theme(text = element_text(size = 45))+
  stat_poly_line(aes(as.numeric(day), thick_cytc_LOD, group = asthma,
                     color = asthma),
                 size = 5, level = 0.9)
tiff(paste("4-14-23_work/",prot,"asthma_over_time_B_0.9.tiff"), units="in", width=15, height=15, res=300)
print(plot)
dev.off()
}

plot <- ggplot(data_for_plot %>% filter(Protein.Name %in% proteins) %>%
                 filter(!is.na(asthma)),
               aes(x = Protein.Name, y = thick_cytc_LOD, fill = asthma)) +
  facet_wrap(~time_draw)+
  geom_boxplot() + 
  xlab("Protein Name")+ylab("Thickness")+
  theme_bw()+
  stat_compare_means(aes(group = asthma), label = "p.signif", hide.ns = FALSE,size = 15,
                     method = "t.test", label.x = 1.5,
                     symnum.args = list(cutpoints = c(0, 0.001, 0.05, 0.1, 1), 
                                        symbols = c("***", "**", "*", "ns"))) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 45),
        axis.text.x = element_text(angle = 90, size = 30, vjust=0.5, hjust=1)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 18))
tiff(paste("may_work/","asthma_vs_not_time.tiff"), units="in", width=30, height=15, res=300)
print(plot)
dev.off()

data_for_plot <- unite(data_for_plot, col='day', c('day_draw', 'day_post_second_dose'), 
                       sep='', na.rm = TRUE)

plot <- ggplot(data_for_plot %>% filter(Protein.Name == "SARS-CoV-2 RBD") %>%
                 filter(!is.na(asthma)) %>% filter(!is.na(time_draw)),
               aes(x = as.numeric(day), y = thick_cytc_LOD, color = asthma)) +
  facet_wrap(~time_draw)+
  geom_point(size = 7) + 
  xlab("day of draw (post second dose)")+
  ylab("thickness")+
  theme(text = element_text(size = 45),
        axis.text.x = element_text(size = 25))
tiff(paste("may_work/","asthma_vs_not_days.tiff"), units="in", width=30, height=15, res=300)
print(plot)
dev.off()
####end####

#dupi vs healthy
####begin####
data_for_plot <- final_total_2  %>% filter(!(ID %in% c(1,9,12))) %>% 
  filter(visit == "B") #%>% filter(!biologic %in% c("X","N","F","none"))
  
data_for_plot <- unite(data_for_plot, col='day', c('day_draw', 'day_post_second_dose'), 
                  sep='', na.rm = TRUE)

proteins <- c("SARS-CoV-2 S1", "Influenza A California H1N1 07-2009",
              "SARS-CoV RBD", "Influenza B Massachusetts", 
              "CoV-229E Spike (S1+S2)", "Influenza A Texas H3N2",
              "SARS-CoV-2 Spike (S1+S2)", "SARS-CoV-2 S1 D614G mutation",
              "SARS-CoV-2 RBD")
proteins <- c("SARS-CoV-2 S1 D614G mutation")
data_for_plot$time_draw <- factor(data_for_plot$time_draw,      # Reordering group factor levels
                                  levels = c("Early", "Mid", "Late"))
df_try <- data_for_plot %>% filter(Protein.Name == "SARS-CoV-2 RBD") %>%
  filter(!is.na(asthma)) %>% filter(as.numeric(day) > 37 & as.numeric(day) < 91)

plot <- ggplot(data_for_plot %>% filter(Protein.Name == "SARS-CoV-2 RBD") %>%
         filter(!is.na(asthma)) %>% filter(as.numeric(day) > 37 & as.numeric(day) < 91) %>%
           filter(!biologic %in% c("X","N","F")),
       aes(as.numeric(day), thick_cytc_LOD, color = ID))+
 # geom_point(shape = 21, size = 8, stroke = 4, aes(color = dupi))+
  geom_point(size = 7)+
  facet_wrap(~time_draw)+xlab("day")+ylab("thickness")+
  theme_bw()+
  theme(text = element_text(size = 45)) 
tiff(paste("may_work/",prot,"dupi_vs_healthy_days_ID_noelse.tiff"), units="in", width=20, height=15, res=300)
print(plot)
dev.off()

ggplot(data_for_plot %>% filter(Protein.Name == "SARS-CoV-2 RBD") %>%
         filter(!is.na(asthma)) %>% filter(as.numeric(day) > 37 & as.numeric(day) < 91),
       aes(as.numeric(day), thick_cytc_LOD, color = ID, shape = asthma))+
  geom_point(size = 7)+
  facet_wrap(~time_draw)
tiff(paste("may_work/",prot,"dupi_vs_healthy_timed.tiff"), units="in", width=20, height=15, res=300)
print(plot)
dev.off()


ggplot(data_for_plot %>% filter(Protein.Name == "SARS-CoV-2 RBD") %>%
         filter(!is.na(asthma)) %>% filter(asthma == "yes"),
       aes(as.numeric(day), thick_cytc_LOD, color = biologic, shape = asthma))+
  geom_point()

#for (prot in proteins){
  plot <- ggplot(data_for_plot %>% filter(Protein.Name %in% proteins),
                 aes(x = Protein.Name, y = thick_cytc_LOD, color = dupi)) +
    geom_boxplot() + ggtitle("dupi vs healthy wilcox ") +
    geom_point(size = 6,position = position_dodge(width = .75)) +
    theme_bw()+
    stat_compare_means(label = "p.signif", hide.ns = FALSE,size = 15,
                       method = "wilcox", label.x = 1.5,
                       symnum.args = list(cutpoints = c(0, 0.001, 0.05, 0.1, 1), 
                                          symbols = c("****", "**", "*", "ns"))) +
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = 45),
          axis.text.x = element_text(size = 20, angle = 90)) 
  tiff(paste("may_work/","dupi_vs_healthy_all_wilcox.tiff"), units="in", width=20, height=15, res=300)
  print(plot)
  dev.off()
#}
####end####
  
  
#for Sandy poster
  