library(ggpubr)         
library(readr)
library(tidyverse)


dat <- read_csv("C:/Users/finni/Desktop/NE_Calulator/example/Ne_ratio_long_for_R_hominid.csv")%>% 
  mutate(chr_type = factor(chr_type, levels = c("X","Y","MT"), ordered = TRUE))

p <- ggplot(dat,
            aes(x = chr_type, y = ne_ratio, fill = species)) +
  geom_col(position = position_dodge(0.7), width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = ne_ratio , ymax = ne_ratio + se),
                position = position_dodge(0.7), width = 0.2, color = "black")+
  labs(x = "Chromosome type",
       y = "Effective population size ratio ",
       ) +
  
  scale_fill_manual(values = c("Bonobo"     = "#E41A1C",
                               "Chimpanzee" = "#4DAF4A",
                               "Human"      = "#377EB8")) +
  
  theme_classic(base_size = 14) +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
ggsave("C:/Users/finni/Desktop/Ne_ratio_hominid.pdf", p, width = 5, height = 4, dpi = 600)