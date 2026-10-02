### Curto core gene analysis for DB ####
library(tidyverse)
setwd("~/SyncFolder/CA_curto/curto_core_genesDB/")

## read in summary table
sumtab <- read.table("final_summary_table.txt", header = T, sep = "\t")

sumtab <- sumtab %>%
  mutate(true_pcnt = ((unique_curto_seqs/total_curto) * 100)) %>%
  arrange(desc(true_pcnt)) %>%
  mutate(n = 1:nrow(sumtab))

DBgenes <- sumtab %>%
  filter(true_pcnt > 50) %>%
  select(protein)

write.table(DBgenes, "DBgenes.txt", row.names = F, col.names = F, quote = F)

ggplot(sumtab, aes(x = true_pcnt, y = n)) +
  geom_line(color = "black", group = 1) +
  geom_point(color = "black", size = 2, shape = 19) +
  #geom_text(aes(label = ifelse(P_Value < 0.05, "*", "")), color = "black", vjust = 0, size = 8) +
  #ylim(0, 0.5) +
  ggtitle(label = "Gene \"uniqueness\" across all 401 curto genomes") +
  geom_vline(xintercept = 75 , color = "red", linetype = "dashed") +
  geom_vline(xintercept = 50 , color = "red", linetype = "dashed") +
  labs(x = "% of unique sequences for a specific gene", y = "cumulative # of genes") +
  #annotate("text", x = 7.5, y = 0.5, label = substitute(paste("* ", italic("p"), " < 0.05")), size = 4) +
  theme_minimal()

ggsave("gene_uniqueness.png",
       width = 35,
       height = 25,
       units = "cm",
       dpi = 300,
       bg = "white")

ggplot(sumtab, aes(x = unique_curto_seqs, y = n)) +
  geom_line(color = "black", group = 1) +
  geom_point(color = "black", size = 2, shape = 19) +
  #geom_text(aes(label = ifelse(P_Value < 0.05, "*", "")), color = "black", vjust = 0, size = 8) +
  #ylim(0, 0.5) +
  #ggtitle(label = "Correlation between ANI distance and phenotypic distance",
  #        subtitle = "Trait: pH performance") +
  #geom_vline(xintercept = 75 , color = "red", linetype = "dashed") +
  #geom_vline(xintercept = 50 , color = "red", linetype = "dashed") +
  labs(x = "unique sequences for a specific gene", y = "cumulative # of genes") +
  #annotate("text", x = 7.5, y = 0.5, label = substitute(paste("* ", italic("p"), " < 0.05")), size = 4) +
  theme_minimal()



ggplot(sumtab, aes(x = true_pcnt, y = pcnt_curto_total)) +
  #geom_line(color = "black", group = 1) +
  geom_point(color = "black", size = 2, shape = 19) +
  #geom_text(aes(label = ifelse(P_Value < 0.05, "*", "")), color = "black", vjust = 0, size = 8) +
  ylim(0, 100) +
  ggtitle(label = "Gene prevalence and \"uniqueness\"") +
  #geom_vline(xintercept = 75 , color = "red", linetype = "dashed") +
  #geom_vline(xintercept = 50 , color = "red", linetype = "dashed") +
  labs(x = "% of unique sequences for a specific gene", y = "% genomes with a specific gene") +
  #annotate("text", x = 7.5, y = 0.5, label = substitute(paste("* ", italic("p"), " < 0.05")), size = 4) +
  theme_minimal()

ggsave("gene_prevalence.png",
       width = 35,
       height = 25,
       units = "cm",
       dpi = 300,
       bg = "white")
