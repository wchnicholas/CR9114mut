#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(ggrepel)
library(gridExtra)
library(data.table)
require(cowplot)

plot_freq <- function(freq_table, top_clones, WT, HA, HA_name){
  WT_new <- paste(WT,' (WT)',sep='')
  clones_level <- c(WT_new, rev(top_clones))
  freq_table_top_clones <- freq_table %>%
			     filter(Variant %in% top_clones)
  setDT(freq_table_top_clones)
  freq_table_top_clones <- melt(freq_table_top_clones, id='Variant') %>%
			     filter(grepl(HA,variable) | grepl('R0', variable)) %>%
			     mutate(round=str_sub(variable,-1,-1)) %>%
                             mutate(Variant=str_replace(Variant, WT, WT_new)) %>%
                             mutate(Variant=factor(Variant, levels=clones_level))
  colorscale_dark   <- c(brewer.pal(9,"Set1"))
  textsize <- 7
  p <-  ggplot(freq_table_top_clones, aes(round, value/10000, group=Variant, color=Variant)) + 
          geom_line(alpha=0.6) +
          xlab('') +
          ylab('Frequency (%)') +
          scale_color_manual(values=colorscale_dark) +
          ggtitle(paste('Selection target:', HA_name)) +
          theme(plot.title=element_text(size=textsize,face="bold"),
                axis.title=element_text(size=textsize,face="bold"),
                axis.text=element_text(size=textsize,face="bold"),
                axis.text.x=element_text(size=textsize,face="bold",angle=90,hjust=1,vjust=0.5),
                legend.title=element_blank(),
                legend.text=element_text(size=textsize,face="bold"),
                legend.key.size = unit(3, 'mm'),
                legend.position='none') +
          scale_y_continuous(limits=c(0,16)) +
          scale_x_discrete(labels=c('Input','Round 1','Round 2','Round 3')) +
          guides(color=guide_legend(ncol=1, override.aes = list(size=0.5)))
  ggsave(paste('graph/Freq_YDisplay_',HA_name,'.png',sep=''), p, height=1.8, width=2)
  }

freq_table <- read_tsv('data/VariantFreqTable.tsv')
num_top_clones <- 4
WT <- 'NPIFY'
WT_R3_top_clones <- freq_table %>% arrange(WT_R3) %>% tail(num_top_clones) %>% .$Variant %>% unique(c(. ,WT))
I45M_R3_top_clones <- freq_table %>% arrange(I45M_R3) %>% tail(num_top_clones) %>% .$Variant %>% unique(c(. ,WT))
I45T_R3_top_clones <- freq_table %>% arrange(I45T_R3) %>% tail(num_top_clones) %>% .$Variant %>% unique(c(. ,WT))
I45F_R3_top_clones <- freq_table %>% arrange(I45F_R3) %>% tail(num_top_clones) %>% .$Variant %>% unique(c(. ,WT))
plot_freq(freq_table, WT_R3_top_clones, WT, 'WT', 'WT')
plot_freq(freq_table, I45M_R3_top_clones, WT, 'I45M', 'I45M')
plot_freq(freq_table, I45T_R3_top_clones, WT, 'I45T', 'I45T')
plot_freq(freq_table, I45F_R3_top_clones, WT, 'I45F', 'I45F')
