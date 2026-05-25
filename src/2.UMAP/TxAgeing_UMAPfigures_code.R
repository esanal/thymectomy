library(ggplot2)
library(stringr)
library(RColorBrewer)
library(writexl)
library(ggsignif)
library(readxl)
library(umap)
library(ggstar)
library(ggpubr)
library(dplyr)

### scaling + run UMAP ###
load("UMAP_dataframe_all_2025.RData")
scale_values <- function(x){   (x - min(x))/(max(x)-min(x)) * (0.9 - 0.1) + 0.1 }

UMAP_dataframe_all_s <-UMAP_dataframe_all
rownames(UMAP_dataframe_all_s) <- UMAP_dataframe_all_s[,'X.1']
colnames(UMAP_dataframe_all_s)
UMAP_dataframe_all_s <- UMAP_dataframe_all_s[,-c(1,15:18,35:38,54:57)]

UMAP_dataframe_all_s <- apply(UMAP_dataframe_all_s, 2, scale_values)

UMAP_all <- umap(UMAP_dataframe_all_s)
UMAP_all_df <- UMAP_all$layout %>%
  as.data.frame()%>%
  mutate(X.1=rownames(UMAP_dataframe_all_s))%>%
  inner_join(UMAP_dataframe_all, by="X.1")
colnames(UMAP_all_df)
UMAP_all_df <- UMAP_all_df[,-c(17:20,37:40)]

save(UMAP_all_df, file='dataframewithUMAPcoordinates_All_2025.RData')

#### UMAP figures ####
load("dataframewithUMAPcoordinates_All_2025.RData")

fig1_d1 <- ggplot(UMAP_all_df, aes(x=V1, y=V2)) + geom_star(size=20,aes(starshape=group.wspecials,fill=CD4_naive))+
  ggtitle('% CD4 Naive T cells')+
  theme_classic()+
  theme(plot.title = element_text(hjust = .5, vjust=5, size=60, face='bold'), 
        axis.title.y = element_blank(),axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=40),legend.title =  element_text(size=40),
        legend.key.size = unit(.9,"line"),legend.spacing.y = unit(.7, 'cm'),
        legend.position = c(1.32,.65),
        plot.margin = margin(100,460,60,60))+
  scale_starshape_manual(values=c(13,15,11,1,1,1))+
  scale_fill_gradient(low = "blue3", high="gold", 
                      guide=guide_colourbar(direction='horizontal', title='%', title.position = 'top', title.hjust = 0.5, barwidth = 20, barheight=2,
                                            order=1))+
  guides(starshape = guide_legend(byrow = T, order=2,
                                  override.aes = list(size = 14)))+
  labs(starshape=NULL)
fig1_d1

ggsave(file='UMAP_all_scaled_CD4Naive_2025.pdf', fig1_d1,width = unit(18, 'cm'), height = unit(14, 'cm'))


fig1_d2 <- ggplot(UMAP_all_df, aes(x=V1, y=V2)) + geom_star(size=20,aes(starshape=group.wspecials,fill=CD8_naive))+
  ggtitle('% CD8 Naive T cells')+
  theme_classic()+
  theme(plot.title = element_text(hjust = .2, vjust=5, size=60, face='bold'), 
        axis.title.y = element_blank(),axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=44),legend.title =  element_text(size=40),
        legend.key.size = unit(.9,"line"),legend.spacing.y = unit(.7, 'cm'),
        legend.position = c(1.32,.65),
        plot.margin = margin(100,480,60,60))+
  scale_starshape_manual(values=c(13,15,11,1,1,1))+
  scale_fill_gradient(low = "blue3", high="gold",
                      guide=guide_colourbar(direction='horizontal', title='%', title.position = 'top', title.hjust = 0.5, barwidth = 20, barheight=2,
                                            order=1))+
  guides(starshape = guide_legend(byrow = T, order=2,
                                  override.aes = list(size = 18)))+
  labs(starshape=NULL)
fig1_d2

ggsave(file='UMAP_all_scaled_CD8Naive_2025.pdf', fig1_d2,width = unit(18, 'cm'), height = unit(14, 'cm'))

fig1_e1 <- ggplot(UMAP_all_df, aes(x=V1, y=V2)) + geom_star(size=20,aes(starshape=group.wspecials,fill=CD4_TrueNaive_CD31plus))+
  ggtitle('% CD31+ Naive CD4 T cells')+
  theme_classic()+
  theme(plot.title = element_text(hjust = .5, vjust=5, size=60, face='bold'), 
        axis.title.y = element_blank(),axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=40),legend.title =  element_text(size=40),
        legend.key.size = unit(.9,"line"),legend.spacing.y = unit(.7, 'cm'),
        legend.position = c(1.32,.65),
        plot.margin = margin(100,460,60,60))+
  scale_starshape_manual(values=c(13,15,11,1,1,1))+
  scale_fill_gradient(low = "blue3", high="gold", 
                      guide=guide_colourbar(direction='horizontal', title='%', title.position = 'top', title.hjust = 0.5, barwidth = 20, barheight=2,
                                            order=1))+
  guides(starshape = guide_legend(byrow = T, order=2,
                                  override.aes = list(size = 14)))+
  labs(starshape=NULL)
fig1_e1

ggsave(file='UMAP_all_scaled_CD31CD4Naive_2025.pdf', fig1_e1,width = unit(18, 'cm'), height = unit(14, 'cm'))

fig1_e2 <- ggplot(UMAP_all_df, aes(x=V1, y=V2)) + geom_star(size=20,aes(starshape=group.wspecials,fill=CD4_CD31plus_IL8))+
  ggtitle('% IL8+ in CD31+ Naive CD4 T cells')+
  theme_classic()+
  theme(plot.title = element_text(hjust = .2, vjust=5, size=60, face='bold'), 
        axis.title.y = element_blank(),axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=40),legend.title =  element_text(size=40),
        legend.key.size = unit(.9,"line"),legend.spacing.y = unit(.7, 'cm'),
        legend.position = c(1.32,.65),
        plot.margin = margin(100,460,60,60))+
  scale_starshape_manual(values=c(13,15,11,1,1,1))+
  scale_fill_gradient(low = "blue3", high="gold",
                      guide=guide_colourbar(direction='horizontal', title='%', title.position = 'top', title.hjust = 0.5, barwidth = 20, barheight=2,
                                            order=1))+
  guides(starshape = guide_legend(byrow = T, order=2,
                                  override.aes = list(size = 14)))+
  labs(starshape=NULL)
fig1_e2

ggsave(file='UMAP_all_scaled_IL8CD31CD4Naive_2025.pdf', fig1_e2,width = unit(18, 'cm'), height = unit(14, 'cm'))

fig3_c1 <- ggplot(UMAP_all_df, aes(x=V1, y=V2)) + geom_star(size=20,aes(starshape=group.wspecials,fill=CD8_Immunosenescent))+
  ggtitle(bquote(bold('% Immunosenescent CD8'^+''*' T cells')))+
  theme_classic()+
  theme(plot.title = element_text(hjust = .2, vjust=5, size=60, face='bold'), 
        axis.title.y = element_blank(),axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=44),legend.title =  element_text(size=40),
        legend.key.size = unit(.9,"line"),legend.spacing.y = unit(.7, 'cm'),
        legend.position = c(1.32,.65),
        plot.margin = margin(100,480,60,60))+
  scale_starshape_manual(values=c(13,15,4,1,1,1))+
  scale_fill_gradient(low = "blue3", high="gold",
                      guide=guide_colourbar(direction='horizontal', title='%', title.position = 'top', title.hjust = 0.5, barwidth = 20, barheight=2,
                                            order=1),
                      breaks=c(20,40,60))+
  guides(starshape = guide_legend(byrow = T, order=2,
                                  override.aes = list(size = 18)))+
  labs(starshape=NULL)
fig3_c1

ggsave(file='UMAP_immunosenescent_2025.pdf', fig3_c1,width = unit(18, 'cm'), height = unit(14, 'cm'))
