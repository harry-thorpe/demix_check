#!/usr/bin/env Rscript

library(tidyverse)
library(cowplot)

theme_set(theme_cowplot())

args <- commandArgs(trailingOnly=TRUE)

ref_d=args[1]
out_f=paste(ref_d, "/ref_dis_plot.pdf", sep="")
ref=basename(ref_d)

dis_ref <- read_tsv(paste(ref_d, "/ref_msh_dis_clu.tsv.gz", sep="")) %>%
  filter(ref_id != met_id) %>%
  mutate(method="mash",
         comparison="ref_ref",
         category=ifelse(met_cluster == ref_cluster, "same cluster", "different cluster"))

nonsin_clu <- dis_ref %>%
  filter(category == "same cluster") %>%
  select(ref_cluster) %>%
  distinct() %>%
  arrange(ref_cluster)

dis_ref_r <- dis_ref %>%
  group_by(ref_cluster, met_cluster, category, method, comparison) %>%
  summarise(distance_med=median(distance),
            distance_max=max(distance),
            distance_min=min(distance)) %>%
  ungroup() %>%
  mutate(ref_cluster_lab=ifelse(ref_cluster %in% nonsin_clu$ref_cluster, ref_cluster, "singleton clusters"),
         ref_cluster_lab=factor(ref_cluster_lab, levels=c(nonsin_clu$ref_cluster, "singleton clusters")))

thr <- read_tsv(paste(ref_d, "/ref_clu_thr.tsv", sep="")) %>%
  mutate(met_cluster=cluster, ref_cluster=cluster) %>%
  mutate(ref_cluster_lab=ifelse(ref_cluster %in% nonsin_clu$ref_cluster, ref_cluster, "singleton clusters"),
         ref_cluster_lab=factor(ref_cluster_lab, levels=c(nonsin_clu$ref_cluster, "singleton clusters")))

plot_ref <- ggplot(dis_ref_r, aes(x=ref_cluster_lab, y=distance_med, ymin=distance_min, ymax=distance_max, colour=category, group=interaction(ref_cluster, met_cluster))) +
  geom_pointrange(position=position_dodge(width=0.5)) +
  geom_hline(data=thr, aes(yintercept=threshold), linetype="solid", colour="black") +
  facet_grid(.~ref_cluster_lab, space="free", scales="free") +
  theme(legend.position="none",
        strip.text=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1)) +
  labs(x="cluster", y="distance", title=ref)

pdf(out_f, height=10, width=20)

print(plot_ref)

dev.off()
