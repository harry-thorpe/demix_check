#!/usr/bin/env Rscript

library(tidyverse)
library(cowplot)

theme_set(theme_cowplot())

args <- commandArgs(trailingOnly=TRUE)

out_d=args[1]
ref_d=args[2]

out_f <- paste(out_d, "/sample_plot.pdf", sep="")

summary_f=paste(out_d, "/clu_score.tsv", sep="")

summary <- read_tsv(summary_f) %>%
  arrange(desc(score), desc(abundance), cluster) %>%
  mutate(idx=1:n())

n <- max(summary$idx)

pl <- vector(mode="list", length=n)

for(i in 1:nrow(summary)){

  clu=summary$cluster[i]
  abun=summary$abundance[i]
  score=summary$score[i]
  idx=summary$idx[i]
  
  dis_ref <- read_tsv(paste(ref_d, "/ref_msh_dis_clu.tsv", sep="")) %>%
    filter(ref_id != met_id) %>%
    filter(met_cluster == clu) %>%
    mutate(method="mash",
           comparison="ref_ref",
           category=ifelse(met_cluster == ref_cluster, "same cluster", "different cluster"))
  
  dis_ref_r <- dis_ref %>%
    group_by(ref_cluster, met_cluster, category, method, comparison) %>%
    summarise(distance_med=median(distance),
              distance_max=max(distance),
              distance_min=min(distance)) %>%
    ungroup()
  
  thr <- read_tsv(paste(ref_d, "/ref_clu_thr.tsv", sep="")) %>%
    mutate(met_cluster=cluster, ref_cluster=cluster)
  
  dis_query_m <- read_tsv(paste(out_d, "/binned_reads_sketches/", clu, "_msh_dis_clu.tsv", sep="")) %>%
    mutate(method="mash",
           comparison="query_ref",
           category=ifelse(clu == ref_cluster, "same cluster", "different cluster"))
  
  dis_query_m_r <- dis_query_m %>%
    group_by(ref_cluster, category, method, comparison) %>%
    summarise(distance_med=median(distance),
              distance_max=max(distance),
              distance_min=min(distance)) %>%
    ungroup()
  
  thr_p <- thr %>%
    select(ref_cluster, threshold) %>%
    filter(ref_cluster == clu)
  
  dis_ref_query_r <- full_join(dis_ref_r, dis_query_m_r, by="ref_cluster", suffix=c("_ref", "_query")) %>%
    left_join(., thr_p) %>%
    mutate(ref_data=ifelse(is.na(distance_med_ref), "no", "yes"),
           distance_med_ref=ifelse(is.na(distance_med_ref), threshold, distance_med_ref),
           clu=clu)
  
  p1 <- ggplot(dis_ref_query_r, aes(x=distance_med_ref, y=distance_med_query, shape=ref_data, colour=category_query, fill=category_query)) +
    geom_abline(intercept=0, slope=1, linetype="dashed", colour="grey") +
    geom_segment(aes(x=0, xend=threshold, y=threshold, yend=threshold), linetype="dashed", colour="red") +
    geom_segment(aes(x=threshold, xend=threshold, y=0, yend=threshold), linetype="dashed", colour="red") +
    geom_point() +
    geom_segment(aes(x=distance_min_ref, xend=distance_max_ref, y=distance_med_query, yend=distance_med_query)) +
    geom_segment(aes(x=distance_med_ref, xend=distance_med_ref, y=distance_min_query, yend=distance_max_query)) +
    coord_fixed() +
    scale_x_continuous(limits=c(0,max(c(dis_ref_query_r$distance_max_ref, dis_ref_query_r$distance_max_query), na.rm=TRUE))) +
    scale_y_continuous(limits=c(0,max(c(dis_ref_query_r$distance_max_ref, dis_ref_query_r$distance_max_query), na.rm=TRUE))) +
    scale_shape_manual(values=c("no"=1, "yes"=21)) +
    theme(legend.position="none",
          plot.title=element_text(size=9)) +
    labs(x="ref distance", y="query distance", title=paste(clu, "\nabundance - ", as.integer(abun*100), "%\nscore - ", score, sep=""))
  
  pl[[idx]] <- p1

}

rows <- 1
if(n <= 5){
  rows <- 1
}else if(n > 5 & n <= 8){
  rows <- 2
}else if(n > 9){
  rows <- 3
}
cols <- ceiling(n/rows)

p_all <- ggdraw() +
  draw_plot(plot_grid(plotlist=pl, nrow=rows, byrow=FALSE), 0, 0, 1, 1)

pdf(out_f, height=rows*3, width=cols*3)

plot(p_all)

dev.off()

