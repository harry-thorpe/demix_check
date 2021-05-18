#!/usr/bin/env Rscript

library(tidyverse)
library(cowplot)

theme_classic_mod <- theme_classic() +
  theme(axis.text=element_text(colour="black"))

theme_set(theme_classic_mod)

args <- commandArgs(trailingOnly=TRUE)

out_d=args[1]

out_f <- paste(out_d, "/summary_plot.pdf", sep="")

summary_f=paste(out_d, "/clu_out_summary.tsv", sep="")

summary <- read_tsv(summary_f)

summary_l0 <- summary %>%
  filter(level == 0) %>%
  arrange(score, desc(abundance), cluster) %>%
  group_by(ref) %>%
  mutate(l0_idx=1:n()) %>%
  ungroup()

summary_l1 <- summary %>%
  filter(level == 1) %>%
  left_join(., summary_l0 %>% select(cluster, l0_idx), by=c("ref"="cluster")) %>%
  arrange(l0_idx, score, desc(abundance)) %>%
  group_by(ref) %>%
  mutate(l1_idx=1:n()) %>%
  ungroup()

summary_l01 <- bind_rows(summary_l0, summary_l1)

n_l0 <- max(summary_l0$l0_idx)
n_l1 <- max(summary_l1$l1_idx)

pl_l0 <- vector(mode="list", length=n_l0)

l <- vector(mode="list", length=n_l1)
pl_l1 <- list()
for(i in 1:n_l0){
  pl_l1[[i]] <- l
}

for(i in 1:nrow(summary_l01)){

  ref_d=summary_l01$ref_d[i]
  out_dr=summary_l01$out_d[i]
  ref=summary_l01$ref[i]
  clu=summary_l01$cluster[i]
  score=summary_l01$score[i]
  abun=summary_l01$abundance[i]
  level=summary_l01$level[i]
  l0_idx=summary_l01$l0_idx[i]
  l1_idx=summary_l01$l1_idx[i]
  
  dis_ref <- read_tsv(paste(ref_d, "/ref_msh_dis_clu.tsv.gz", sep="")) %>%
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
  
  dis_query_m <- read_tsv(paste(out_dr, "/binned_reads_sketches/", clu, "_msh_dis_clu.tsv.gz", sep="")) %>%
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
           level=level,
           ref=ref,
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
    labs(x="ref-ref distance", y="query-ref distance", title=paste(clu, "\nabundance - ", as.integer(abun*100), "%\nscore - ", score, sep=""))
  
  if(level == 0){
    pl_l0[[l0_idx]] <- p1
  }else if(level == 1){
    pl_l1[[l0_idx]][[l1_idx]] <- p1
  } 
}

l1_rows <- 1
if(n_l1 <= 5){
  l1_rows <- 1
}else if(n_l1 > 5 & n_l1 <= 8){
  l1_rows <- 2
}else if(n_l1 > 9){
  l1_rows <- 3
}
l1_cols <- ceiling(n_l1/l1_rows)

l0_rows <- l1_rows
l0_cols <- l1_rows

h <- n_l0*l0_rows
w <- l0_rows+l1_cols

p_all <- ggdraw() +
  draw_line(x=c(l0_cols/w, l0_cols/w), y=c(0.01, 0.99), linetype="dashed", colour="red") +
  draw_plot_label(label=c("level 0", "level 1"), x=c(0, l0_cols/w), y=c(1, 1)) +
  draw_plot(plot_grid(plotlist=pl_l0, ncol=1), 0, 0, l0_cols/w, 1)

for(i in 1:n_l0){
  p_all <- p_all + draw_plot(plot_grid(plotlist=pl_l1[[i]], nrow=l1_rows, byrow=FALSE), l0_cols/w, 1-(i*l1_rows/h), l1_cols/w, l1_rows/h)
}

pdf(out_f, height=h*3*1.25, width=w*3)

plot(p_all)

dev.off()
