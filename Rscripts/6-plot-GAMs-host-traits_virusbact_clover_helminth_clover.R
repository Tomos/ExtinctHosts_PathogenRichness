# This code was originally written by Olival et al. (2017)
# and was adapted (lightly) by Liam Shaw 2019 (liam.philip.shaw at gmail dot com)
# for this project.

# See: https://zenodo.org/record/807517 for the original code repository this code was sourced from

# This file is based on: scripts/06-make-Figure02-all-gams.R from that repository
# The difference here is we compare effects for variables for viral richness per host and bacterial richness per host (rather than zoonotic proportion as the second row of plots as in Olival et al.)


# I am grateful to Olival et al. for making their original code available under an MIT License, which also applies here. 
# https://opensource.org/licenses/MIT
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

library(ggplot2)
library(tidyr)
library(purrr)
library(stringi)
library(cowplot)
library(viridis)
library(dplyr)
library(svglite)
library(mgcv)
library(magrittr)
library(car)
library(gratia)
library(mgcv)
library(gratia)
library(gsubfn)
library(grid)
library(countreg)
library(stringr)

set.seed(0)
rm(list = ls())
if(!is.null(dev.list())) dev.off()
base_dir <- "C:/Users/tomos/OneDrive/PhD/Host-parasite-paper/"
source(paste0(base_dir, "Olival-functions2.R"))
SHOW_DEV_EXPL = FALSE

partials_theme = theme(text = element_text(family="Helvetica", size=9),
                       panel.border=element_blank(),
                       panel.background=element_blank(),
                       axis.title.y = element_blank(),
                       panel.grid = element_blank(),
                       axis.ticks.x = element_line(size=0.3),
                       axis.ticks.y = element_blank(),
                       axis.text = element_text(color="black"),
                       axis.title.x = element_text(lineheight = 1.2),
                       legend.position="none"
                       #plot.margin=margin(l=0)
)

blankPlot <- ggplot()+geom_blank(aes(1,1)) +
  cowplot::theme_nothing()

model_dist <- 'tw'
# viruses:

bgam_ne <- readRDS(paste0(base_dir, 'intermediates/noExtinction_clover_virus_models_', model_dist, '.rds'))$model[[1]]
summary(bgam_ne)
anova.gam(bgam_ne)
# Check whether there are relationships between residulas using gam.check():
gam.check(bgam_ne, rep = 1000)
# Plot residuals
appraise(bgam_ne, method = 'simulate')
# Concurvity:
concurvity(bgam_ne, full=T)
con <- concurvity(bgam_ne, full=F)
con_worst <-format(con$worst, scientific = F)
#con_worst[con_worst < 0.8] <- 0

#root_bgam_ne <- rootogram(bgam_ne, style = "hanging")

model_dist <- 'tw'

bgam = readRDS(paste0(base_dir, 'intermediates/clover_virus_models_', model_dist, '.rds'))$model[[1]]
summary(bgam)
# Generates table of the best model's summary stats:
virus_summary <- summary(bgam)
vs_table <- as.data.frame(virus_summary$s.table)
vs_table$p_value_nonsci <- format(vs_table$`p-value`, scientific = F)
vs_table$input_variable <- rownames(vs_table)
vs_table <- vs_table[, c(6,1,2,3,4,5)]
rownames(vs_table) <- NULL
vs_form <- virus_summary$formula
vs_form <- gsub("\"", "", vs_form)
vs_form <- sub("\n", "", vs_form)
vs_form <- sub("\t", "", vs_form)
names(vs_table)[names(vs_table) == 'input_variable'] <- paste0(as.character(vs_form)[2], " ", as.character(vs_form)[1], " ", as.character(vs_form)[3])
vs_table$edf <- formatC(vs_table$edf, format = 'f', digits = 2)
vs_table$F <- formatC(vs_table$F, format = 'f', digits = 2)
write.csv(vs_table, paste0(base_dir, "Figures/clover_virus_bgam_summary.csv"))

#anova.gam(bgam)
aic_scores <- as.data.frame(AIC(bgam_ne, bgam))
write.csv(aic_scores, paste0(base_dir, "Figures/clover_virus_bgam_aic.csv"))

#comp_gam <- anova.gam(bgam_ne, bgam, test = "LRT") #Chisq or F
#write.csv(comp_gam, paste0(base_dir, "Figures/clover_virus_gam_chisq.csv"))


# Check whether there are relationships between residulas using gam.check():
gam.check(bgam, rep = 1000)
# Plot residuals
appraise(bgam, method = 'simulate')
# Concurvity:
concurvity(bgam, full=T)
con <- concurvity(bgam, full=F)
con_worst <-format(con$worst, scientific = F)
#con_worst[con_worst < 0.8] <- 0

#root_bgam <- rootogram(bgam, style = "Hanging")

de_bgam =  get_relative_contribs(bgam) %>%
  mutate(dev_explained = rel_deviance_explained * summary(bgam)$dev.expl) %>%
  mutate(dev_explained = paste0(stri_trim_both(formatC(dev_explained*100, format = "fg", digits=2)), "%"))
#de_bgam <- arrange(de_bgam, desc(as.numeric(sub("%", "", de_bgam$dev_explained))))
print(de_bgam)


# Generates table of the best 5 model's:
de_bgam1 <- de_bgam
bgam1_sum <- as.data.frame(summary(bgam)$s.table)
bgam1_sum_sig <- bgam1_sum[bgam1_sum$`p-value`<=0.05,]
bgam1_names <- gsub("s\\(", "", rownames(bgam1_sum_sig))
bgam1_names <- gsub("\\)", "", bgam1_names)
de_bgam1 =  get_relative_contribs(bgam) %>%
  mutate(dev_explained = rel_deviance_explained * summary(bgam)$dev.expl) %>%
  mutate(dev_explained = paste0(stri_trim_both(formatC(dev_explained*100, format = "fg", digits=2)), "%")) %>%
  filter(term %in% bgam1_names)

bgam2 = readRDS(paste0(base_dir, 'intermediates/clover_virus_models_', model_dist, '.rds'))$model[[2]]
bgam2_sum <- as.data.frame(summary(bgam2)$s.table)
bgam2_sum_sig <- bgam2_sum[bgam2_sum$`p-value`<=0.05,]
bgam2_names <- gsub("s\\(", "", rownames(bgam2_sum_sig))
bgam2_names <- gsub("\\)", "", bgam2_names)
de_bgam2 =  get_relative_contribs(bgam2) %>%
  mutate(dev_explained = rel_deviance_explained * summary(bgam2)$dev.expl) %>%
  mutate(dev_explained = paste0(stri_trim_both(formatC(dev_explained*100, format = "fg", digits=2)), "%")) %>%
  filter(term %in% bgam2_names)

bgam3 = readRDS(paste0(base_dir, 'intermediates/clover_virus_models_', model_dist, '.rds'))$model[[3]]
bgam3_sum <- as.data.frame(summary(bgam3)$s.table)
bgam3_sum_sig <- bgam3_sum[bgam3_sum$`p-value`<=0.05,]
bgam3_names <- gsub("s\\(", "", rownames(bgam3_sum_sig))
bgam3_names <- gsub("\\)", "", bgam3_names)
de_bgam3 =  get_relative_contribs(bgam3) %>%
  mutate(dev_explained = rel_deviance_explained * summary(bgam3)$dev.expl) %>%
  mutate(dev_explained = paste0(stri_trim_both(formatC(dev_explained*100, format = "fg", digits=2)), "%")) %>%
  filter(term %in% bgam3_names)

bgam4 = readRDS(paste0(base_dir, 'intermediates/clover_virus_models_', model_dist, '.rds'))$model[[4]]
bgam4_sum <- as.data.frame(summary(bgam4)$s.table)
bgam4_sum_sig <- bgam4_sum[bgam4_sum$`p-value`<=0.05,]
bgam4_names <- gsub("s\\(", "", rownames(bgam4_sum_sig))
bgam4_names <- gsub("\\)", "", bgam4_names)
de_bgam4 =  get_relative_contribs(bgam4) %>%
  mutate(dev_explained = rel_deviance_explained * summary(bgam4)$dev.expl) %>%
  mutate(dev_explained = paste0(stri_trim_both(formatC(dev_explained*100, format = "fg", digits=2)), "%")) %>%
  filter(term %in% bgam4_names)


bgam5 = readRDS(paste0(base_dir, 'intermediates/clover_virus_models_', model_dist, '.rds'))$model[[5]]
bgam5_sum <- as.data.frame(summary(bgam5)$s.table)
bgam5_sum_sig <- bgam5_sum[bgam5_sum$`p-value`<=0.05,]
bgam5_names <- gsub("s\\(", "", rownames(bgam5_sum_sig))
bgam5_names <- gsub("\\)", "", bgam5_names)
de_bgam5 =  get_relative_contribs(bgam5) %>%
  mutate(dev_explained = rel_deviance_explained * summary(bgam5)$dev.expl) %>%
  mutate(dev_explained = paste0(stri_trim_both(formatC(dev_explained*100, format = "fg", digits=2)), "%")) %>%
  filter(term %in% bgam5_names)

de_bgam1 <- as.data.frame(de_bgam1)
de_bgam1 <- cbind(de_bgam1, rep("1", nrow(de_bgam1)))
colnames(de_bgam1)[4] <- "Model_rank"
de_bgam2 <- as.data.frame(de_bgam2)
de_bgam2 <- cbind(de_bgam2, rep("2", nrow(de_bgam2)))
colnames(de_bgam2)[4] <- "Model_rank"
de_bgam3 <- as.data.frame(de_bgam3)
de_bgam3 <- cbind(de_bgam3, rep("3", nrow(de_bgam3)))
colnames(de_bgam3)[4] <- "Model_rank"
de_bgam4 <- as.data.frame(de_bgam4)
de_bgam4 <- cbind(de_bgam4, rep("4", nrow(de_bgam4)))
colnames(de_bgam4)[4] <- "Model_rank"
de_bgam5 <- as.data.frame(de_bgam5)
de_bgam5 <- cbind(de_bgam5, rep("5", nrow(de_bgam5)))
colnames(de_bgam5)[4] <- "Model_rank"
top5_model <- rbind(de_bgam1, de_bgam2, de_bgam3, de_bgam4, de_bgam5)
top5_model_results <- top5_model %>%
  group_by(term) %>%
  summarize(n = n())
write.csv(top5_model_results, paste0(base_dir, "Figures/virus_top5gam_variable_breakdown.csv"))
rel_deviance_perc <- as.numeric(de_bgam$rel_deviance_explained*100)
de_bgam$rel_deviance_explained_perc <- paste0(formatC(rel_deviance_perc, format = 'f', digits = 2), "%")
de_bgam$deviance_explained_numeric <- as.numeric(sub("%", "", de_bgam$dev_explained))
write.csv(de_bgam, paste0(base_dir, "Figures/clover_virus_gam_deviance.csv"))

ex_v_de <- sum(de_bgam$deviance_explained_numeric[de_bgam$term == "extinct.min.phylo.dist"])

bio_v_de <- sum(de_bgam$deviance_explained_numeric[de_bgam$term == "extinct.min.phylo.dist"],
                de_bgam$deviance_explained_numeric[de_bgam$term == "hOrderCETARTIODACTYLA"],
                de_bgam$deviance_explained_numeric[de_bgam$term == "hOrderEULIPOTYPHLA"],
                de_bgam$deviance_explained_numeric[de_bgam$term == "hOrderPRIMATES"],
                de_bgam$deviance_explained_numeric[de_bgam$term == "hOrderRODENTIA"],
                de_bgam$deviance_explained_numeric[de_bgam$term == "LnAreaHost"],
                de_bgam$deviance_explained_numeric[de_bgam$term == "S40"],
                de_bgam$deviance_explained_numeric[de_bgam$term == "Sympat_n_Phylacine"])

paste0("Relative deviance (of CLOVER richness per host) across all extinction variables: ", ex_v_de, " %")
paste0("Relative deviance (of CLOVER richness per host) for biologically relevant variables: ", bio_v_de, " %")

# All Viruses Plot
binary_vars = c("hOrder")

preds <- predict(bgam, type="iterms", se.fit=TRUE)
intercept <- attributes(preds)$constant

preds = preds %>% map(as.data.frame) %>%
  map(~setNames(., stri_replace_first_regex(names(.) ,"s\\(([^\\)]+)\\)", "$1")))
model_data_vir = bgam$model
gterms = attributes(terms(bgam))
ilfun <- bgam$family$linkinv
lfun <- bgam$family$linkfun

binary_terms = which(stri_detect_regex(names(model_data_vir), paste0("(", paste0(binary_vars, collapse="|"), ")")))

smooth_data_vir = model_data_vir[, -c(binary_terms, gterms$response, gterms$offset)]
smooth_ranges = data.frame(map(smooth_data_vir, ~seq(min(.), max(.), length.out = 100)))

binary_data = model_data_vir[, binary_terms]
binary_ranges = setNames(as.data.frame(diag(length(binary_terms))), names(model_data_vir)[binary_terms])

offset_name = stri_replace_first_regex(names(model_data_vir)[gterms$offset], "offset\\(([^\\)]+)\\)", "$1")

smooth_ranges = cbind(smooth_ranges,
                      setNames(as.data.frame(lapply(c(names(binary_data), offset_name), function(x) rep(0, nrow(smooth_ranges)))), c(names(binary_data), offset_name)))
binary_ranges = cbind(binary_ranges,
                      setNames(as.data.frame(lapply(c(names(smooth_data_vir), offset_name), function(x) rep(0, nrow(binary_ranges)))), c(names(smooth_data_vir), offset_name)))

smooth_preds <- predict(bgam, newdata=smooth_ranges, type="iterms", se.fit=TRUE)  %>% map(as.data.frame) %>%
  map(~setNames(., stri_replace_first_regex(names(.) ,"s\\(([^\\)]+)\\)", "$1")))
binary_preds <- predict(bgam, newdata=binary_ranges, type="iterms", se.fit=TRUE) %>% map(as.data.frame) %>%
  map(~setNames(., stri_replace_first_regex(names(.) ,"s\\(([^\\)]+)\\)", "$1")))

partials <- as.data.frame(lapply(1:ncol(preds$fit), function(cl) {
  x = lfun(model_data_vir[[gterms$response]]) - rowSums(preds$fit[,-cl]) - intercept
  # (rowSums(preds$fit[,-cl]))
  # x = (lfun(model_data_vir[[gterms$response]]) - model_data_vir[[gterms$offset]])
  #  x = model_data_vir[[gterms$response]]
  x
}))
names(partials) <- names(preds$fit)

#smooth_titles = list("Extinct\nhome range", "Extinct minimum\nphylgenetic distance", "Disease\ncitations (log)", "Mammal\nsympatry") #"PVR body mass", bquote('Range (log' ~ km^{2} ~ ')'),
#names(smooth_titles) = names(smooth_data_vir)
smooth_plots_vir = map(names(smooth_data_vir), function(smooth_term_vir) {
  pl =  ggplot() +
    geom_hline(yintercept = 0, size=0.1, col="grey50") +
    geom_point(mapping = aes(x=model_data_vir[[smooth_term_vir]], y = (partials[[smooth_term_vir]])),
               shape=21, fill="red", col="red", alpha=0.25, size=2, stroke=0.1) +
    geom_ribbon(mapping = aes(x = smooth_ranges[[smooth_term_vir]],
                              ymin = (smooth_preds$fit[[smooth_term_vir]] - 2 * smooth_preds$se.fit[[smooth_term_vir]]),
                              ymax = (smooth_preds$fit[[smooth_term_vir]] + 2 * smooth_preds$se.fit[[smooth_term_vir]])),
                alpha = 0.5, fill=viridis(5)[4]) +
    geom_line(mapping = aes(x = smooth_ranges[[smooth_term_vir]], y = (smooth_preds$fit[[smooth_term_vir]])), size=0.3) +
    #xlab(smooth_titles[[smooth_term_vir]]) +
    scale_y_continuous(limits=c(-2.2,2.2), oob=scales::rescale_none) +
    theme_bw() + partials_theme
  
  #if (SHOW_DEV_EXPL) pl <- pl + annotate("label", x = max(model_data_vir[[smooth_term_vir]]), y = -2, label = paste0("DE = ", de_bgam$dev_explained[de_bgam$term == smooth_term_vir]), hjust = 1, size=1.5,  label.size=0, fill="#FFFFFF8C")
  #  geom_rug(mapping = aes(x =model_data_vir[[smooth_term]]), alpha=0.3) +
  
  return(pl)
  
})

xlabel <- expression()

#smooth_plots_vir[[1]] = smooth_plots_vir[[1]] + ylab() +
#  theme(axis.title.y=element_text(angle=90, lineheight=1.2, margin=margin(r=3))) 
#smooth_plots_vir[[1]] <- smooth_plots_vir[[1]] +
#  theme(axis.line.y.left = element_line(size = 0.5, 
#                                        color = 'grey50'),
#        axis.text.x = element_text(size = 10)) +
#  xlab(bquote(atop("Median home range of", "extinct hosts" ~ (log(Km^2)))))
  #xlab(bquote('Mean home range of extinct hosts '(Km^2)))
  #xlab(paste("Logged mean home\nrange of extinct\nhosts ", expression((Km)^2)))
  #xlab(bquote(atop("Mean home range of", "extinct hosts" ~ (Km^2))))
  #xlab("No. endangered hosts\n(Phylacine v1.2.1)")

#smooth_plots_vir[[1]]

smooth_plots_vir[[1]] = smooth_plots_vir[[1]] + theme(axis.title.y = element_blank()) +
  theme(axis.line.y.left = element_line(size = 0.5, color = 'grey50'),
        axis.text.x = element_text(size = 10))
smooth_plots_vir[[1]] <- smooth_plots_vir[[1]] +
  xlab("Median phylogenetic\ndistance to extinct\nhosts (log(Ma))")
  #xlab("No. extinct hosts")
smooth_plots_vir[[1]]

#smooth_plots_vir[[2]] = smooth_plots_vir[[2]] + theme(axis.title.y = element_blank()) +
#  theme(axis.line.y.left = element_line(size = 0.5, color = 'grey50'),
#        axis.text.x = element_text(size = 10))
#smooth_plots_vir[[2]] <- smooth_plots_vir[[2]] +
#  xlab("No. extinct hosts\nsame family")
  #xlab(bquote(atop("Mean home range of", "extinct hosts" ~ (Km^2))))
#smooth_plots_vir[[2]]

smooth_plots_vir[[2]] = smooth_plots_vir[[2]] + theme(axis.title.y = element_blank()) +
  theme(axis.line.y.left = element_line(size = 0.5, color = 'grey50'),
        axis.text.x = element_text(size = 10))
smooth_plots_vir[[2]] <- smooth_plots_vir[[2]] +
  xlab("No. extinct hosts\nsame family")
  #xlab("Min phylogenetic\ndistance to\nextinct hosts (Ma)")
smooth_plots_vir[[2]]

smooth_plots_vir[[3]] = smooth_plots_vir[[3]] + theme(axis.title.y = element_blank()) +
  theme(axis.line.y.left = element_line(size = 0.5, color = 'grey50'),
        axis.text.x = element_text(size = 10))
smooth_plots_vir[[3]] <- smooth_plots_vir[[3]] +
  xlab("Disease citations\n(log)")
  #xlab("No. extinct hosts\nsame order")
  #xlab(bquote('Range '(log(Km^2))))
smooth_plots_vir[[3]]

smooth_plots_vir[[4]] = smooth_plots_vir[[4]] + theme(axis.title.y = element_blank()) +
  theme(axis.line.y.left = element_line(size = 0.5, color = 'grey50'),
        axis.text.x = element_text(size = 10))
smooth_plots_vir[[4]] <- smooth_plots_vir[[4]] +
  xlab(bquote('Range '(log(Km^2))))
  #xlab("Disease citations\n(log)")
  #xlab("Plourde et al. (2017)\nPC 1")
smooth_plots_vir[[4]]

smooth_plots_vir[[5]] = smooth_plots_vir[[5]] + theme(axis.title.y = element_blank()) +
  theme(axis.line.y.left = element_line(size = 0.5, color = 'grey50'),
        axis.text.x = element_text(size = 10))
smooth_plots_vir[[5]] <- smooth_plots_vir[[5]] +
  #xlab("Plourde et al. (2017)\nPC 2")
  xlab("Mammal sympatry")
  #xlab("No disturbance\nsympatry")
  #xlab("PVR body mass")
smooth_plots_vir[[5]]


smooth_plots_vir[[6]] = smooth_plots_vir[[6]] + theme(axis.title.y = element_blank()) +
  theme(axis.line.y.left = element_line(size = 0.5, color = 'grey50'),
        axis.text.x = element_text(size = 10))
smooth_plots_vir[[6]] <- smooth_plots_vir[[6]] +
  xlab("No. extant hosts lost\nlocally (log)")
smooth_plots_vir[[6]]

pdplot_citation <- smooth_plots_vir[[3]]
pdplot_range <- smooth_plots_vir[[4]]
pdplot_sym <- smooth_plots_vir[[5]]
pdplot_ex <- smooth_plots_vir[[6]]

smooth_plots_vir[[2]] <- pdplot_citation
smooth_plots_vir[[3]] <- pdplot_range
smooth_plots_vir[[4]] <- pdplot_sym
smooth_plots_vir[[5]] <- pdplot_ex
smooth_plots_vir[[6]] <- blankPlot

bin_vir_data = binary_preds %>% map(function(x) {
  x = x[, stri_detect_regex(names(x), paste0("(", paste0(binary_vars, collapse="|"), ")"))]
  n = names(x)
  x = rowSums(x)
  data_frame(response=x, variable=n)
})


bin_vir_data$fit$se = bin_vir_data$se.fit$response
bin_vir_data = bin_vir_data$fit
bin_vir_data$response = bin_vir_data$response
bin_vir_data$labels = stri_replace_first_regex(bin_vir_data$variable, "hOrder", "")
#bin_vir_data$labels = stri_replace_first_regex(bin_vir_data$labels, "hHuntedIUCN", "Hunted")
bin_vir_data$signif = summary(bgam)$s.table[stri_detect_regex(rownames(summary(bgam)$s.table), paste0("(", paste0(binary_vars, collapse="|"), ")")), "p-value"] < 0.05
bin_vir_data = bin_vir_data %>%
  arrange(desc(signif), response) %>%
  mutate(no = 1:nrow(bin_vir_data))

bin_vir_partials = lapply(binary_terms, function(x) {
  vals = partials[as.logical(model_data_vir[[x]]), x-1]
  variable = names(model_data_vir)[x]
  variable_rep = rep(variable, length(vals))
  no_rep = rep(bin_vir_data$no[bin_vir_data$variable == variable], length(vals))
  data_frame(variable=variable_rep, partial=vals, no=no_rep)
}) %>% bind_rows

bin_vir_data = bin_vir_partials %>%
  group_by(variable) %>%
  summarize(minval = min(partial)) %>%
  inner_join(bin_vir_data, by="variable") %>%
  mutate(minval = pmin(minval, response - 2*se)) %>%
  left_join(de_bgam, by=c('variable' = 'term'))

bin_plot_vir = ggplot() +
  geom_hline(yintercept = 0, size=0.1, col="grey50") +
  geom_point(data=bin_vir_partials, mapping=aes(x=no, y=(partial)), position=position_jitter(width=0.25),
             shape=21, fill="red", col="red", alpha=0.25, size=2, stroke=0.1) +
  geom_rect(data = bin_vir_data, mapping=aes(xmin = no - 0.35, xmax  = no + 0.35, ymin=(response-2*se), ymax=(response+2*se), fill=signif), alpha = 0.5) +
  geom_segment(data = bin_vir_data, mapping=aes(x=no - 0.35, xend = no + 0.35, y=(response), yend=(response)), col="black", size=0.3) +
  scale_fill_manual(breaks = c(TRUE, FALSE), values=c(viridis(5)[4], viridis(5)[4])) +
  scale_x_continuous(breaks = bin_vir_data$no, labels = stri_trans_totitle(bin_vir_data$labels)) +
  scale_y_continuous(limits=c(-2,2.2), breaks=seq(-2,2, by=1), oob=scales::rescale_none, name="") +
  theme_bw() + partials_theme +
  theme(axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.x=element_text(colour="black", angle=90, hjust=1), legend.position="none",
        axis.title.y=element_blank(), axis.line.y = element_line(size = 0.5)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.line.y.left = element_line(size = 0.5, color = 'grey50'))

#if (SHOW_DEV_EXPL) bin_plot_vir = bin_plot_vir + geom_label(data = bin_vir_data, mapping=aes(x=no, y = response + 2*se + 0.4, label=dev_explained), color="black", family="Helvetica", size=1.5, label.size=0, fill="#FFFFFF8C") +


vir_plots <- plot_grid(textGrob("Strength of effect on\nviral richness", 
                                gp=gpar(fontface="bold", col="black", fontsize=10), rot=90),
                       plot_grid(plotlist = c(smooth_plots_vir[1:5]), 
                                 nrow=1, 
                                 align="h", 
                                 rel_widths = c(1,1,1,1,1),
                                 labels=c("(a)", "(b)", "(c)", "(d)", "(e)"), 
                                 label_size=10, 
                                 hjust=0),
                       bin_plot_vir, 
                       nrow= 1, 
                       rel_widths = c(0.5, 5, 1.3), 
                       labels=c("", "", "(f)"), 
                       label_size=10, 
                       hjust=0)


plot(vir_plots)
ggsave(paste0(base_dir, 'Figures/virus_plot.png'), width = 14, height=5)












#---
# For bacteria
model_dist <- 'tw'
bgam_ne <- readRDS(paste0(base_dir, 'intermediates/noExtinction_clover_bacteria_models_', model_dist, '.rds'))$model[[1]]
summary(bgam_ne)
anova.gam(bgam_ne)
# Check whether there are relationships between residulas using gam.check():
gam.check(bgam_ne, rep = 1000)
# Plot residuals
appraise(bgam_ne, method = 'simulate')


# Concurvity
concurvity(bgam_ne, full=T)
con <- concurvity(bgam_ne, full=F)
con_worst <-format(con$worst, scientific = F)
con_worst <- as.data.frame(con_worst)
#con_worst[con_worst < 0.8] <- 0


model_dist <- 'tw'
bgam = readRDS(paste0(base_dir, 'intermediates/clover_bacteria_models_', model_dist, '.rds'))$model[[1]]
summary(bgam)
# Generates table of the best model's summary stats:
bacteria_summary <- summary(bgam)
bac_table <- as.data.frame(bacteria_summary$s.table)
bac_table$p_value_nonsci <- format(bac_table$`p-value`, scientific = F)
bac_table$input_variable <- rownames(bac_table)
bac_table <- bac_table[, c(6,1,2,3,4,5)]
rownames(bac_table) <- NULL
bac_form <- bacteria_summary$formula
bac_form <- gsub("\"", "", bac_form)
bac_form <- sub("\n", "", bac_form)
bac_form <- sub("\t", "", bac_form)
names(bac_table)[names(bac_table) == 'input_variable'] <- paste0(as.character(bac_form)[2], " ", as.character(bac_form)[1], " ", as.character(bac_form)[3])
bac_table$edf <- formatC(bac_table$edf, format = 'f', digits = 2)
bac_table$F <- formatC(bac_table$F, format = 'f', digits = 2)
write.csv(bac_table, paste0(base_dir, "Figures/clover_bacteria_bgam_summary.csv"))

#anova.gam(bgam)
aic_scores <- as.data.frame(AIC(bgam_ne, bgam))
write.csv(aic_scores, paste0(base_dir, "Figures/clover_bacteria_bgam_aic.csv"))

#comp_gam <- anova.gam(bgam_ne, bgam, test = "F") #Chisq or F
#write.csv(comp_gam, paste0(base_dir, "Figures/clover_bacteria_gam_chisq.csv"))
# Concurvity
concurvity(bgam, full=T)
con <- concurvity(bgam, full=F)
con_worst <-format(con$worst, scientific = F)
#con_worst[con_worst < 0.8] <- 0


# Check whether there are relationships between residulas using gam.check():
gam.check(bgam, rep = 1000)
# Plot residuals
appraise(bgam, method = 'simulate')



de_bgam =  get_relative_contribs(bgam) %>%
  mutate(dev_explained = rel_deviance_explained * summary(bgam)$dev.expl) %>%
  mutate(dev_explained = paste0(stri_trim_both(formatC(dev_explained*100, format = "fg", digits=2)), "%"))
#de_bgam <- arrange(de_bgam, desc(as.numeric(sub("%", "", de_bgam$dev_explained))))
print(de_bgam)


# Generates table of the best 5 model's:
de_bgam1 <- de_bgam
bgam1_sum <- as.data.frame(summary(bgam)$s.table)
bgam1_sum_sig <- bgam1_sum[bgam1_sum$`p-value`<=0.05,]
bgam1_names <- gsub("s\\(", "", rownames(bgam1_sum_sig))
bgam1_names <- gsub("\\)", "", bgam1_names)
de_bgam1 =  get_relative_contribs(bgam) %>%
  mutate(dev_explained = rel_deviance_explained * summary(bgam)$dev.expl) %>%
  mutate(dev_explained = paste0(stri_trim_both(formatC(dev_explained*100, format = "fg", digits=2)), "%")) %>%
  filter(term %in% bgam1_names)

bgam2 = readRDS(paste0(base_dir, 'intermediates/clover_bacteria_models_', model_dist, '.rds'))$model[[2]]
bgam2_sum <- as.data.frame(summary(bgam2)$s.table)
bgam2_sum_sig <- bgam2_sum[bgam2_sum$`p-value`<=0.05,]
bgam2_names <- gsub("s\\(", "", rownames(bgam2_sum_sig))
bgam2_names <- gsub("\\)", "", bgam2_names)
de_bgam2 =  get_relative_contribs(bgam2) %>%
  mutate(dev_explained = rel_deviance_explained * summary(bgam2)$dev.expl) %>%
  mutate(dev_explained = paste0(stri_trim_both(formatC(dev_explained*100, format = "fg", digits=2)), "%")) %>%
  filter(term %in% bgam2_names)

bgam3 = readRDS(paste0(base_dir, 'intermediates/clover_bacteria_models_', model_dist, '.rds'))$model[[3]]
bgam3_sum <- as.data.frame(summary(bgam3)$s.table)
bgam3_sum_sig <- bgam3_sum[bgam3_sum$`p-value`<=0.05,]
bgam3_names <- gsub("s\\(", "", rownames(bgam3_sum_sig))
bgam3_names <- gsub("\\)", "", bgam3_names)
de_bgam3 =  get_relative_contribs(bgam3) %>%
  mutate(dev_explained = rel_deviance_explained * summary(bgam3)$dev.expl) %>%
  mutate(dev_explained = paste0(stri_trim_both(formatC(dev_explained*100, format = "fg", digits=2)), "%")) %>%
  filter(term %in% bgam3_names)

bgam4 = readRDS(paste0(base_dir, 'intermediates/clover_bacteria_models_', model_dist, '.rds'))$model[[4]]
bgam4_sum <- as.data.frame(summary(bgam4)$s.table)
bgam4_sum_sig <- bgam4_sum[bgam4_sum$`p-value`<=0.05,]
bgam4_names <- gsub("s\\(", "", rownames(bgam4_sum_sig))
bgam4_names <- gsub("\\)", "", bgam4_names)
de_bgam4 =  get_relative_contribs(bgam4) %>%
  mutate(dev_explained = rel_deviance_explained * summary(bgam4)$dev.expl) %>%
  mutate(dev_explained = paste0(stri_trim_both(formatC(dev_explained*100, format = "fg", digits=2)), "%")) %>%
  filter(term %in% bgam4_names)


bgam5 = readRDS(paste0(base_dir, 'intermediates/clover_bacteria_models_', model_dist, '.rds'))$model[[5]]
bgam5_sum <- as.data.frame(summary(bgam5)$s.table)
bgam5_sum_sig <- bgam5_sum[bgam5_sum$`p-value`<=0.05,]
bgam5_names <- gsub("s\\(", "", rownames(bgam5_sum_sig))
bgam5_names <- gsub("\\)", "", bgam5_names)
de_bgam5 =  get_relative_contribs(bgam5) %>%
  mutate(dev_explained = rel_deviance_explained * summary(bgam5)$dev.expl) %>%
  mutate(dev_explained = paste0(stri_trim_both(formatC(dev_explained*100, format = "fg", digits=2)), "%")) %>%
  filter(term %in% bgam5_names)

de_bgam1 <- as.data.frame(de_bgam1)
de_bgam1 <- cbind(de_bgam1, rep("1", nrow(de_bgam1)))
colnames(de_bgam1)[4] <- "Model_rank"
de_bgam2 <- as.data.frame(de_bgam2)
de_bgam2 <- cbind(de_bgam2, rep("2", nrow(de_bgam2)))
colnames(de_bgam2)[4] <- "Model_rank"
de_bgam3 <- as.data.frame(de_bgam3)
de_bgam3 <- cbind(de_bgam3, rep("3", nrow(de_bgam3)))
colnames(de_bgam3)[4] <- "Model_rank"
de_bgam4 <- as.data.frame(de_bgam4)
de_bgam4 <- cbind(de_bgam4, rep("4", nrow(de_bgam4)))
colnames(de_bgam4)[4] <- "Model_rank"
de_bgam5 <- as.data.frame(de_bgam5)
de_bgam5 <- cbind(de_bgam5, rep("5", nrow(de_bgam5)))
colnames(de_bgam5)[4] <- "Model_rank"
top5_model <- rbind(de_bgam1, de_bgam2, de_bgam3, de_bgam4, de_bgam5)
top5_model_results <- top5_model %>%
  group_by(term) %>%
  summarize(n = n())
write.csv(top5_model_results, paste0(base_dir, "Figures/bacteria_top5gam_variable_breakdown.csv"))
rel_deviance_perc <- as.numeric(de_bgam$rel_deviance_explained*100)
de_bgam$rel_deviance_explained_perc <- paste0(formatC(rel_deviance_perc, format = 'f', digits = 2), "%")
de_bgam$deviance_explained_numeric <- as.numeric(sub("%", "", de_bgam$dev_explained))
write.csv(de_bgam, paste0(base_dir, "Figures/clover_bacteria_gam_deviance.csv"))

ex_v_de <- sum(de_bgam$deviance_explained_numeric[de_bgam$term == "extinct.n.family"],
               de_bgam$deviance_explained_numeric[de_bgam$term == "extinct.host.no"])

bio_v_de <- sum(de_bgam$deviance_explained_numeric[de_bgam$term == "extinct.n.family"],
                de_bgam$deviance_explained_numeric[de_bgam$term == "extinct.host.no"],
                de_bgam$deviance_explained_numeric[de_bgam$term == "LnAreaHost"],
                de_bgam$deviance_explained_numeric[de_bgam$term == "Plourde_PC1"],
                de_bgam$deviance_explained_numeric[de_bgam$term == "S20"],
                de_bgam$deviance_explained_numeric[de_bgam$term == "hOrderCETARTIODACTYLA"],
                de_bgam$deviance_explained_numeric[de_bgam$term == "hOrderLAGOMORPHA"],
                de_bgam$deviance_explained_numeric[de_bgam$term == "hOrderRODENTIA"],
                de_bgam$deviance_explained_numeric[de_bgam$term == "S20"])

paste0("Relative deviance (of CLOVER richness per host) across all extinction variables: ", ex_v_de, " %")
paste0("Relative deviance (of CLOVER richness per host) for biologically relevant variables: ", bio_v_de, " %")


# All bacteria Plot
binary_vars = c("hOrder")

preds <- predict(bgam, type="iterms", se.fit=TRUE)
intercept <- attributes(preds)$constant

preds = preds %>% map(as.data.frame) %>%
  map(~setNames(., stri_replace_first_regex(names(.) ,"s\\(([^\\)]+)\\)", "$1")))
model_data_bac = bgam$model
gterms = attributes(terms(bgam))
ilfun <- bgam$family$linkinv
lfun <- bgam$family$linkfun

binary_terms = which(stri_detect_regex(names(model_data_bac), paste0("(", paste0(binary_vars, collapse="|"), ")")))

smooth_data_bac = model_data_bac[, -c(binary_terms, gterms$response, gterms$offset)]
smooth_ranges = data.frame(map(smooth_data_bac, ~seq(min(.), max(.), length.out = 100)))

binary_data = model_data_bac[, binary_terms]
binary_ranges = setNames(as.data.frame(diag(length(binary_terms))), names(model_data_bac)[binary_terms])

offset_name = stri_replace_first_regex(names(model_data_bac)[gterms$offset], "offset\\(([^\\)]+)\\)", "$1")

smooth_ranges = cbind(smooth_ranges,
                      setNames(as.data.frame(lapply(c(names(binary_data), offset_name), function(x) rep(0, nrow(smooth_ranges)))), c(names(binary_data), offset_name)))
binary_ranges = cbind(binary_ranges,
                      setNames(as.data.frame(lapply(c(names(smooth_data_bac), offset_name), function(x) rep(0, nrow(binary_ranges)))), c(names(smooth_data_bac), offset_name)))

smooth_preds <- predict(bgam, newdata=smooth_ranges, type="iterms", se.fit=TRUE)  %>% map(as.data.frame) %>%
  map(~setNames(., stri_replace_first_regex(names(.) ,"s\\(([^\\)]+)\\)", "$1")))
binary_preds <- predict(bgam, newdata=binary_ranges, type="iterms", se.fit=TRUE) %>% map(as.data.frame) %>%
  map(~setNames(., stri_replace_first_regex(names(.) ,"s\\(([^\\)]+)\\)", "$1")))

partials <- as.data.frame(lapply(1:ncol(preds$fit), function(cl) {
  x = lfun(model_data_bac[[gterms$response]]) - rowSums(preds$fit[,-cl]) - intercept
  # (rowSums(preds$fit[,-cl]))
  # x = (lfun(model_data_bac[[gterms$response]]) - model_data_bac[[gterms$offset]])
  #  x = model_data_bac[[gterms$response]]
  x
}))
names(partials) <- names(preds$fit)

#smooth_titles = list("Extinct\nfecal diffusivity", "Extinct minimum\nphylgenetic distance", "Disease\ncitations (log)", "Mammal\nsympatry") #bquote('Range (log' ~ km^{2} ~ ')'),
#names(smooth_titles) = names(smooth_data_bac)
smooth_plots_bac = map(names(smooth_data_bac), function(smooth_term_bac) {
  pl =  ggplot() +
    geom_hline(yintercept = 0, size=0.1, col="grey50") +
    geom_point(mapping = aes(x=model_data_bac[[smooth_term_bac]], y = (partials[[smooth_term_bac]])),
               shape=21, fill="deepskyblue3", col="deepskyblue3", alpha=0.25, size=2, stroke=0.1) +
    geom_ribbon(mapping = aes(x = smooth_ranges[[smooth_term_bac]],
                              ymin = (smooth_preds$fit[[smooth_term_bac]] - 2 * smooth_preds$se.fit[[smooth_term_bac]]),
                              ymax = (smooth_preds$fit[[smooth_term_bac]] + 2 * smooth_preds$se.fit[[smooth_term_bac]])),
                alpha = 0.5, fill=viridis(5)[4]) +
    geom_line(mapping = aes(x = smooth_ranges[[smooth_term_bac]], y = (smooth_preds$fit[[smooth_term_bac]])), size=0.3) +
    #xlab(smooth_titles[[smooth_term_bac]]) +
    scale_y_continuous(limits=c(-3,2.2), oob=scales::rescale_none) +
    theme_bw() + partials_theme
  
  #if (SHOW_DEV_EXPL) pl <- pl + annotate("label", x = max(model_data_bac[[smooth_term_bac]]), y = -2, label = paste0("DE = ", de_bgam$dev_explained[de_bgam$term == smooth_term_bac]), hjust = 1, size=1.5,  label.size=0, fill="#FFFFFF8C")
  #  geom_rug(mapping = aes(x =model_data_bac[[smooth_term]]), alpha=0.3) +
  
  return(pl)
  
})

smooth_plots_bac[[1]] = smooth_plots_bac[[1]] + ylab("") +
  theme(axis.title.y=element_text(angle=90, lineheight=1.2, margin=margin(r=3))) 
smooth_plots_bac[[1]] <- smooth_plots_bac[[1]] +
  theme(axis.line.y.left = element_line(size = 0.5, 
                                        color = 'grey50'),
        axis.text.x = element_text(size = 10)) +
  xlab("No. extinct hosts")
  #xlab("Median phylogenetic\ndistance to extinct\nhosts (log(Ma))")
  #xlab("No. extinct hosts\nsame family")
  #xlab(bquote(atop("Median home range of", "extinct hosts" ~ (log(Km^2)))))
smooth_plots_bac[[1]]

smooth_plots_bac[[2]] = smooth_plots_bac[[2]] + theme(axis.title.y = element_blank()) +
  theme(axis.line.y.left = element_line(size = 0.5, color = 'grey50'),
        axis.text.x = element_text(size = 10))
smooth_plots_bac[[2]] <- smooth_plots_bac[[2]] +
  #xlab("Mean home range of\nextinct hosts (Km)")
  xlab("Disease citations\n(log)")
  #xlab("Median phylogenetic\ndistance to extinct\nhosts (log(Ma))")
smooth_plots_bac[[2]]

smooth_plots_bac[[3]] = smooth_plots_bac[[3]] + theme(axis.title.y = element_blank()) +
  theme(axis.line.y.left = element_line(size = 0.5, color = 'grey50'),
        axis.text.x = element_text(size = 10))
smooth_plots_bac[[3]] <- smooth_plots_bac[[3]] +
  #xlab("Mean phylogenetic\ndistance to\nextinct hosts (Ma)")
  #xlab("No. extinct hosts\nsame family")
  #xlab("PVR body mass")
  xlab("Plourde et al. (2017)\nPC 1")
  #xlab(bquote('Range '(log(Km^2))))
smooth_plots_bac[[3]]

smooth_plots_bac[[4]] = smooth_plots_bac[[4]] + theme(axis.title.y = element_blank()) +
  theme(axis.line.y.left = element_line(size = 0.5, color = 'grey50'),
        axis.text.x = element_text(size = 10))
smooth_plots_bac[[4]] <- smooth_plots_bac[[4]] +
  #xlab("No. extinct hosts\nsame family")
  #xlab("Disease citations\n(log)")
  xlab("Plourde et al. (2017)\nPC 2")
  
smooth_plots_bac[[4]]

smooth_plots_bac[[5]] = smooth_plots_bac[[5]] + theme(axis.title.y = element_blank()) +
  theme(axis.line.y.left = element_line(size = 0.5, color = 'grey50'),
        axis.text.x = element_text(size = 10))
smooth_plots_bac[[5]] <- smooth_plots_bac[[5]] +
  #xlab("No. extinct hosts\nsame family")
  #xlab("Disease citations\n(log)")
  #xlab("Plourde et al. (2017)\nPC 2")
  #xlab("Threatened sympatric\nspecies (Phylacine)")
  xlab("Mammal sympatry")
smooth_plots_bac[[5]]

smooth_plots_bac[[6]] = smooth_plots_bac[[6]] + theme(axis.title.y = element_blank()) +
  theme(axis.line.y.left = element_line(size = 0.5, color = 'grey50'),
        axis.text.x = element_text(size = 10))
smooth_plots_bac[[6]] <- smooth_plots_bac[[6]] +
  #xlab(bquote('Range '(log(Km^2))))
  xlab("No. extant hosts lost\nlocally (log)")
smooth_plots_bac[[6]]


pdplot_sym <- smooth_plots_bac[[5]]
pdplot_ex <- smooth_plots_bac[[6]]
smooth_plots_bac[[4]] <- pdplot_sym
smooth_plots_bac[[5]] <- pdplot_ex
smooth_plots_bac[[6]] <- blankPlot


bin_bac_data = binary_preds %>% map(function(x) {
  x = x[, stri_detect_regex(names(x), paste0("(", paste0(binary_vars, collapse="|"), ")"))]
  n = names(x)
  x = rowSums(x)
  data_frame(response=x, variable=n)
})


bin_bac_data$fit$se = bin_bac_data$se.fit$response
bin_bac_data = bin_bac_data$fit
bin_bac_data$response = bin_bac_data$response
bin_bac_data$labels = stri_replace_first_regex(bin_bac_data$variable, "hOrder", "")
#bin_bac_data$labels = stri_replace_first_regex(bin_bac_data$labels, "hHuntedIUCN", "Hunted")
bin_bac_data$signif = summary(bgam)$s.table[stri_detect_regex(rownames(summary(bgam)$s.table), paste0("(", paste0(binary_vars, collapse="|"), ")")), "p-value"] < 0.05
bin_bac_data = bin_bac_data %>%
  arrange(desc(signif), response) %>%
  mutate(no = 1:nrow(bin_bac_data))

bin_bac_partials = lapply(binary_terms, function(x) {
  vals = partials[as.logical(model_data_bac[[x]]), x-1]
  variable = names(model_data_bac)[x]
  variable_rep = rep(variable, length(vals))
  no_rep = rep(bin_bac_data$no[bin_bac_data$variable == variable], length(vals))
  data_frame(variable=variable_rep, partial=vals, no=no_rep)
}) %>% bind_rows

bin_bac_data = bin_bac_partials %>%
  group_by(variable) %>%
  summarize(minval = min(partial)) %>%
  inner_join(bin_bac_data, by="variable") %>%
  mutate(minval = pmin(minval, response - 2*se)) %>%
  left_join(de_bgam, by=c('variable' = 'term'))

bin_plot_bac = ggplot() +
  geom_hline(yintercept = 0, size=0.1, col="grey50") +
  geom_point(data=bin_bac_partials, mapping=aes(x=no, y=(partial)), position=position_jitter(width=0.25),
             shape=21, fill="deepskyblue3", col="deepskyblue3", alpha=0.25, size=2, stroke=0.1) +
  geom_rect(data = bin_bac_data, mapping=aes(xmin = no - 0.35, xmax  = no + 0.35, ymin=(response-2*se), ymax=(response+2*se), fill=signif), alpha = 0.5) +
  geom_segment(data = bin_bac_data, mapping=aes(x=no - 0.35, xend = no + 0.35, y=(response), yend=(response)), col="black", size=0.3) +
  scale_fill_manual(breaks = c(TRUE, FALSE), values=c(viridis(5)[4], viridis(5)[4])) +
  scale_x_continuous(breaks = bin_bac_data$no, labels = stri_trans_totitle(bin_bac_data$labels)) +
  scale_y_continuous(limits=c(-3,2.2), breaks=seq(-2,2, by=1), oob=scales::rescale_none, name="") +
  theme_bw() + partials_theme +
  theme(axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.x=element_text(colour="black", angle=90, hjust=1), legend.position="none",
        axis.title.y=element_blank()) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.line.y.left = element_line(size = 0.5, color = 'grey50'))

bac_plots <- plot_grid(textGrob("Strength of effect on\nbacterial richness", 
                                gp=gpar(fontface="bold", col="black", fontsize=10), rot=90),
                       plot_grid(plotlist = c(smooth_plots_bac[1:5]), 
                                 nrow=1, 
                                 align="h", 
                                 rel_widths = c(1,1,1,1,1),
                                 labels=c("(g)", "(h)", "(i)", "(j)", "(k)"), 
                                 label_size=10, 
                                 hjust=0),
                       bin_plot_bac, 
                       nrow= 1, 
                       rel_widths = c(0.5, 5, 1.3), 
                       labels=c("", "", "(l)"), 
                       label_size=10, 
                       hjust=0)



plot(bac_plots)
ggsave(paste0(base_dir, 'Figures/bacteria_plot.png'), width = 14, height=5)









model_dist <- 'tw'
bgam_ne <- readRDS(paste0(base_dir, 'intermediates/noExtinction_clover_helminth_models_', model_dist, '.rds'))$model[[1]]
summary(bgam_ne)
anova.gam(bgam_ne)
# Check whether there are relationships between residulas using gam.check():
gam.check(bgam_ne, rep = 1000)
# Plot residuals
appraise(bgam_ne, method = 'simulate')
# Concurvity:
concurvity(bgam_ne, full=T)
con <- concurvity(bgam_ne, full=F)
con_worst <-format(con$worst, scientific = F)



model_dist <- "tw"
bgam = readRDS(paste0(base_dir, 'intermediates/clover_helminth_models_', model_dist, '.rds'))$model[[1]]
summary(bgam)
# Generates table of the best model's summary stats:
helminth_summary <- summary(bgam)
helminth_table <- as.data.frame(helminth_summary$s.table)
helminth_table$p_value_nonsci <- format(helminth_table$`p-value`, scientific = F)
helminth_table$input_variable <- rownames(helminth_table)
helminth_table <- helminth_table[, c(6,1,2,3,4,5)]
rownames(helminth_table) <- NULL
h_form <- helminth_summary$formula
h_form <- gsub("\"", "", h_form)
h_form <- sub("\n", "", h_form)
h_form <- sub("\t", "", h_form)
names(helminth_table)[names(helminth_table) == 'input_variable'] <- paste0(as.character(h_form)[2], " ", as.character(h_form)[1], " ", as.character(h_form)[3])
helminth_table$edf <- formatC(helminth_table$edf, format = 'f', digits = 2)
helminth_table$F <- formatC(helminth_table$F, format = 'f', digits = 2)
write.csv(helminth_table, paste0(base_dir, "Figures/clover_helminth_bgam_summary.csv"))

#anova.gam(bgam)
aic_scores <- as.data.frame(AIC(bgam_ne, bgam))
write.csv(aic_scores, paste0(base_dir, "Figures/clover_helminth_bgam_aic.csv"))

comp_gam <- anova.gam(bgam_ne, bgam, test = "Chisq")
write.csv(comp_gam, paste0(base_dir, "Figures/clover_helminth_gam_chisq.csv"))

# Check whether there are relationships between residulas using gam.check():
gam.check(bgam, rep = 1000)
# Plot residuals
appraise(bgam, method = 'simulate')
# Concurvity:
concurvity(bgam, full=T)
con <- concurvity(bgam, full=F)
con_worst <-format(con$worst, scientific = F)
#con_worst[con_worst < 0.8] <- 0


de_bgam =  get_relative_contribs(bgam) %>%
  mutate(dev_explained = rel_deviance_explained * summary(bgam)$dev.expl) %>%
  mutate(dev_explained = paste0(stri_trim_both(formatC(dev_explained*100, format = "fg", digits=2)), "%"))
#de_bgam <- arrange(de_bgam, desc(as.numeric(sub("%", "", de_bgam$dev_explained))))
print(de_bgam)

# Generates table of the best 5 model's:
de_bgam1 <- de_bgam
bgam1_sum <- as.data.frame(summary(bgam)$s.table)
bgam1_sum_sig <- bgam1_sum[bgam1_sum$`p-value`<=0.05,]
bgam1_names <- gsub("s\\(", "", rownames(bgam1_sum_sig))
bgam1_names <- gsub("\\)", "", bgam1_names)
de_bgam1 =  get_relative_contribs(bgam) %>%
  mutate(dev_explained = rel_deviance_explained * summary(bgam)$dev.expl) %>%
  mutate(dev_explained = paste0(stri_trim_both(formatC(dev_explained*100, format = "fg", digits=2)), "%")) %>%
  filter(term %in% bgam1_names)

bgam2 = readRDS(paste0(base_dir, 'intermediates/clover_helminth_models_', model_dist, '.rds'))$model[[2]]
bgam2_sum <- as.data.frame(summary(bgam2)$s.table)
bgam2_sum_sig <- bgam2_sum[bgam2_sum$`p-value`<=0.05,]
bgam2_names <- gsub("s\\(", "", rownames(bgam2_sum_sig))
bgam2_names <- gsub("\\)", "", bgam2_names)
de_bgam2 =  get_relative_contribs(bgam2) %>%
  mutate(dev_explained = rel_deviance_explained * summary(bgam2)$dev.expl) %>%
  mutate(dev_explained = paste0(stri_trim_both(formatC(dev_explained*100, format = "fg", digits=2)), "%")) %>%
  filter(term %in% bgam2_names)

bgam3 = readRDS(paste0(base_dir, 'intermediates/clover_helminth_models_', model_dist, '.rds'))$model[[3]]
bgam3_sum <- as.data.frame(summary(bgam3)$s.table)
bgam3_sum_sig <- bgam3_sum[bgam3_sum$`p-value`<=0.05,]
bgam3_names <- gsub("s\\(", "", rownames(bgam3_sum_sig))
bgam3_names <- gsub("\\)", "", bgam3_names)
de_bgam3 =  get_relative_contribs(bgam3) %>%
  mutate(dev_explained = rel_deviance_explained * summary(bgam3)$dev.expl) %>%
  mutate(dev_explained = paste0(stri_trim_both(formatC(dev_explained*100, format = "fg", digits=2)), "%")) %>%
  filter(term %in% bgam3_names)

bgam4 = readRDS(paste0(base_dir, 'intermediates/clover_helminth_models_', model_dist, '.rds'))$model[[4]]
bgam4_sum <- as.data.frame(summary(bgam4)$s.table)
bgam4_sum_sig <- bgam4_sum[bgam4_sum$`p-value`<=0.05,]
bgam4_names <- gsub("s\\(", "", rownames(bgam4_sum_sig))
bgam4_names <- gsub("\\)", "", bgam4_names)
de_bgam4 =  get_relative_contribs(bgam4) %>%
  mutate(dev_explained = rel_deviance_explained * summary(bgam4)$dev.expl) %>%
  mutate(dev_explained = paste0(stri_trim_both(formatC(dev_explained*100, format = "fg", digits=2)), "%")) %>%
  filter(term %in% bgam4_names)


bgam5 = readRDS(paste0(base_dir, 'intermediates/clover_helminth_models_', model_dist, '.rds'))$model[[5]]
bgam5_sum <- as.data.frame(summary(bgam5)$s.table)
bgam5_sum_sig <- bgam5_sum[bgam5_sum$`p-value`<=0.05,]
bgam5_names <- gsub("s\\(", "", rownames(bgam5_sum_sig))
bgam5_names <- gsub("\\)", "", bgam5_names)
de_bgam5 =  get_relative_contribs(bgam5) %>%
  mutate(dev_explained = rel_deviance_explained * summary(bgam5)$dev.expl) %>%
  mutate(dev_explained = paste0(stri_trim_both(formatC(dev_explained*100, format = "fg", digits=2)), "%")) %>%
  filter(term %in% bgam5_names)

de_bgam1 <- as.data.frame(de_bgam1)
de_bgam1 <- cbind(de_bgam1, rep("1", nrow(de_bgam1)))
colnames(de_bgam1)[4] <- "Model_rank"
de_bgam2 <- as.data.frame(de_bgam2)
de_bgam2 <- cbind(de_bgam2, rep("2", nrow(de_bgam2)))
colnames(de_bgam2)[4] <- "Model_rank"
de_bgam3 <- as.data.frame(de_bgam3)
de_bgam3 <- cbind(de_bgam3, rep("3", nrow(de_bgam3)))
colnames(de_bgam3)[4] <- "Model_rank"
de_bgam4 <- as.data.frame(de_bgam4)
de_bgam4 <- cbind(de_bgam4, rep("4", nrow(de_bgam4)))
colnames(de_bgam4)[4] <- "Model_rank"
de_bgam5 <- as.data.frame(de_bgam5)
de_bgam5 <- cbind(de_bgam5, rep("5", nrow(de_bgam5)))
colnames(de_bgam5)[4] <- "Model_rank"
top5_model <- rbind(de_bgam1, de_bgam2, de_bgam3, de_bgam4, de_bgam5)
top5_model_results <- top5_model %>%
  group_by(term) %>%
  summarize(n = n())
write.csv(top5_model_results, paste0(base_dir, "Figures/helminth_top5gam_variable_breakdown.csv"))
rel_deviance_perc <- as.numeric(de_bgam$rel_deviance_explained*100)
de_bgam$rel_deviance_explained_perc <- paste0(formatC(rel_deviance_perc, format = 'f', digits = 2), "%")
de_bgam$deviance_explained_numeric <- as.numeric(sub("%", "", de_bgam$dev_explained))
write.csv(de_bgam, paste0(base_dir, "Figures/clover_helminth_gam_deviance.csv"))

ex_v_de <- sum(de_bgam$deviance_explained_numeric[de_bgam$term == "extinct.mean.phylo.dist"])

bio_v_de <- sum(de_bgam$deviance_explained_numeric[de_bgam$term == "extinct.mean.phylo.dist"],
                de_bgam$deviance_explained_numeric[de_bgam$term == "hMassGramsPVR"],
                de_bgam$deviance_explained_numeric[de_bgam$term == "LnAreaHost"])

paste0("Relative deviance (of CLOVER richness per host) across all extinction variables: ", ex_v_de, " %")
paste0("Relative deviance (of CLOVER richness per host) for biologically relevant variables: ", bio_v_de, " %")

# All Viruses Plot
binary_vars = c("hOrder")

preds <- predict(bgam, type="iterms", se.fit=TRUE)
intercept <- attributes(preds)$constant

preds = preds %>% map(as.data.frame) %>%
  map(~setNames(., stri_replace_first_regex(names(.) ,"s\\(([^\\)]+)\\)", "$1")))
model_data_hel = bgam$model
gterms = attributes(terms(bgam))
ilfun <- bgam$family$linkinv
lfun <- bgam$family$linkfun

binary_terms = which(stri_detect_regex(names(model_data_hel), paste0("(", paste0(binary_vars, collapse="|"), ")")))

smooth_data_hel = model_data_hel[, -c(binary_terms, gterms$response, gterms$offset)]
smooth_ranges = data.frame(map(smooth_data_hel, ~seq(min(.), max(.), length.out = 100)))

binary_data = as.data.frame(model_data_hel[, binary_terms])
colnames(binary_data) <- names(model_data_hel)[binary_terms]
binary_ranges = setNames(as.data.frame(diag(length(binary_terms))), names(model_data_hel)[binary_terms])

offset_name = stri_replace_first_regex(names(model_data_hel)[gterms$offset], "offset\\(([^\\)]+)\\)", "$1")

smooth_ranges = cbind(smooth_ranges,
                      setNames(as.data.frame(lapply(c(names(binary_data), offset_name), function(x) rep(0, nrow(smooth_ranges)))), c(names(binary_data), offset_name)))
binary_ranges = cbind(binary_ranges,
                      setNames(as.data.frame(lapply(c(names(smooth_data_hel), offset_name), function(x) rep(0, nrow(binary_ranges)))), c(names(smooth_data_hel), offset_name)))

smooth_preds <- predict(bgam, newdata=smooth_ranges, type="iterms", se.fit=TRUE)  %>% map(as.data.frame) %>%
  map(~setNames(., stri_replace_first_regex(names(.) ,"s\\(([^\\)]+)\\)", "$1")))
binary_preds <- predict(bgam, newdata=binary_ranges, type="iterms", se.fit=TRUE) %>% map(as.data.frame) %>%
  map(~setNames(., stri_replace_first_regex(names(.) ,"s\\(([^\\)]+)\\)", "$1")))

partials <- as.data.frame(lapply(1:ncol(preds$fit), function(cl) {
  x = lfun(model_data_hel[[gterms$response]]) - rowSums(preds$fit[,-cl]) - intercept
  # (rowSums(preds$fit[,-cl]))
  # x = (lfun(model_data_hel[[gterms$response]]) - model_data_hel[[gterms$offset]])
  #  x = model_data_hel[[gterms$response]]
  x
}))
names(partials) <- names(preds$fit)

#smooth_titles = list("Extinct\nhome range", "Extinct minimum\nphylgenetic distance", "Disease\ncitations (log)", "Mammal\nsympatry") #"PVR body mass", bquote('Range (log' ~ km^{2} ~ ')'),
#names(smooth_titles) = names(smooth_data_hel)
smooth_plots_helminth = map(names(smooth_data_hel), function(smooth_term_hel) {
  pl =  ggplot() +
    geom_hline(yintercept = 0, size=0.1, col="grey50") +
    geom_point(mapping = aes(x=model_data_hel[[smooth_term_hel]], y = (partials[[smooth_term_hel]])),
               shape=21, fill="purple", col="purple", alpha=0.25, size=2, stroke=0.1) +
    geom_ribbon(mapping = aes(x = smooth_ranges[[smooth_term_hel]],
                              ymin = (smooth_preds$fit[[smooth_term_hel]] - 2 * smooth_preds$se.fit[[smooth_term_hel]]),
                              ymax = (smooth_preds$fit[[smooth_term_hel]] + 2 * smooth_preds$se.fit[[smooth_term_hel]])),
                alpha = 0.5, fill=viridis(5)[4]) +
    geom_line(mapping = aes(x = smooth_ranges[[smooth_term_hel]], y = (smooth_preds$fit[[smooth_term_hel]])), size=0.3) +
    #xlab(smooth_titles[[smooth_term_hel]]) +
    scale_y_continuous(limits=c(-2.2,2.2), oob=scales::rescale_none) +
    theme_bw() + partials_theme
  
  #if (SHOW_DEV_EXPL) pl <- pl + annotate("label", x = max(model_data_hel[[smooth_term_hel]]), y = -2, label = paste0("DE = ", de_bgam$dev_explained[de_bgam$term == smooth_term_hel]), hjust = 1, size=1.5,  label.size=0, fill="#FFFFFF8C")
  #  geom_rug(mapping = aes(x =model_data_hel[[smooth_term]]), alpha=0.3) +
  
  return(pl)
  
})


smooth_plots_helminth[[1]] = smooth_plots_helminth[[1]] + ylab() +
  theme(axis.title.y=element_text(angle=90, lineheight=1.2, margin=margin(r=3))) 
smooth_plots_helminth[[1]] <- smooth_plots_helminth[[1]] +
  theme(axis.line.y.left = element_line(size = 0.5, 
                                        color = 'grey50'),
        axis.text.x = element_text(size = 10)) +
  #xlab(bquote('Mean home range of extinct hosts '(Km^2)))
  #xlab(paste("Mean home range of\nextinct hosts ", expression((Km)^2)))
  #bquote(atop("Mean home range of", "extinct hosts" ~ (Km^2)))
  #xlab("No extinct hosts")
  xlab("Mean phylogenetic\ndistance to extinct\nhosts (log(Ma))")
smooth_plots_helminth[[1]]

smooth_plots_helminth[[2]] = smooth_plots_helminth[[2]] + theme(axis.title.y = element_blank()) +
  theme(axis.line.y.left = element_line(size = 0.5, color = 'grey50'),
        axis.text.x = element_text(size = 10))
smooth_plots_helminth[[2]] <- smooth_plots_helminth[[2]] +
  #xlab("Mean phylogenetic\ndistance to\nextinct hosts (Ma)")
  xlab("Disease citations\n(log)")
smooth_plots_helminth[[2]]

smooth_plots_helminth[[3]] = smooth_plots_helminth[[3]] + theme(axis.title.y = element_blank()) +
  theme(axis.line.y.left = element_line(size = 0.5, color = 'grey50'),
        axis.text.x = element_text(size = 10))
smooth_plots_helminth[[3]] <- smooth_plots_helminth[[3]] +
  #xlab(bquote('Range '(log(Km^2))))
  xlab("PVR body mass")
#xlab("Median day range of\nextinct hosts (Km)")
smooth_plots_helminth[[3]]

smooth_plots_helminth[[4]] = smooth_plots_helminth[[4]] + theme(axis.title.y = element_blank()) +
  theme(axis.line.y.left = element_line(size = 0.5, color = 'grey50'),
        axis.text.x = element_text(size = 10))
smooth_plots_helminth[[4]] <- smooth_plots_helminth[[4]] +
  #xlab("Plourde et al. (2017)\nPC 1")
  #xlab("No. extinct hosts\nsame order")
  xlab(bquote('Range '(log(Km^2))))
smooth_plots_helminth[[4]]


mass_pdplot <- smooth_plots_helminth[[3]]
range_pdplot <- smooth_plots_helminth[[4]]
smooth_plots_helminth[[3]] = range_pdplot
smooth_plots_helminth[[4]] = mass_pdplot 
smooth_plots_helminth[[5]] = blankPlot
smooth_plots_helminth[[6]] = blankPlot


bin_hel_data = binary_preds %>% map(function(x) {
  x = x[, stri_detect_regex(names(x), paste0("(", paste0(binary_vars, collapse="|"), ")"))]
  n = names(x)
  x = rowSums(x)
  data_frame(response=x, variable=n)
})


bin_hel_data$fit$se = bin_hel_data$se.fit$response
bin_hel_data = bin_hel_data$fit
bin_hel_data$response = bin_hel_data$response
bin_hel_data$labels = stri_replace_first_regex(bin_hel_data$variable, "hOrder", "")
#bin_hel_data$labels = stri_replace_first_regex(bin_hel_data$labels, "hHuntedIUCN", "Hunted")
bin_hel_data$signif = summary(bgam)$s.table[stri_detect_regex(rownames(summary(bgam)$s.table), paste0("(", paste0(binary_vars, collapse="|"), ")")), "p-value"] < 0.05
bin_hel_data = bin_hel_data %>%
  arrange(desc(signif), response) %>%
  mutate(no = 1:nrow(bin_hel_data))

bin_hel_partials = lapply(binary_terms, function(x) {
  vals = partials[as.logical(model_data_hel[[x]]), x-1]
  variable = names(model_data_hel)[x]
  variable_rep = rep(variable, length(vals))
  no_rep = rep(bin_hel_data$no[bin_hel_data$variable == variable], length(vals))
  data_frame(variable=variable_rep, partial=vals, no=no_rep)
}) %>% bind_rows

bin_hel_data = bin_hel_partials %>%
  group_by(variable) %>%
  summarize(minval = min(partial)) %>%
  inner_join(bin_hel_data, by="variable") %>%
  mutate(minval = pmin(minval, response - 2*se)) %>%
  left_join(de_bgam, by=c('variable' = 'term'))

bin_plot_hel = ggplot() +
  geom_hline(yintercept = 0, size=0.1, col="grey50") +
  geom_point(data=bin_hel_partials, mapping=aes(x=no, y=(partial)), position=position_jitter(width=0.25),
             shape=21, fill="purple", col="purple", alpha=0.25, size=2, stroke=0.1) +
  geom_rect(data = bin_hel_data, mapping=aes(xmin = no - 0.35, xmax  = no + 0.35, ymin=(response-2*se), ymax=(response+2*se), fill=signif), alpha = 0.5) +
  geom_segment(data = bin_hel_data, mapping=aes(x=no - 0.35, xend = no + 0.35, y=(response), yend=(response)), col="black", size=0.3) +
  scale_fill_manual(breaks = c(TRUE, FALSE), values=c(viridis(5)[4], viridis(5)[4])) +
  scale_x_continuous(breaks = bin_hel_data$no, labels = stri_trans_totitle(bin_hel_data$labels)) +
  scale_y_continuous(limits=c(-2,2.2), breaks=seq(-2,2, by=1), oob=scales::rescale_none, name="") +
  theme_bw() + partials_theme +
  theme(axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.x=element_text(colour="black", angle=90, hjust=1), legend.position="none",
        axis.title.y=element_blank(), axis.line.y = element_line(size = 0.5)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.line.y.left = element_line(size = 0.5, color = 'grey50'))

#if (SHOW_DEV_EXPL) bin_plot_vir = bin_plot_vir + geom_label(data = bin_hel_data, mapping=aes(x=no, y = response + 2*se + 0.4, label=dev_explained), color="black", family="Helvetica", size=1.5, label.size=0, fill="#FFFFFF8C") +


helminth_plots <- plot_grid(textGrob("Strength of effect on\nhelminth richness", 
                                     gp=gpar(fontface="bold", col="black", fontsize=10), rot=90),
                            plot_grid(plotlist = c(smooth_plots_helminth[1:5]), 
                                      nrow=1, 
                                      align="h", 
                                      rel_widths = c(1,1,1,1,1),
                                      labels=c("(m)", "(n)", "(o)", "(p)", ""), 
                                      label_size=10, 
                                      hjust=0),
                            bin_plot_hel,
                            nrow= 1, 
                            rel_widths = c(0.5, 5, 1.3), 
                            labels=c("", "", "(q)"), 
                            label_size=10, 
                            hjust=0)

plot(helminth_plots)
ggsave(paste0(base_dir, 'Figures/helminth_plot.png'), width = 14, height=5)





# COMBINE THE PLOTS
allplots = cowplot::plot_grid(vir_plots, bac_plots, helminth_plots, nrow=3)
plot(allplots)
ggsave(paste0(base_dir, 'Figures/GAMs-total-richness-bact-virus-helminth.png'), width = 14, height=10)

#plot(allplots)
# Save pdf
pdf(file= paste0(base_dir, 'Figures/GAMs-total-richness-bact-virus-helminth.pdf'), width = 14, height= 10)
allplots
dev.off()

