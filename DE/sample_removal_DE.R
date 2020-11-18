# Sample removing differential expression analysis

# Imports
library(tidyverse)
library(edgeR)
library(scales)
library(reshape2)
library(GO.db)
library(goseq)
library(GGally)

# Style
source("~/Git/ggplot_theme_Publication/ggplot_theme_Publication-2.R")

setwd("~/Git/Victor_MS_2018/Master/DE/")
load("sugarcane_counts_go.RData")
load("sample_removal_dataframe.RData")

differential_expression = function(selected_samples, strategy, complete = FALSE) {
  new_experiment = experiment %>%
    filter(replicate == strategy) %>%
    filter(samples %in% selected_samples)
  
  y = DGEList(counts = txi$counts[,which(colnames(txi$counts) %in% new_experiment$samples)],
              group = new_experiment$brix)
  keep = rowSums(cpm(y) > 1) >= 2
  y = y[keep, , keep.lib.sizes = FALSE]
  y = calcNormFactors(y, method = "TMM")
  
  model = model.matrix(~ new_experiment$brix + 0, data = y$samples)
  y = estimateDisp(y, model, robust = TRUE)
  fit = glmFit(y, design = model)
  
  VLBxVHB.HB.LB = makeContrasts(VLB - 1/3*(VHB + HB + LB), levels = factor(new_experiment$brix))
  VHBxHB.LB = makeContrasts(VHB - 1/2*(HB + LB), levels = factor(new_experiment$brix))
  HBxLB = makeContrasts(HB - LB, levels = factor(new_experiment$brix))
  
  contrasts = list("VLBxVHB.HB.LB" = VLBxVHB.HB.LB, "VHBxHB.LB" = VHBxHB.LB, "HBxLB" = HBxLB)
  
  result = c()
  result2 = data.frame()
  for (con in list("VLBxVHB.HB.LB", "VHBxHB.LB", "HBxLB")) {
    de = rep(-2, nrow(txi$length))
    
    lrt = glmLRT(fit, contrast = contrasts[[con]])
    contrast_de = decideTestsDGE(lrt, adjust.method = "BH")

    if (complete) {
      lrt$table$de = contrast_de
      table = as.data.frame(matrix(NA, nrow = nrow(txi$length), ncol = 5))
      rownames(lrt$table) = c()
      table[which(rownames(txi$length) %in% rownames(contrast_de)), ] = lrt$table
      result2 = rbind(result2, table)
    } else {
      de[which(rownames(txi$length) %in% rownames(contrast_de))] = contrast_de
      result = c(result, de)
    }
  }
  
  if (complete) {
    return(result2[2:nrow(result2),])
  } else {
    return(result)
  }
}

create_combinations = function(number, strategy){
  VLB = as.character(experiment$samples[which(experiment$brix == "VLB" & experiment$replicate == strategy)])
  LB = as.character(experiment$samples[which(experiment$brix == "LB" & experiment$replicate == strategy)])
  HB = as.character(experiment$samples[which(experiment$brix == "HB" & experiment$replicate == strategy)])
  VHB = as.character(experiment$samples[which(experiment$brix == "VHB" & experiment$replicate == strategy)])
  comb = expand.grid(2:3, 2:3, 2:3, 2:3)
  comb = comb[which(rowSums(comb) == number),]
  
  result = c()
  for (i in 1:nrow(comb)) {
    result = cbind(result, apply(expand.grid(
      combn(VLB, comb[i, 1], simplify = FALSE),
      combn(LB, comb[i, 2], simplify = FALSE),
      combn(HB, comb[i, 3], simplify = FALSE),
      combn(VHB, comb[i, 4], simplify = FALSE)
    ), 1, unlist))
  }
  return(t(result))
}

n = nrow(txi$counts)
n_deg = list()
reproduc_df = data.frame()

for (rep in list("SBC", "SBDG")) {
  original_data = suppressWarnings(differential_expression(create_combinations(12, rep), rep))
  n_deg[[rep]] = sum(abs(original_data) == 1)

  for (i in 11:8) {
    reproduc_data = apply(create_combinations(i, rep), 1, function(x) suppressWarnings(differential_expression(x, rep)))
    
    temp = as.data.frame(reproduc_data) %>%
      mutate(gene = rep(rownames(txi$counts), 3),
             up = rowSums(reproduc_data == 1),
             down = rowSums(reproduc_data == -1),
             notsig = rowSums(reproduc_data == 0),
             filtered = rowSums(reproduc_data == -2)) %>%
      dplyr::select(gene, up, down, notsig, filtered) %>%
      mutate(contrast = rep(c("VLBxVHB.HB.LB", "VHBxHB.LB", "HBxLB"), each = n), 
             removed = 12 - i, 
             strategy = rep) %>%
      mutate(matchesDE = ifelse(original_data == 1, up, 0) + ifelse(original_data == -1, down, 0),
             invertedDE = ifelse(original_data == -1, up, 0) + ifelse(original_data == 1, down, 0),
             notsigDE = ifelse(abs(original_data) == 1, notsig, 0),
             filteredDE = ifelse(abs(original_data) == 1, filtered, 0),
             matchesNS = ifelse(original_data == 0, notsig, 0),
             degNS = ifelse(original_data == 0, up + down, 0),
             filteredNS = ifelse(original_data == 0, filtered, 0)) %>%
      dplyr::select(gene, contrast, removed, strategy, matchesDE, invertedDE, 
             notsigDE, filteredDE, matchesNS, degNS, filteredNS) %>%
      mutate(isDE = abs(original_data) == 1,
             isNS = original_data == 0)
    
    reproduc_df = rbind(reproduc_df, temp)  
  }
}

n_deg = list()
for (rep in list("SBC", "SBDG")) {
  original_data = suppressWarnings(differential_expression(create_combinations(12, rep), rep))
  n_deg[[rep]] = sum(abs(original_data) == 1)
}

save(reproduc_df, file = "sample_removal_dataframe.RData")

# if the contrast VLB x VHB.HB.LB is called VLB x ALL, run the code below:
# reproduc_df$contrast[which(reproduc_df$contrast == "VLBxALL")] = "VLBxVHB.HB.LB"

# plot A
plotA_df = reproduc_df %>%
  filter(isDE | isNS) %>%
  mutate(combinations = matchesDE + invertedDE + notsigDE + filteredDE + matchesNS + degNS + filteredNS) %>%
  group_by(removed, strategy, isDE) %>%
  summarise(Matches = mean(ifelse(isDE, matchesDE / combinations, matchesNS / combinations)), 
            Inverted = mean(ifelse(isDE, invertedDE / combinations, NA)), 
            NotSig = mean(ifelse(isDE, notsigDE / combinations, NA)), 
            Filtered = mean(ifelse(isDE, filteredDE / combinations, filteredNS / combinations)),
            DEGs = mean(ifelse(isDE, NA, degNS / combinations))) %>%
  melt(id.vars = 1:3) %>%
  mutate(isDE = ifelse(isDE, "Originally DEG", "Originally NotSig")) %>%
  mutate(variable_strategy = paste0(variable, " ", strategy))

plotA_aux = data.frame(removed = rep(0, 20), 
		       value = rep(0, 20), 
		       strategy = rep(c("SBC", "SBDG", "SBDG", "SBC"), 5),  
		       isDE = rep(c("Originally DEG", "Originally NotSig"), 10),
		       variable = rep(c("Matches", "Inverted", "NotSig", "Filtered", "DEGs"), each = 4)) %>%
	    mutate(variable_strategy = paste0(variable, " ", strategy))
plotA_aux$value[which(plotA_aux$variable == "Matches")] = 1
plotA_df = rbind(plotA_aux, plotA_df)
plotA_df$variable = factor(plotA_df$variable, levels = c("Matches", "Inverted", "NotSig", "Filtered", "DEGs"))

ggplot(plotA_df) +
  geom_line(aes(x = removed, y = value, color = variable, group = variable_strategy, linetype = strategy), na.rm = TRUE) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Number of removed samples",
       y = "Percentage of genes",
       linetype = "Strategy",
       color = "Variables") +
  scale_color_manual(values = c("#d64a3b", "#e29d3e", "#80c34f", "#30acec", "#a666e1")) +
  theme_Publication(base_family = "Helvetica") +
  theme(legend.position = "right", legend.direction = "vertical") +
  facet_wrap(~isDE)

ggsave("removal.eps", width = 5.5, height = 4.5)

# plot B
plotB_df = reproduc_df %>%
  filter(isDE) %>%
  mutate(combinations = matchesDE + invertedDE + notsigDE + filteredDE) %>%
  group_by(matchesDE, strategy, removed, combinations) %>%
  summarise(n = n()) %>%
  mutate(n = ifelse(strategy == "SBC", n / n_deg[["SBC"]], n / n_deg[["SBDG"]]))

for (i in 1:4) {
  for (rep in list("SBC", "SBDG")) {
    plotB_df[sort(which(plotB_df$strategy == rep & plotB_df$removed == i),
                  decreasing = TRUE), 5] =
      cumsum(plotB_df$n[sort(which(plotB_df$strategy == rep &plotB_df$removed == i),
                             decreasing = TRUE)])
  }
  plotB_df$removed[which(plotB_df$removed == i)] = paste0(i, " sample(s) removed")
}

ggplot(plotB_df) + 
  geom_line(aes(x = (combinations - matchesDE)/combinations, y = n, group = strategy, color = strategy)) +
  scale_x_continuous(labels = percent_format(accuracy = 1), breaks = c(0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_color_manual(values = c("SBC" = "#50b9e1", "SBDG" = "#45b833")) +
  labs(x = "Mismatches in the identification of DEGs",
       y = "Cumulative percentage of original DEGs",
       color = "Strategy") +
  theme_Publication(base_family = "Helvetica") +
  theme(aspect.ratio = 1) +
  facet_wrap(~removed, ncol = 2)

ggsave("cumulative.eps", width = 6, height = 6)

plotB_df %>%
  filter(matchesDE == combinations)

# plot C
plotC_df = reproduc_df %>%
  filter(isDE) %>%
  mutate(combinations = matchesDE + invertedDE + notsigDE + filteredDE) %>%
  group_by(removed, strategy, contrast) %>%
  summarise(Matches = mean(matchesDE / combinations), 
            Inverted = mean(invertedDE / combinations), 
            NotSig = mean(notsigDE / combinations), 
            Filtered = mean(filteredDE / combinations)) %>%
  melt(id.vars = 1:3) %>%
  mutate(variable_strategy = paste0(variable, " ", strategy))

plotC_df$contrast = factor(plotC_df$contrast, levels = c("VLBxVHB.HB.LB", "VHBxHB.LB", "HBxLB"))

ggplot(plotC_df) +
  geom_line(aes(x = removed, y = value, color = variable, group = variable_strategy, linetype = strategy)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), breaks = pretty_breaks()) +
  labs(x = "Number of removed samples",
       y = "Percentage of genes",
       linetype = "Strategy",
       color = "Variables") +
  scale_color_manual(values = c("#d64a3b", "#e29d3e", "#80c34f", "#30acec")) +
  theme(aspect.ratio = 2) +
  guides(color = guide_legend(order = 1)) +
  facet_wrap(~contrast)
  
# plot D
sbc_de = reproduc_df %>%
  filter(strategy == "SBC") %>%
  mutate(sbc = matchesDE + invertedDE + degNS) %>%
  filter(sbc > 0) %>%
  group_by(gene, contrast, sbc, removed) %>% 
  summarise() %>%
  ungroup()

sbdg_de = reproduc_df %>%
  filter(strategy == "SBDG") %>%
  mutate(sbdg = matchesDE + invertedDE + degNS) %>%
  filter(sbdg > 0) %>%
  group_by(gene, contrast, sbdg, removed) %>%
  summarise() %>%
  ungroup()

plotD_df = merge(sbc_de, sbdg_de, by = c("gene", "contrast", "removed")) %>%
  group_by(sbc, sbdg, contrast, removed) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(sbc = as.integer(sbc), sbdg = as.integer(sbdg)) %>%
  complete(sbc, sbdg, contrast, removed) %>%
  mutate(n = replace_na(n, 0))

combinations = function(i) {nrow(create_combinations(12-i, "SBC"))}
for (i in 1:4) {
  plotD_df$sbc[which(plotD_df$removed == i)] = plotD_df$sbc[which(plotD_df$removed == i)] / combinations(i)
  plotD_df$sbdg[which(plotD_df$removed == i)] = plotD_df$sbdg[which(plotD_df$removed == i)] / combinations(i)
  sbc_de$sbc[which(sbc_de$removed == i)] = sbc_de$sbc[which(sbc_de$removed == i)] / combinations(i)
  sbdg_de$sbdg[which(sbdg_de$removed == i)] = sbdg_de$sbdg[which(sbdg_de$removed == i)] / combinations(i)
}

plotD_df$contrast = factor(plotD_df$contrast, levels = c("VLBxVHB.HB.LB", "VHBxHB.LB", "HBxLB"))
plotD_df$removed = paste0(plotD_df$removed, " sample(s) removed")

ggplot(plotD_df) +
  geom_raster(aes(x = sbc, y = sbdg, fill = log10(n+1)), hjust = 0, vjust = 0) + 
  scale_fill_gradient(low = "#FBFBFB", high = "#E52839", 
                      breaks = c(0, log10(11), log10(101), log10(1001)), 
                      labels = c(1, 10, 100, 1000)) + 
  labs(fill = "Number of DEGs", 
       x = "Fraction of significant DE tests in SBC combinations", 
       y = "Fraction of significant DE tests in SBDG combinations") +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) + 
  theme(panel.background = element_blank(),
        legend.position = "bottom") + 
  facet_grid(removed ~ contrast, scales = "free")

ggsave("heatmap.eps", width = 8.3, height = 11.7)

# plot D2 - Functional enrichment
colnames(sbc_de)[3] = "n"
colnames(sbdg_de)[3] = "n"
sbc_de$strategy = "SBC"
sbdg_de$strategy = "SBDG"

gene_de = rbind(sbc_de, sbdg_de) %>%
  group_by(gene, contrast, strategy) %>%
  summarise(filter = sum(n > 0.95)) %>%
  filter(filter == 4) %>%
  dplyr::select(gene, contrast, strategy)

gene_keep = reproduc_df %>%
  group_by(gene, strategy) %>%
  summarise(n = sum(matchesDE, invertedDE, notsigDE, matchesNS, degNS)) %>%
  filter(n > 0.05 * max(n)) %>%
  dplyr::select(gene, strategy)
gene_keep$gene = as.character(gene_keep$gene)

conserved_df = list()
for (con in list("VLBxVHB.HB.LB", "VHBxHB.LB", "HBxLB")) {
  conserved_df[[con]] = data.frame()
}

for (rep in list("SBC", "SBDG")) {
  cols = which(colnames(txi$counts) %in% experiment$samples[which(experiment$replicate == rep)])
  rows = which(rownames(txi$counts) %in% gene_keep$gene[gene_keep$strategy == rep])
  gene_counts = txi$counts[rows, cols]
  gene_lengths = rowSums(txi$length[rows, cols] * gene_counts / rowSums(gene_counts))
  
  for (con in list("VLBxVHB.HB.LB", "VHBxHB.LB", "HBxLB")) {
    genes = as.integer(names(gene_lengths) %in% gene_de$gene[which(gene_de$contrast == con & gene_de$strategy == rep)])
    names(genes) = names(gene_lengths)
    
    print(paste0(sum(gene_de$gene[which(gene_de$contrast == con & gene_de$strategy == rep)] %in% sugarcane_GO$gene), con, rep))
    pwf = nullp(genes, bias.data = gene_lengths, plot.fit = F)
    
    go_table = goseq(pwf, gene2cat = sugarcane_GO)
    go_table = go_table[which(!is.na(go_table[,6])),]
    go_FDR = p.adjust(go_table$over_represented_pvalue, method = "BH")
    
    enriched = go_FDR < 0.01
    if (sum(enriched) != 0) {
      conserved_df[[con]] = rbind(
        conserved_df[[con]],
        data.frame(Term = go_table$term[enriched],
                   Strategy = rep(rep, length(sum(enriched))),
                   Number = go_table$numDEInCat[enriched],
                   Ontology = go_table$ontology[enriched]))
    }
  }
}

plot_d2 = list()
for (con in list("VLBxVHB.HB.LB", "VHBxHB.LB", "HBxLB")) {
  plot_d2[[con]] = ggplot(conserved_df[[con]], aes(x = Term, y = Number, fill = Strategy)) +
    geom_bar(stat = "identity", position = position_dodge2(preserve = "single"), 
             width = length(levels(conserved_df[[con]]$Strategy))/3) + 
    coord_flip() + ylab("Number of genes") + ylim(c(0, 600)) +
    scale_fill_manual(values = c("SBC" = "#50b9e1", "SBDG" = "#45b833")) + xlab("") +
    theme_Publication()
}
ggmatrix(plot_d2, 3, 1, yProportions = c(length(unique(conserved_df[[1]]$Term)), 
                                         length(unique(conserved_df[[2]]$Term)), 
                                         length(unique(conserved_df[[3]]$Term))), 
         yAxisLabels = c("VLBxVHB.HB.LB", "VHBxHB.LB", "HBxLB"),
         legend = grab_legend(plot_d2[[1]])) + theme_Publication(base_family = "Helvetica")
ggsave("goremoval.eps", width = 6, height = 8)

# plot E
sbc_de = reproduc_df %>%
  filter(strategy == "SBC") %>%
  mutate(sbc = matchesDE + invertedDE + degNS) %>%
  group_by(gene, contrast, removed) %>% 
  summarise(sbc = sbc > 0)

sbdg_de = reproduc_df %>%
  filter(strategy == "SBDG") %>%
  mutate(sbdg = matchesDE + invertedDE + degNS) %>%
  group_by(gene, contrast, removed) %>%
  summarise(sbdg = sbdg > 0)

plotE_df = merge(sbc_de, sbdg_de, by = c("gene", "contrast", "removed")) %>%
  filter(sbc | sbdg) %>%
  group_by(sbc, sbdg, contrast, removed) %>%
  summarise(n = n()) %>%
  mutate(common = sbc + 2 * sbdg)

plotE_df$common = factor(plotE_df$common, labels = c("SBC only", "SBDG only", "Common"))

ggplot(plotE_df) +
  geom_bar(aes(x = removed, y = n, group = common, fill = common), stat = "identity", position = "dodge") +
  geom_text(aes(label = n, x = removed + (as.numeric(common)-1)*0.3, y = n), vjust = 4.5, hjust = 1.2 , size = 3) +
  scale_y_continuous(trans = "log10", labels = trans_format("log10", math_format(10^.x))) +
  coord_flip() +
  labs(fill = "DEGs groupment", 
       y = "Number of DEGs", 
       x = "Number of removed samples") +
  scale_fill_manual(values = c("SBC only" = "#50b9e1", "SBDG only" = "#45b833", "Common" = "#e25575")) +
  theme_Publication() +
  theme(aspect.ratio = 3) +
  facet_wrap(~contrast)

################################################################################################
#                                                                                              #
#                       Comparing only the different genotypes                                 #
#                                                                                              #
################################################################################################
library(VennDiagram)
library(grid)
library(GGally)

venn_triple = function(group1, group2, group3, labels, colors, two_layers = TRUE){
  grid.newpage()
  group1 = replace_na(group1, -2)
  group2 = replace_na(group2, -2)
  group3 = replace_na(group3, -2)
  draw.triple.venn(
    area1 = sum(abs(group1) == 1),
    area2 = sum(abs(group2) == 1),
    area3 = sum(abs(group3) == 1),
    n12 = sum(group1 == group2 & abs(group1) == 1),
    n13 = sum(group3 == group1 & abs(group3) == 1),
    n23 = sum(group2 == group3 & abs(group2) == 1),
    n123 = sum(group1 == group2 & group1 == group3 & abs(group1) == 1),
    category = labels,
    fill = colors,
    cex = rep(0.8, 7))
  if (two_layers) {
    draw.triple.venn(
      area1 = sum(abs(group1) == 1),
      area2 = sum(abs(group2) == 1),
      area3 = sum(abs(group3) == 1),
      n12 = sum(group1 == group2 & abs(group1) == 1),
      n13 = sum(group3 == group1 & abs(group3) == 1),
      n23 = sum(group2 == group3 & abs(group2) == 1),
      n123 = sum(group1 == group2 & group1 == group3 & abs(group1) == 1),
      category = labels,
      cex = rep(0.8, 7))
  }
}

exclusive_samples = c("6B_S21_L003_R1", "17B_S23_L003_R1", "1A_S16_L003_R1", "2B_S20_L003_R1", 
                      "10A_S18_L003_R1", "7A_S17_L003_R1",  "13A_S19_L003_R1", "16B_S22_L003_R1")

exclusive = suppressWarnings(differential_expression(exclusive_samples, "SBDG", TRUE))
sbdg = suppressWarnings(differential_expression(create_combinations(12, "SBDG"), "SBDG", TRUE))
sbc = suppressWarnings(differential_expression(create_combinations(12, "SBC"), "SBC", TRUE))

labels = c("", "", "")
colors = c("#50b9e1", "#45b833", "#ffa500")
for (i in 1:3) {
  pdf(paste0(i,".pdf"), width = 2, height = 2)
  venn_triple(sbc$V5[((i-1)*n):(i*n)], 
              sbdg$V5[((i-1)*n):(i*n)], 
              exclusive$V5[((i-1)*n):(i*n)], 
              labels, 
              colors)
  dev.off()
}

labels = c("SBC", "SBDG", "Exclusive")
na_index = which(!is.na(exclusive[,1]) & !is.na(sbdg[,1]) & !is.na(sbc[,1]) & 1:nrow(sbc) < n)
length(na_index)
exclusive = exclusive[na_index,]
sbdg = sbdg[na_index,]
sbc = sbc[na_index,]

plot_list = list()
rep = list(sbc, sbdg, exclusive)

for (i in 1:3) {
  hist_plot = hist(rep[[i]]$V1[which(abs(rep[[i]]$V5) == 1)], breaks = 20, plot = FALSE)
  hist_df = data.frame(x = hist_plot$mids, y = hist_plot$counts, de = "DEGs")
  hist_plot = hist(rep[[i]]$V1[which(rep[[i]]$V5 == 0)], breaks = hist_plot$breaks, plot = FALSE)
  hist_df = rbind(hist_df, data.frame(x = hist_plot$mids, y = hist_plot$counts, de = "NotSig"))
  hist_plot = ggplot(hist_df) +
    geom_bar(aes(x = x, y = y, fill = de), stat = "identity") +
    scale_fill_manual(values = c(colors[i], "#707070")) + 
    scale_x_continuous(breaks = c(-10, 0, 10), limits = c(-14, 14)) + 
    scale_y_continuous(limits = c(0, 14000)) +
    theme(aspect.ratio = 1)
  plot_list[[i*4-3]] = hist_plot
}

for (i in 1:2) {
  for (j in (i+1):3) {
    logfc_df = data.frame(x = rep[[i]]$V1, y = rep[[j]]$V1)
    logfc_df$de = "NotSig"
    logfc_df$de[which(abs(rep[[i]]$V5) == 1)] = labels[i]
    logfc_df$de[which(abs(rep[[j]]$V5) == 1)] = labels[j]
    logfc_df$de[which(abs(rep[[i]]$V5) == 1 & abs(rep[[j]]$V5) == 1)] = "Common"
    logfc_plot = ggplot(logfc_df) +
      geom_point(aes(x = x, y = y, color = de), size = 0.2, alpha = 0.25) +
      labs(x = labels[i], y = labels[j]) +
      scale_x_continuous(breaks = c(-10, 0, 10), limits = c(-14, 14)) + 
      scale_y_continuous(breaks = c(-10, 0, 10), limits = c(-14, 14)) + 
      scale_color_manual(values = c("SBC" = colors[1],
                                    "SBDG" = colors[2],
                                    "Exclusive" = colors[3],
                                    "Common" = "#e25575",
                                    "NotSig" = "#707070"))
    
    plot_list[[i+j*3-3]] = logfc_plot
    correlation = round(cor(rep[[i]]$V1, rep[[j]]$V1), digits = 4)
    text = paste0("Corr: ", correlation, "\nGenotypes in common: ", 8*i -4*j + 4)
    plot_list[[i*3+j-3]] = ggally_text(text) + 
      theme_grey() + theme(aspect.ratio = 1)
  }
}

DEGs_factor = factor(c("SBC", "SBDG", "Exclusive", "Common", "NotSig"), 
                     levels = c("SBC", "SBDG", "Exclusive", "Common", "NotSig"))

legend = ggplot(data.frame(x = 1:5, y = 1:5, DEGs = DEGs_factor)) + 
  geom_bar(aes(x = x, y = y, fill = DEGs), stat = "identity") + 
  scale_fill_manual(values = c("SBC" = colors[1],
                             "SBDG" = colors[2],
                             "Exclusive" = colors[3],
                             "Common" = "#e25575",
                             "NotSig" = "#707070")) + 
  theme_Publication(base_family = "Helvetica")

ggmatrix(plot_list, 3, 3, labels, labels, legend = grab_legend(legend), xlab = "log(FC)") +
  theme_Publication(base_family = "Helvetica") +
  theme(panel.grid.major = element_blank())

ggsave("cor.eps", width = 7, height = 7.5, device = cairo_ps)
