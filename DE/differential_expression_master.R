# Genic differential expression analysis

# Imports
library(tximport)
library(edgeR)
library(ggplot2)
library(goseq)
library(GO.db)
library(VennDiagram)
library(dplyr)
library(rlist)
library(GGally)

# Style
extrafont::font_import()
extrafont::loadfonts(device = "postscript")
source("~/Git/ggplot_theme_Publication/ggplot_theme_Publication-2.R")

# Getting the data
folders = list.files("../Quantification/Default_useEM/quants_out/")
load("sugarcane_counts_go.RData")
  
genotypes = factor(
  c("RB92579","White Mauritius","F36-819","F36-819",
    "RB835486","F36-819","Krakatau","SP80-3280",
    "R570","Criolla Rayada","SP80-3280","SP80-3280",
    "R570","IN84-58","IN84-58","R570",
    "R570","IN84-58","IJ76-317","IN84-58",
    "F36-819","SES205A","White Transparent","SP80-3280"),
  levels = c("IN84-58","F36-819","R570","SP80-3280",
    "Krakatau","Criolla Rayada","White Transparent","White Mauritius",
    "SES205A","IJ76-317","RB92579","RB835486"))

brix = factor(
  c("HB", "VHB", "LB","LB", "VHB", "LB","VLB", "VHB", "HB","LB", "VHB", "VHB",
    "HB", "VLB", "VLB","HB", "HB", "VLB","LB", "VLB", "LB","VLB", "HB", "VHB"), 
  levels = c("VLB", "LB", "HB", "VHB"))

replicate = rep("SBC", 24)
replicate[which(nchar(folders) > 10)] = "SBDG"
replicate = factor(replicate)

brix_palette = c("dodgerblue", "limegreen", "goldenrod2", "brown1")

experiment = data.frame(samples = folders, genotypes = genotypes, brix = brix, replicate = replicate)

files = file.path("../Quantification/Default_useEM/quants_out", folders, "quant.sf")
names(files) = folders
tx2gene = read.table("../Assembly/Assembly_default/default.Trinity.fasta.map")
txi = tximport(files, type = "salmon", tx2gene = tx2gene)

# Physiological data
read.csv("../Physiologic/data.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE) %>%
  filter(Geno %in% casefold(experiment$genotypes, upper = TRUE)) %>%
  mutate(Geno = factor(Geno)) %>%
  group_by(Geno) %>%
  summarise(mean = mean(BRIX_Lab), sd = sd(BRIX_Lab), n = n()) %>%
  arrange(mean, desc = T)

# Delete low expression genes and normalization
y = list()
keep = list()
for (rep in list("SBC", "SBDG")) {
  y[[rep]] = DGEList(counts = txi$counts[,which(replicate == rep)],
                    group = brix[which(replicate == rep)])
  keep[[rep]] = rowSums(cpm(y[[rep]]) > 1) >= 3
  y[[rep]] = y[[rep]][keep[[rep]], , keep.lib.sizes = FALSE]
  y[[rep]] = calcNormFactors(y[[rep]])
}

# Common genes after filter
common_kept = rownames(y$SBC$counts)[rownames(y$SBC$counts) %in% rownames(y$SBDG$counts)]
pdf("venn.pdf", width = 4, height = 3.5)
grid.newpage()
draw.pairwise.venn(area1 = nrow(y$SBC$counts), area2 = nrow(y$SBDG$counts),
                   cross.area = length(common_kept), 
                   fill = c("#50b9e1", "#45b833"),
                   ext.dist = rep(0.05, 2), 
                   cex = rep(0.8, 3))
dev.off()

# MDS 
mds_plot = list()
for (rep in list("SBC", "SBDG")) {
  mds = plotMDS(y[[rep]], plot = F, top = 2000)
  mds_df = data.frame(x = mds$x, y = mds$y, genotypes = genotypes[which(replicate == rep)], 
                      brix = brix[which(replicate == rep)])
  mds_plot[[rep]] = ggplot(mds_df, aes(x, y, size = 5, shape = genotypes, color = brix)) + 
    geom_point() +
    guides(size = F, shape = guide_legend(order = 1)) + 
    scale_shape_manual(values=seq(65,65+12)) + 
    scale_color_manual(values = brix_palette) +
    lims(x = c(-2, 4), y = c(-2, 4)) +
    labs(color = "Brix", shape = "Genotypes")
}

ggmatrix(mds_plot, 2, 1, yAxisLabels = c("SBC", "SBDG"), 
         legend = grab_legend(mds_plot[[2]])) + 
  theme_Publication(base_family = "Helvetica") +
  theme(legend.position = "right") +
  labs(x = "logFC dim 1", y = "logFC dim 2")

ggsave(filename = "mds.eps", width = 7, height = 6)


# Estimate dispersion
model = list()
fit = list()
par(mfrow = c(1,2))
for (rep in list("SBC", "SBDG")) {
  model[[rep]] = model.matrix(~ experiment$brix[which(replicate == rep)] + 0, data=y[[rep]]$samples)
  y[[rep]] = estimateDisp(y[[rep]], model[[rep]], robust = T)
  plotBCV(y[[rep]])
  fit[[rep]] = glmFit(y[[rep]], design = model[[rep]])
}

# Comparison of BCVs
sbc_disp = y$SBC$tagwise.dispersion[which(rownames(y$SBC$counts) %in% common_kept)]
sbdg_disp = y$SBDG$tagwise.dispersion[which(rownames(y$SBDG$counts) %in% common_kept)]
mean(sbdg_disp > sbc_disp)

# Creating contrasts 
VLBxVHB.HB.LB = makeContrasts(VLB - 1/3*(VHB + HB + LB), levels = factor(brix))
VHBxHB.LB = makeContrasts(VHB - 1/2*(HB + LB), levels = factor(brix))
HBxLB = makeContrasts(HB - LB, levels = factor(brix))

contrasts = list("VLBxVHB.HB.LB" = VLBxVHB.HB.LB, "VHBxHB.LB" = VHBxHB.LB, "HBxLB" = HBxLB)

# Differential expression
DE = list()
de_df = data.frame(gene = c(), logFC = c(), logCPM = c(), DE = c(), contrast = c(), strategy = c())

for (rep in list("SBC", "SBDG")) {
  for (con in list("VLBxVHB.HB.LB", "VHBxHB.LB", "HBxLB")) {
  lrt = glmLRT(fit[[rep]], contrast = contrasts[[con]])
  contrast_de = decideTestsDGE(lrt, adjust.method = "BH")
  #print(topTags(lrt, n = 5, adjust.method = "BH", sort.by = "p.value"))
  DE[[paste(rep, con)]] = contrast_de
  
  lrt_df = as.data.frame(lrt$table) %>%
    mutate(gene = rownames(.),
           DE = as.numeric(contrast_de),
           contrast = con,
           strategy = rep) %>%
    dplyr::select(-c(LR, PValue))
  
  lrt_df$DE[which(lrt_df$DE == 0)] = "NotSig"
  lrt_df$DE[which(lrt_df$DE == 1)] = "Up"
  lrt_df$DE[which(lrt_df$DE == -1)] = "Down"
  de_df = rbind(de_df, lrt_df)
  
  print(c(sum(contrast_de == 1), sum(contrast_de == -1), length(contrast_de)))
  }
}

de_df$contrast = factor(de_df$contrast, levels = c("VLBxVHB.HB.LB", "VHBxHB.LB", "HBxLB"))
de_df$DE = factor(de_df$DE, levels = c("Up", "NotSig", "Down"))

ggplot(de_df) +
  geom_point(aes(x = logCPM, y = logFC, color = DE), size = 0.3) +
  scale_color_manual(values = c("red", "black", "green")) +
  facet_grid(strategy ~ contrast) +
  theme_Publication(base_family = "Helvetica")

ggsave(filename = "ma.tiff", width = 8, height = 5)
de_df = de_df[,c(3,2,1,4,5,6)]
write.csv(de_df, "FileS1.csv", row.names = F, col.names = T)

# Gene Ontology annotation
sugarcane_GO_file = read.table("../Annotation/go_annotation.txt", col.names = c("gene", "GO"), stringsAsFactors = F)
sugarcane_GO = data.frame()

for (i in 1:nrow(sugarcane_GO_file)) {
  go_list = sugarcane_GO_file[i,2]
  for (go in strsplit(go_list, ",")) {
    sugarcane_GO = rbind(sugarcane_GO, data.frame(gene = sugarcane_GO_file[i,1], GO = go))
  }
}

annotation_GO = sugarcane_GO %>%
  group_by(GO) %>%
  summarise(n = n()) %>%
  mutate(Term = Term(as.character(GO)),
         Ontology = Ontology(as.character(GO))) %>%
  arrange(-n)

annotation_GO %>%
  group_by(Ontology) %>%
  summarise(q = sum(n)) %>%
  ggplot() + geom_bar(aes(x = "", y = q/sum(q), fill = Ontology), width = 1, stat = "identity") +
  coord_polar(theta = "y", start = 3) + theme_void()


# Terms in assembly
length(unique(sugarcane_GO$gene))
length(keep[["SBC"]])
sum(unique(sugarcane_GO$gene) %in% names(which(keep[["SBC"]])))
sum(unique(sugarcane_GO$gene) %in% names(which(keep[["SBDG"]])))


# Calculating weighted lenghts
length_weighted = list()
for (rep in list("SBC", "SBDG")) {
  cols = which(colnames(txi$counts) %in% experiment$samples[which(experiment$replicate == rep)])
  count_keep = txi$counts[keep[[rep]], cols]
  length_weighted[[rep]] = txi$length[keep[[rep]], cols]
  length_weighted[[rep]] = rowSums(length_weighted[[rep]] * count_keep / rowSums(count_keep))
}

# Funcional enrichment
enrichment_df = list()
go = data.frame(GO = c(), pvalue = c(), numDEInCat = c(), numInCat = c(), term = c(), contrast = c(), strategy = c())
for (con in list("VLBxVHB.HB.LB", "VHBxHB.LB", "HBxLB")) {
  enrichment_df[[con]] = data.frame()
}

for (rep in list("SBC", "SBDG")) {
  for (con in list("VLBxVHB.HB.LB", "VHBxHB.LB", "HBxLB")) {
    genes = as.integer(rownames(y[[rep]]) %in% rownames(y[[rep]])[as.logical(DE[[paste(rep, con)]])])
    names(genes) = rownames(y[[rep]])
    
    pwf = nullp(genes, bias.data = length_weighted[[rep]], plot.fit = F)
    
    go_table = goseq(pwf, gene2cat = sugarcane_GO)
    go_table = go_table[which(!is.na(go_table[,6])),]
    go_FDR = p.adjust(go_table$over_represented_pvalue, method = "BH")
    
    enriched = go_FDR < 0.01
    if (sum(enriched) != 0) {
      enrichment_df[[con]] = rbind(
        enrichment_df[[con]],
        data.frame(Term = go_table$term[enriched],
                   Strategy = rep(rep, length(sum(enriched))),
                   Number = go_table$numDEInCat[enriched],
                   Ontology = go_table$ontology[enriched]))
    }
    
    go_table$FDR = go_FDR
    go_table$contrast = con
    go_table$strategy = rep
    go_table = go_table[,-c(2,3,7)]
    
    go = rbind(go, go_table)
  }
}

write.csv(go, "FileS2.csv")

enrichment_plot = list()
for (con in list("VLBxVHB.HB.LB", "VHBxHB.LB", "HBxLB")) {
  enrichment_plot[[con]] = ggplot(enrichment_df[[con]], aes(x = Term, y = Number, fill = Strategy)) +
    geom_bar(stat = "identity", position = position_dodge2(preserve = "single"), 
             width = length(levels(enrichment_df[[con]]$Strategy))/3) + 
    coord_flip() + theme(legend.position="top") + ylab("Number of genes") + ylim(c(0, 450)) +
    scale_fill_manual(values = c("SBC" = "#50b9e1", "SBDG" = "#45b833")) + xlab("")
}

ggmatrix(enrichment_plot, 3, 1, yProportions = c(length(unique(enrichment_df[[1]]$Term)), 
                                                 length(unique(enrichment_df[[2]]$Term)), 
                                                 length(unique(enrichment_df[[3]]$Term))), 
         yAxisLabels = c("VLBxVHB.HB.LB", "VHBxHB.LB", "HBxLB"), 
         legend = grab_legend(enrichment_plot[[1]])) + 
  theme_Publication(base_family = "Helvetica")
  
ggsave(filename = "go.eps", height = 8, width = 6)

# Common DE genes
for (con in list("VLBxVHB.HB.LB", "VHBxHB.LB", "HBxLB")) {
  pdf(paste0(con, ".pdf"), width = 2, height = 2)
  grid.newpage()
  deg_sbc <- rep(0, nrow(txi$counts)) -> deg_sbdg
  deg_sbc[keep[["SBC"]]] = DE[[paste("SBC", con)]] != 0
  deg_sbdg[keep[["SBDG"]]] = DE[[paste("SBDG", con)]] != 0
  
  draw.pairwise.venn(area1 = sum(deg_sbc), area2 = sum(deg_sbdg),
                     cross.area = sum(deg_sbc & deg_sbdg), 
                     fill = c("#50b9e1", "#45b833"),
                     ext.dist = rep(0.05, 2), 
                     cex = rep(0.8, 3))
  dev.off()
}

sum(DE[["SBC VLBxVHB.HB.LB"]] != 0 | DE[["SBC VHBxHB.LB"]] != 0 | DE[["SBC HBxLB"]] != 0)
sum(keep[["SBC"]])
sum(DE[["SBDG VLBxVHB.HB.LB"]] != 0 | DE[["SBDG VHBxHB.LB"]] != 0 | DE[["SBDG HBxLB"]] != 0)
sum(keep[["SBDG"]])


sbc_genes = rownames(y$SBC$counts)
sbdg_genes = rownames(y$SBDG$counts)
common_genes = sbc_genes[which(sbc_genes %in% sbdg_genes)]

var_df = data.frame(genes = common_genes,
                    all_sbc = apply(cpm(y$SBC)[which(sbc_genes %in% common_genes),], 1, function(x) var(x)/mean(x)),
                    all_sbdg = apply(cpm(y$SBDG)[which(sbdg_genes %in% common_genes),], 1, function(x) var(x)/mean(x)),
                    group_sbc = apply(cpmByGroup(y$SBC)[which(sbc_genes %in% common_genes),], 1, function(x) var(x)/mean(x)),
                    group_sbdg = apply(cpmByGroup(y$SBDG)[which(sbdg_genes %in% common_genes),], 1, function(x) var(x)/mean(x)))



var_df$de = "Not DE"
sbc_DE = rownames(DE[[1]])[which(DE[["SBC VLBxVHB.HB.LB"]] != 0 | DE[["SBC VHBxHB.LB"]] != 0 | DE[["SBC HBxLB"]] != 0)]
sbdg_DE = rownames(DE[[1]])[which(DE[["SBDG VLBxVHB.HB.LB"]] != 0 | DE[["SBDG VHBxHB.LB"]] != 0 | DE[["SBDG HBxLB"]] != 0)]
common_DE = sbc_DE[which(sbc_DE %in% sbdg_DE)]

var_df$de[which(common_genes %in% sbc_DE)] = "SBC DE"
var_df$de[which(common_genes %in% sbdg_DE)] = "SBDG DE"
var_df$de[which(common_genes %in% common_DE)] = "Both DE"

ggplot(var_df) +
  geom_point(aes(x = log(all_sbc), y = log(all_sbdg)), size = 0.2) +
  geom_abline(color = "red") + xlab("var(SBDG)") + ylab("var(SBC)") + 
  theme_Publication() +
  theme(aspect.ratio = 1)

var_df_sbc = var_df
var_df_sbc$de[which(var_df_sbc$de == "SBDG DE")] = "Not DE"
var_df_sbc$de[which(var_df_sbc$de == "Both DE")] = "SBC DE"

ggplot(var_df_sbc) +
  geom_point(aes(x = group_sbc/all_sbc, y = group_sbdg/all_sbdg, color = de), size = 0.2) + facet_wrap(~ de)

var_df_sbdg = var_df
var_df_sbdg$de[which(var_df_sbdg$de == "SBC DE")] = "Not DE"
var_df_sbdg$de[which(var_df_sbdg$de == "Both DE")] = "SBDG DE"

ggplot(var_df_sbdg) +
  geom_point(aes(x = group_sbc/all_sbc, y = group_sbdg/all_sbdg, color = de), size = 0.2) + facet_wrap(~ de)

