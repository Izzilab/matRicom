### Performance evaluation script ###

library(dplyr)
library(nichenetr)
library(matRicom)
library(tidyverse)
library(purrr)
library(ggsci)

#since the matricom ligand-target matrix cannot be optimized via mlrMBO (see the "l" list below), we have to create an unoptimized ligand-target matrix from NicheNet
df1 <- as.data.frame(nichenetr::lr_network)
df2 <- as.data.frame(nichenetr::sig_network)
df3 <- as.data.frame(nichenetr::gr_network)
weighted_networks = construct_weighted_networks(lr_network = df1, sig_network = df2, gr_network = df3, source_weights_df = optimized_source_weights_df)
weighted_networks = apply_hub_corrections(weighted_networks = weighted_networks, lr_sig_hub = hyperparameter_list$lr_sig_hub, gr_hub = hyperparameter_list$gr_hub)
ligands <- lapply(unique(df1$from),as.list)
ligand_target_matrix <- construct_ligand_target_matrix(weighted_networks = weighted_networks, ligands = ligands, algorithm = "PPR", damping_factor = hyperparameter_list$damping_factor, ltf_cutoff = hyperparameter_list$ltf_cutoff)

settings = expression_settings_validation %>% lapply(convert_expression_settings_evaluation)
settings = settings %>% purrr::discard(~length(.$from) > 1)

performances = settings %>% lapply(evaluate_target_prediction, ligand_target_matrix) %>% bind_rows()
performances = performances %>% dplyr::select(-aupr, -auc_iregulon,-pearson_log_pval,-spearman_log_pval ,-sensitivity_roc, -specificity_roc) %>% gather(key = scorename, value = scorevalue, auroc:spearman)
performances$ligandmatrix <- "nichenet"
scorelabels = c(auroc="AUROC", aupr_corrected="AUPR (corrected)", auc_iregulon_corrected = "AUC-iRegulon (corrected)",pearson = "Pearson correlation", spearman = "Spearman's rank correlation",mean_rank_GST_log_pval = "Mean-rank gene-set enrichment")
scorerandom = c(auroc=0.5, aupr_corrected=0, auc_iregulon_corrected = 0, pearson = 0, spearman = 0,mean_rank_GST_log_pval = 0) %>% data.frame(scorevalue=.) %>% rownames_to_column("scorename")

aggregate(performances$scorevalue,list(performances$scorename),mean)
#auc_iregulon_corrected 0.04559323
#aupr_corrected 0.01873371
#auroc 0.54980611
#mean_rank_GST_log_pval 7.38866929
#pearson 0.06782953
#spearman 0.02784096

performances %>%
  mutate(model = "NicheNet") %>%
  ggplot() +
  geom_violin(aes(model, scorevalue, group=model, fill = model)) +
  geom_boxplot(aes(model, scorevalue, group = model),width = 0.05) +
  scale_fill_npg() +
  scale_y_continuous("Score target prediction") +
  facet_wrap(~scorename, scales = "free", labeller=as_labeller(scorelabels)) +
  geom_hline(aes(yintercept=scorevalue), data=scorerandom, linetype = 2, color = "red") +
  theme_bw() + ggtitle("NicheNet")

new_ligand_target_matrix <- matricom_obj$new_ligand_target_matrix
l <- list() #RETURNS EMPTY
for(i in names(settings)){
  z <- settings[[i]]
  if(z$from %in% colnames(new_ligand_target_matrix)){
    l[[i]] <- z}else{
      next
    }
}

comparative_ligand_target_matrix <- readRDS("comparative_ligand_target_matrix.RDS")
performances2 = settings %>% lapply(evaluate_target_prediction, comparative_ligand_target_matrix) %>% bind_rows()
performances2 = performances2 %>% dplyr::select(-aupr, -auc_iregulon,-pearson_log_pval,-spearman_log_pval ,-sensitivity_roc, -specificity_roc) %>% gather(key = scorename, value = scorevalue, auroc:spearman)
performances2$ligandmatrix <- "matricom"
scorelabels = c(auroc="AUROC", aupr_corrected="AUPR (corrected)", auc_iregulon_corrected = "AUC-iRegulon (corrected)",pearson = "Pearson correlation", spearman = "Spearman's rank correlation",mean_rank_GST_log_pval = "Mean-rank gene-set enrichment")
scorerandom = c(auroc=0.5, aupr_corrected=0, auc_iregulon_corrected = 0, pearson = 0, spearman = 0,mean_rank_GST_log_pval = 0) %>% data.frame(scorevalue=.) %>% rownames_to_column("scorename")

aggregate(performances2$scorevalue,list(performances2$scorename),mean)
#auc_iregulon_corrected 0.04066423
#aupr_corrected 0.01848569
#auroc 0.58758742
#mean_rank_GST_log_pval 16.82192424
#pearson 0.06543646
#spearman 0.03709467

performances2 %>%
  mutate(model = "matricom") %>%
  ggplot() +
  geom_violin(aes(model, scorevalue, group=model, fill = model)) +
  geom_boxplot(aes(model, scorevalue, group = model),width = 0.05) +
  scale_fill_jco() +
  scale_y_continuous("Score target prediction") +
  facet_wrap(~scorename, scales = "free", labeller=as_labeller(scorelabels)) +
  geom_hline(aes(yintercept=scorevalue), data=scorerandom, linetype = 2, color = "red") +
  theme_bw() + ggtitle("matRicom")

x <- data.frame(source=c("NicheNet","matRicom","NicheNet","matRicom"),
                variable=c("mean_rank_GST_log_pval","mean_rank_GST_log_pval","Spearman","Spearman"),
                value=c(7.388,16.821,0.027,0.037))
x$source <- factor(x$source,levels=c("NicheNet","matRicom"))
ggplot(x,aes(source,value,fill=source)) +
  geom_bar(stat="Identity") +
  scale_fill_manual(values=c("lavenderblush2","lightblue2")) +
  theme_bw() + xlab("") + ylab("performance indicator") +
  facet_wrap(~x$variable,scales="free")


performances3 = settings %>% lapply(evaluate_target_prediction, new_ligand_target_matrix) %>% bind_rows() #just for checking, Error in FUN(X[[i]], ...) : ligand should be in ligand_target_matrix

### AN ADDITIONAL EXAMPLE USING THE p-EMT SIGNATURE ###
library(nichenetr)
library(tidyverse)
library(matRicom)
library(eulerr)
library(ggplot2)
library(ggsci)

hnscc_expression <- readRDS(url("https://zenodo.org/record/3260758/files/hnscc_expression.rds"))
expression <- hnscc_expression$expression
sample_info <- hnscc_expression$sample_info 

tumors_remove <- c("HN10","HN","HN12", "HN13", "HN24", "HN7", "HN8","HN23")
CAF_ids <- sample_info %>% filter(`Lymph node` == 0 & !(tumor %in% tumors_remove) & `non-cancer cell type` == "CAF") %>% pull(cell)
malignant_ids <- sample_info %>% filter(`Lymph node` == 0 & !(tumor %in% tumors_remove) & `classified  as cancer cell` == 1) %>% pull(cell)
expressed_genes_CAFs <- expression[CAF_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()
expressed_genes_malignant <- expression[malignant_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()

#nichenet version, using the unoptimized matrix
df1 <- as.data.frame(nichenetr::lr_network) 
df2 <- as.data.frame(nichenetr::gr_network)
df3 <- as.data.frame(nichenetr::sig_network)
weighted_networks <- construct_weighted_networks(lr_network = df1, sig_network = df2, gr_network = df3, source_weights_df = optimized_source_weights_df)
weighted_networks <- apply_hub_corrections(weighted_networks = weighted_networks, lr_sig_hub = hyperparameter_list$lr_sig_hub, gr_hub = hyperparameter_list$gr_hub)
ligands <- lapply(unique(df1$from),as.list)
ligand_target_matrix <- construct_ligand_target_matrix(weighted_networks = weighted_networks, ligands = ligands, algorithm = "PPR", damping_factor = hyperparameter_list$damping_factor, ltf_cutoff = hyperparameter_list$ltf_cutoff)

pemt_geneset <- readr::read_tsv(url("https://zenodo.org/record/3260758/files/pemt_signature.txt"), col_names = "gene") %>% pull(gene) %>% .[. %in% rownames(ligand_target_matrix)] # only consider genes also present in the NicheNet model - this excludes genes from the gene list for which the official HGNC symbol was not used by Puram et al.
background_expressed_genes <- expressed_genes_malignant %>% .[. %in% rownames(ligand_target_matrix)]

lr_network <- nichenetr::lr_network
ligands <- lr_network %>% pull(from) %>% unique()
expressed_ligands <- intersect(ligands,expressed_genes_CAFs)
receptors <- lr_network %>% pull(to) %>% unique()
expressed_receptors <- intersect(receptors,expressed_genes_malignant)
potential_ligands <- lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

ligand_activities <- predict_ligand_activities(geneset = pemt_geneset, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
best_upstream_ligands <- ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)

k <- 10
n <- 3
pemt_gene_predictions_top20_list <- seq(n) %>% lapply(assess_rf_class_probabilities, folds = k, geneset = pemt_geneset, background_expressed_genes = background_expressed_genes, ligands_oi = best_upstream_ligands, ligand_target_matrix = ligand_target_matrix)
target_prediction_performances_discrete_fisher <- pemt_gene_predictions_top20_list %>% lapply(calculate_fraction_top_predicted_fisher, quantile_cutoff = 0.95) 
target_prediction_performances_discrete_fisher %>% unlist() %>% mean() #0.0397151

#matRicom version
#the matrix is available online at https://doi.org/10.5281/zenodo.7385244
comparative_ligand_target_matrix <- readRDS("comparative_ligand_target_matrix.RDS")
lr_network <- matricom_obj$new_lr
ligands <- lr_network %>% pull(from) %>% unique()
expressed_ligands <- intersect(ligands,expressed_genes_CAFs)
receptors <- lr_network %>% pull(to) %>% unique()
expressed_receptors <- intersect(receptors,expressed_genes_malignant)
potential_ligands <- lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

ligand_activities2 <- predict_ligand_activities(geneset = pemt_geneset, background_expressed_genes = background_expressed_genes, ligand_target_matrix = comparative_ligand_target_matrix, potential_ligands = potential_ligands)
best_upstream_ligands2 <- ligand_activities2 %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)

k <- 10
n <- 3
pemt_gene_predictions_top20_list2 <- seq(n) %>% lapply(assess_rf_class_probabilities, folds = k, geneset = pemt_geneset, background_expressed_genes = background_expressed_genes, ligands_oi = best_upstream_ligands2, ligand_target_matrix = comparative_ligand_target_matrix)
target_prediction_performances_discrete_fisher2 <- pemt_gene_predictions_top20_list2 %>% lapply(calculate_fraction_top_predicted_fisher, quantile_cutoff = 0.95) 
target_prediction_performances_discrete_fisher2 %>% unlist() %>% mean() #5.586506e-15

top_predicted_genes <- seq(length(pemt_gene_predictions_top20_list)) %>% lapply(get_top_predicted_genes,pemt_gene_predictions_top20_list) %>% reduce(full_join, by = c("gene","true_target"))
top_predicted_genes <- top_predicted_genes %>% filter(true_target)
top_predicted_genes <- unique(top_predicted_genes$gene)

top_predicted_genes2 <- seq(length(pemt_gene_predictions_top20_list)) %>% lapply(get_top_predicted_genes,pemt_gene_predictions_top20_list2) %>% reduce(full_join, by = c("gene","true_target"))
top_predicted_genes2 <- top_predicted_genes2 %>% filter(true_target)
top_predicted_genes2 <- unique(top_predicted_genes2$gene)

#graphs
eu <- euler(list(nichenet=top_predicted_genes,
                 matricom=top_predicted_genes2))
plot(eu,
     quantities=T,
     shape="ellipse",
     fills=c("lavenderblush2","lightblue2"),
     labels=c("p-EMT predictors in NicheNet","p-EMT predictors in matRicom"))

d1 <- ifelse(top_predicted_genes %in% mat$gene,"matrisome","not matrisome")
d2 <- ifelse(top_predicted_genes2 %in% mat$gene,"matrisome","not matrisome")
d1 <- unlist(table(d1))
d2 <- unlist(table(d2))
f <- data.frame(v=c(d1,d2),
                s=c(rep("NicheNet",2),rep("matRicom",2)),
                cl=rep(c("matrisome","not matrisome"),2))
f$s <- factor(f$s,levels = c("NicheNet","matRicom"))
ggplot(f,aes(cl,v,fill=s)) +
  geom_bar(stat="identity",position="fill") +
  scale_fill_manual(values = c("lavenderblush2","lightblue2")) +
  theme_bw() + 
  xlab("") + ylab("% composition of top-20 \np-EMT predicted genes")

p1 <- -log10(target_prediction_performances_discrete_fisher %>% unlist() %>% mean())
p2 <- -log10(target_prediction_performances_discrete_fisher2 %>% unlist() %>% mean())
f <- data.frame(v=c(p1,p2),
                s=c("NicheNet","matRicom"))
f$s <- factor(f$s,levels = c("NicheNet","matRicom"))
ggplot(f,aes(s,v,fill=s)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("lavenderblush2","lightblue2")) +
  ylim(c(0,15)) +
  theme_bw() + 
  xlab("") + ylab("top-20 p-EMT predicted \ntarget genes (-log10 p-value)")
