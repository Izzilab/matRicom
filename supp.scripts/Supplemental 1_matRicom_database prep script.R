### Network creation script ###

library(dplyr)
library(data.table)
library(KEGGgraph)
library(nichenetr)
library(OmnipathR)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dorothea)
library(igraph)

df1 <- as.data.frame(nichenetr::sig_network) #original nichenet sig_network

df2 <- as.data.frame(interactions <- import_omnipath_interactions() %>% as_tibble())
df2 <- df2[df2$is_stimulation==1,c(3,4)]
df2$source <- "ominpathr_stimulation"
df2$database <- "ominpathr"
names(df2)[1:2] <- names(df1)[1:2]

df3 <- as.data.frame(fread("all_data_14_10_22.tsv")) #data from SIGNOR V3.0 database, https://signor.uniroma2.it/downloads.php
df3 <- df3[df3$EFFECT %in% c("up-regulates activity",
                             "up-regulates quantity by expression",
                             "up-regulates quantity",
                             "up-regulates",
                             "up-regulates quantity by stabilization"),]
xx <- unique(names(as.list(org.Hs.egALIAS2EG)))
df3 <- df3[df3$ENTITYA%in%xx & df3$ENTITYB%in%xx,]
df3 <- df3[,c(1,5)]
df3$source <- "signor_stimulation"
df3$database <- "signor"
names(df3)[1:2] <- names(df2)[1:2]

new_sig <- distinct(bind_rows(df1,df2,df3))
saveRDS(new_sig,"new_sig_2.RDS")

df1 <- as.data.frame(nichenetr::lr_network) #original lr nichenet network

df4 <- as.data.frame(fread("matrixdb_FULL.tab",header = T)) #data from http://matrixdb.univ-lyon1.fr/, IMEX-extended database
df4 <- df4[df4$`Host organism(s)`=="taxid:9606(Homo sapiens)", ]
names(df4)[1] <- "ID(s) interactor A"
df4$`ID(s) interactor A` <- gsub("uniprotkb:","",df4$`ID(s) interactor A`)
df4$`ID(s) interactor B` <- gsub("uniprotkb:","",df4$`ID(s) interactor B`)
df4 <- df4[,c(1,2)]

df4.1 <- bitr(df4$`ID(s) interactor A`, fromType="UNIPROT", toType="SYMBOL", OrgDb="org.Hs.eg.db")
df4 <- distinct(merge(df4,df4.1,by.x="ID(s) interactor A",by.y="UNIPROT"))
df4.2 <- bitr(df4$`ID(s) interactor B`, fromType="UNIPROT", toType="SYMBOL", OrgDb="org.Hs.eg.db")
df4 <- distinct(merge(df4,df4.2,by.x="ID(s) interactor B",by.y="UNIPROT"))
df4 <- df4[,c(3,4)]
df4$source <- "matrixDB"
df4$database <- "matrixDB"
names(df4)[1:2] <- names(df2)[1:2]

p1 <- parseKGML2Graph("hsa04510.xml",expandGenes=TRUE) #data from KEGG focal adhesion, https://www.genome.jp/pathway/hsa04510
p2 <- parseKGML2Graph("hsa04512.xml",expandGenes=TRUE) #data from KEGG ECM-receptor interaction, https://www.genome.jp/pathway/hsa04512
p1 <- as.data.frame(as_edgelist(graph_from_graphnel(p1)))
p2 <- as.data.frame(as_edgelist(graph_from_graphnel(p2)))
kg <- bind_rows(p1,p2)
kg$V1 <- gsub("hsa:","",kg$V1)
kg$V2 <- gsub("hsa:","",kg$V2)
kg$order <- c(1:nrow(kg))
kg1.1 <- bitr(kg$V1, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
kg <- distinct(merge(kg,kg1.1,by.x="V1",by.y="ENTREZID"))
names(kg)[4] <- "SYMBOL.X"
kg1.1 <- bitr(kg$V2, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
kg <- distinct(merge(kg,kg1.1,by.x="V2",by.y="ENTREZID"))
kg <- distinct(kg[,c(4:5)])
names(kg) <- c("from","to")
kg$source <- "selected_KEGG"
kg$database <- "kegg"

new_lr <- distinct(bind_rows(df1,df4,kg))
saveRDS(new_lr,"new_lr_2.RDS")

df1 <- as.data.frame(nichenetr::gr_network) #original gr nichenet network

df2 <- as.data.frame(dorothea::dorothea_hs)
df2 <- df2[df2$mor==1,]
df2 <- df2[,c(1,3)]
df2 <- distinct(df2)
df2$source <- "dorothea"
df2$database <- "dorothea"
names(df2)[1:2] <- names(df1)[1:2]

new_gr <- distinct(bind_rows(df1,df2))
saveRDS(new_gr,"new_gr_2.RDS")

v <- nichenetr::optimized_source_weights_df
mean(v$weight) #average signaling network weight
v[v$source=="trrust",] #for dorothea, same weight as TRRUST

new_network_weights_df = data.frame(source =
                                    c("ominpathr","signor","matrixDB","dorothea"),
                                    weight = c(0.504,0.504,0.504,0.309))
new_source_weights_df = optimized_source_weights_df %>% bind_rows(new_network_weights_df)

weighted_networks = construct_weighted_networks(lr_network = new_lr, sig_network = new_sig, gr_network = new_gr, source_weights_df = new_source_weights_df)
weighted_networks = apply_hub_corrections(weighted_networks = weighted_networks, lr_sig_hub = hyperparameter_list$lr_sig_hub, gr_hub = hyperparameter_list$gr_hub)

saveRDS(weighted_networks,"new_weighted_networks_2.RDS")

ligands <- lapply(unique(new_lr$from),as.list)
comparative_ligand_target_matrix <- construct_ligand_target_matrix(weighted_networks = weighted_networks, ligands = ligands, algorithm = "PPR", damping_factor = hyperparameter_list$damping_factor, ltf_cutoff = hyperparameter_list$ltf_cutoff)

saveRDS(comparative_ligand_target_matrix,"comparative_ligand_target_matrix_2.RDS")

mat <- as.data.frame(fread("matrisome_complete.txt")) #This is from MatrisomeDB
mat <- subset(mat,mat$category=="CORE MATRISOME")
new_ligand_target_matrix <- comparative_ligand_target_matrix[,colnames(comparative_ligand_target_matrix)%in%unique(mat$gene)]

saveRDS(new_ligand_target_matrix,"new_ligand_target_matrix_2.RDS")

matricom_obj <- list()
matricom_obj$new_lr <- new_lr
matricom_obj$new_sig <- new_sig
matricom_obj$new_gr <- new_gr
matricom_obj$new_ligand_target_matrix <- new_ligand_target_matrix
save(matricom_obj,"matricom_obj.rda")
