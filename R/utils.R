#' matricom: cell-ECM communications in scRNAseq data

#' @param seurat.obj character. A Seurat or Seurat-like object.
#' @param group.column character. Name of the column in metadata with cell identities to be used.
#' @param min.pct numeric. Minimum % of cells with non-zero expression of a ligand/receptor/target. All ligands/receptors/targets less present will be ignored. Only for ligands, a total of min.pct is also allowed from the total, on top of cell-specific quantities. Default is 0.3 (30%).
#' @param expr.filter numeric. Mean expression of a ligand/receptor/target. All ligands/receptors/targets less expressed will be ignored. Default is 1.
#' @param target character. A receiver cell type to focus on. If selected, only communications affecting that cell type are considered. Default is NULL, all cell types are considered.
#' @param geneset character. A response geneset (a vector of gene names matching those of the counts in the object) to test ligands against. If NULL (default), a recursive cell type-vs.-all DEG search is launched using Seurat functions.
#' @param max.ligands numeric. Maximum number of ligands to test for response. Default is 10, in decreasing Pearson order.
#' @param max.targets numeric. Maximum number of targets to continue with after having tested the response. Default is 10, in decreasing gene expression order.
#' @param explainable logical. Whether to apply differential correlation analysis to targets' regulators. In the analysis, all known regulators of each target gene per receiving population are correlated to their targets with Spearman correlation. In parallel, the same tests are repeated across the rest of sample (non-receiving population). Results are confronted and the regulator with the highest differential (positive) correlation in the receiving population is chosen as the master regulator and reported. Default is TRUE.
#' @param add.RF logical. Whether to add a random forest run to the results from the correlation analysis. Features are ordered by their VarImp and as many as max.ligands (up to twice that value, once for receptors and once for TF) are chosen. Default is TRUE.
#' @param scale.results logical. Whether to rescale the gene expression values of ligands/receptors/targets to the [0.01-1] interval by cell type to ease cross-comparisons. Default is TRUE.
#' @param use.ensembl logical. Whether to try and convert all matricom objects and complex annotation to ENSEMBLE GENE IDs on the fly to extend compatibility. It is not guaranteed to work.
#' @param verbose logical. Progress messages are printed to the console. Default is TRUE.

#' @return a dataframe with ligand/receptor/TF/target quadruplets, the sender and receiver cells for each, and the mean expression values of each triplet.
#' @export

#' @examples matricom(obj,"group.column")

matricom <- function(seurat.obj=NULL,
                     group.column=NULL,
                     min.pct=0.30,
                     expr.filter=1,
                     target=NULL,
                     geneset=NULL,
                     max.ligands=10,
                     max.targets=10,
                     explainable=TRUE,
                     add.RF=TRUE,
                     scale.results=TRUE,
                     use.ensembl=FALSE,
                     force.flushing=FALSE,
                     verbose=T){

  #error handling
  if(is.null(seurat.obj)){
    cat(crayon::red("no data provided, execution stops \n"))
    stop()
  }
  if(class(seurat.obj)!="Seurat"){
    cat(crayon::red("data are not a valid Seurat object, execution stops \n"))
    stop()
  }
  if(nrow(seurat.obj@meta.data)<1){
    cat(crayon::red("data are not a valid Seurat object, execution stops \n"))
    stop()
  }
  if(nrow(seurat.obj@assays$RNA@scale.data)<1){
    cat(crayon::red("data must be scaled, execution stops \n"))
    stop()
  }
  if(is.null(group.column)){
    cat(crayon::red("no group assignment column provided, execution stops \n"))
    stop()
  }
  if(length(intersect(colnames(seurat.obj@meta.data),group.column))<1){
    cat(crayon::red("group assignment column missing in object metadata slot, execution stops \n"))
    stop()
  }

  #read in the matricom objects
  if(isTRUE(verbose)){cat(crayon::white("(1/9) loading matRicom objects... "))}
  data(matricom_obj, package = "matRicom")
  data(ncpl, package = "matRicom")
  if(isFALSE(use.ensembl)){
    new_lr <- matricom_obj$new_lr
    new_gr <- matricom_obj$new_gr
    new_ligand_target_matrix <- matricom_obj$new_ligand_target_matrix}else{

      new_lr <- matricom_obj$new_lr
      new_gr <- matricom_obj$new_gr
      new_ligand_target_matrix <- matricom_obj$new_ligand_target_matrix

      require("AnnotationDbi")
      require("org.Hs.eg.db")

      new_lr$from.E <- suppressMessages(mapIds(org.Hs.eg.db, keys=new_lr$from, column="ENSEMBL", keytype="SYMBOL", multiVals="first"))
      new_lr$to.E <- suppressMessages(mapIds(org.Hs.eg.db, keys=new_lr$to, column="ENSEMBL", keytype="SYMBOL", multiVals="first"))
      new_lr$from <- NULL
      new_lr$to <- NULL
      new_lr <- new_lr[,c(3,4,1,2)]
      names(new_lr)[1:2] <- c("from","to")
      new_lr <- na.omit(new_lr)

      new_gr$from.E <- suppressMessages(mapIds(org.Hs.eg.db, keys=new_gr$from, column="ENSEMBL", keytype="SYMBOL", multiVals="first"))
      new_gr$to.E <- suppressMessages(mapIds(org.Hs.eg.db, keys=new_gr$to, column="ENSEMBL", keytype="SYMBOL", multiVals="first"))
      new_gr$from <- NULL
      new_gr$to <- NULL
      new_gr <- new_gr[,c(3,4,1,2)]
      names(new_gr)[1:2] <- c("from","to")
      new_gr <- na.omit(new_gr)

      vcn <- suppressMessages(mapIds(org.Hs.eg.db, keys=colnames(new_ligand_target_matrix), column="ENSEMBL", keytype="SYMBOL", multiVals="first"))
      vrn <- suppressMessages(mapIds(org.Hs.eg.db, keys=rownames(new_ligand_target_matrix), column="ENSEMBL", keytype="SYMBOL", multiVals="first"))
      colnames(new_ligand_target_matrix) <- vcn
      rownames(new_ligand_target_matrix) <- vrn
      new_ligand_target_matrix <- new_ligand_target_matrix[!is.na(rownames(new_ligand_target_matrix)),
                                                           !is.na(colnames(new_ligand_target_matrix))]
    }

  if(isTRUE(verbose)){cat(crayon::green("done \n"))}

  #isolation of valid ECM genes (should be expressed at least by 30% of a population but can be changed via the min.pct parameter)
  if(isTRUE(verbose)){cat(crayon::white("(2/9) evaluating active matrisome genes in cells... "))}

  obj <- seurat.obj
  seurat.obj <- NULL

  Idents(obj) <- obj@meta.data[,group.column]

  if(length(intersect(colnames(new_ligand_target_matrix),
                      rownames(obj@assays$RNA@counts)))<1){
    cat(crayon::red("\nno gene in data match with possible ligands, execution stops \n"))
    cat(crayon::red("\ntry checking the format of genes in your data \n"))
    stop()
  }

  sender_celltypes <- unique(obj@meta.data[,group.column])
  list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, obj, min.pct)
  names(list_expressed_genes_sender) <- sender_celltypes
  for(i in names(list_expressed_genes_sender)){
    s <- obj@meta.data[,group.column]
    names(s) <- rownames(obj@meta.data)
    s <- names(s[s%in%i])
    z <- list_expressed_genes_sender[i][[1]]
    k <- as.matrix(obj@assays$RNA@counts[rownames(obj@assays$RNA@counts)%in%z,
                                         colnames(obj@assays$RNA@counts)%in%s])
    if(nrow(k)<2){next}else{
      if(ncol(k)<2){next}else{
        k <- base::rowMeans(k)
        k <- k[k>=expr.filter]
        list_expressed_genes_sender[[i]] <- names(k)}}
  }
  expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
  expressed_genes_sender <- expressed_genes_sender[expressed_genes_sender%in%colnames(new_ligand_target_matrix)]
  common_genes_sender <- apply(obj@assays$RNA@counts[rownames(obj@assays$RNA@counts)%in%colnames(new_ligand_target_matrix),],1,function(x){
    ifelse((length(x[x>0])/length(x))>=min.pct,"1","NA")
  })
  common_genes_sender <- names(common_genes_sender[common_genes_sender!="NA"])
  k <- as.data.frame(obj@assays$RNA@counts[rownames(obj@assays$RNA@counts)%in%common_genes_sender,])
  k <- rowMeans(as.matrix(k))
  k <- k[k>=expr.filter]
  common_genes_sender <- common_genes_sender[common_genes_sender%in%names(k)]

  expressed_genes_sender <- unique(c(expressed_genes_sender,common_genes_sender))

  if(length(expressed_genes_sender)<1){
    cat(crayon::red("\nno gene passes the threshold, execution stops \n"))
    cat(crayon::red("\ntry changing the min.pct or the expr.filter parameter \n"))
    stop()
  }

  ligands <- new_lr %>% pull(from) %>% unique()
  expressed_ligands <- intersect(ligands,expressed_genes_sender)

  if(length(expressed_ligands)<1){
    cat(crayon::red("\nno gene passes the threshold, execution stops \n"))
    cat(crayon::red("\ntry changing the min.pct or the expr.filter parameter \n"))
    stop()
  }
  if(isTRUE(verbose)){cat(crayon::green("done \n"))}

  #isolation of valid ECM receptors (should be expressed at least by 30% of a population but can be changed via the min.pct parameter)
  if(isTRUE(verbose)){cat(crayon::white("(3/9) evaluating matrisome receptor genes in target population(s)... "))}

  if(!is.null(target)){receiver_celltypes <- target}else{
    receiver_celltypes <- unique(obj@meta.data[,group.column])}

  list_expressed_genes_receiver <- receiver_celltypes %>% unique() %>% lapply(get_expressed_genes, obj, min.pct)
  names(list_expressed_genes_receiver) <- receiver_celltypes
  for(i in names(list_expressed_genes_receiver)){
    s <- obj@meta.data[,group.column]
    names(s) <- rownames(obj@meta.data)
    s <- names(s[s%in%i])
    z <- list_expressed_genes_receiver[i][[1]]
    k <- obj@assays$RNA@counts[rownames(obj@assays$RNA@counts)%in%z,
                               colnames(obj@assays$RNA@counts)%in%s]
    k <- rowMeans(as.matrix(k))
    k <- k[k>=expr.filter]
    list_expressed_genes_receiver[[i]] <- names(k)
  }
  background_expressed_genes <- lapply(list_expressed_genes_receiver,function(x){
    x[x%in%rownames(new_ligand_target_matrix)]
  })
  nnn <- unique(unlist(background_expressed_genes))

  if(length(nnn)<1){
    cat(crayon::red("\nno gene passes the threshold, execution stops \n"))
    cat(crayon::red("\ntry changing the min.pct or the expr.filter parameter \n"))
    stop()
  }

  receptors <- new_lr %>% pull(to) %>% unique()
  expressed_receptors <- lapply(background_expressed_genes,function(x){
    intersect(receptors,x)})
  nnn <- unique(unlist(expressed_receptors))

  if(length(nnn)<1){
    cat(crayon::red("\nno gene passes the threshold, execution stops \n"))
    cat(crayon::red("\ntry changing the min.pct or the expr.filter parameter \n"))
    stop()
  }
  if(isTRUE(verbose)){cat(crayon::green("done \n"))}

  if(isTRUE(force.flushing)){
    gc()
    if(isTRUE(verbose)){cat(crayon::yellow("memory flushed, now restarting \n"))}
  }

  #filtering potentially active ligands
  if(isTRUE(verbose)){cat(crayon::white("(4/9) extracting potential matrisome ligands based on receptors... "))}
  potential_ligands <- new_lr %>% filter(from %in% expressed_ligands & to %in% nnn) %>% pull(from) %>% unique()

  if(length(potential_ligands)<1){
    cat(crayon::red("\nno cell-specific pair was found, execution stops \n"))
    cat(crayon::red("\ntry changing the min.pct or the expr.filter parameter \n"))
    stop()
  }

  lndf <- list()
  for(i in names(expressed_receptors)){
    y <- expressed_receptors[i][[1]]
    nr <- new_lr[new_lr$from%in%potential_ligands &
                   new_lr$to%in%y,]
    if(nrow(nr)<1){next}else{
      nr$target.cell <- i
      lndf[[i]] <- nr
    }
  }
  nnn <- bind_rows(lndf)

  if(nrow(nnn)<1){
    cat(crayon::red("\nno path explains ligand activities in cells, execution stops \n"))
    cat(crayon::red("\ntry changing the min.pct or the expr.filter parameter \n"))
    stop()
  }else{
    potential_ligands <- unique(nnn$from)
    expressed_receptors <- expressed_receptors[unique(nnn$target.cell)]
  }
  if(isTRUE(verbose)){cat(crayon::green("done \n"))}

  #extracting response signatures
  if(!is.null(geneset)){
    if(isTRUE(verbose)){cat(crayon::white("(5/9) applying the response geneset... "))}
    ls_of_interest <- list()
    ls_of_interest[[1]] <- geneset
  }else{
    if(isTRUE(verbose)){cat(crayon::white("(5/9) creating population-specific response geneset(s)..."))}

    ls_of_interest <- list()
    for(i in names(expressed_receptors)){
      z <- expressed_receptors[i][[1]]
      if(length(z)<1){next}else{
        g <- background_expressed_genes[i][[1]]
        g <- g[!g%in%z]
        ln <- table(obj@meta.data[,group.column])
        ln <- ln[names(ln)%in%i]
        if(ln<10){next}else{
          seurat_obj_receiver <- SetIdent(obj,value=ifelse(Idents(obj)%in%i,"1","0"))
          seurat_obj_receiver <- subset(seurat_obj_receiver,features=g)
          DE_table_receiver <- FindMarkers(object = seurat_obj_receiver, ident.1 = "1", ident.2 = "0", min.pct = min.pct, max.cells.per.ident = ln,verbose = F)
          DE_table_receiver$gene <- rownames(DE_table_receiver)
          geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & avg_log2FC >= 0.25) %>% pull(gene)
          geneset_oi <- geneset_oi %>% .[. %in% rownames(new_ligand_target_matrix)]
          lll <- list_expressed_genes_receiver[i][[1]]
          geneset_oi <- geneset_oi[geneset_oi%in%lll]
          ls_of_interest[[i]] <- geneset_oi
        }
      }
    }
  }
  nnn <- unique(unlist(ls_of_interest))

  if(length(nnn)<1){
    cat(crayon::red("no gene for geneset found, execution stops \n"))
    cat(crayon::red("try changing the min.pct or the expr.filter parameter \n"))
    stop()
  }
  if(isTRUE(verbose)){cat(crayon::green("done \n"))}

  #filtering potentially active ligands
  if(isTRUE(verbose)){cat(crayon::white("(6/9) testing potential ligand(s) against response genesets (clipping of ligands might occurr if > max.ligands)... "))}

  ligand_activities <- list()
  for(i in names(ls_of_interest)){
    tr <- ls_of_interest[[i]]
    potential_ligands <- potential_ligands[potential_ligands%in%colnames(new_ligand_target_matrix)]
    if(length(tr)<1){
      next
    }else{
      act <- suppressWarnings(predict_ligand_activities(geneset = ls_of_interest[i][[1]], background_expressed_genes = background_expressed_genes[i][[1]], ligand_target_matrix = new_ligand_target_matrix, potential_ligands = potential_ligands))
      act <- act %>% arrange(-pearson)
      act <- act[act$pearson>0,]
      if(nrow(act)<1){next}else{
        if(nrow(act)>max.ligands){
          act<-act[1:max.ligands,]
        }else{act<-act}
        act$target.cell <- i
        ligand_activities[[i]] <- act}}
  }
  nnn <- bind_rows(ligand_activities)
  nnn <- unique(nnn$test_ligand)
  nnn <- as.character(na.omit(nnn))

  if(length(nnn)<1){
    cat(crayon::red("\nno ligand explains geneset(s), execution stops \n"))
    cat(crayon::red("\ntry changing the min.pct or the expr.filter parameter \n"))
    stop()
  }
  if(isTRUE(verbose)){cat(crayon::green("done \n"))}

  #combining the results
  if(isTRUE(verbose)){cat(crayon::white("(7/9) forming triplets (clipping of targets might occurr if > max.targets)... "))}

  best_upstream_ligands <- as.character(na.omit(unique(bind_rows(ligand_activities)$test_ligand)))

  sends <- list()
  for(i in names(list_expressed_genes_sender)){
    z <- list_expressed_genes_sender[i][[1]]
    z <- z[z%in%best_upstream_ligands]
    if(length(z)<1){next}else{
      sends[[i]] <- data.frame(origin.population=i,ECM.component=z)
    }
  }
  sends <- bind_rows(sends)
  avg_expression_ligands <- suppressWarnings(as.data.frame(t(unlist(AverageExpression(obj, features = best_upstream_ligands)[[1]]))))
  nm <- rownames(avg_expression_ligands)
  if(length(avg_expression_ligands)>1){
    avg_expression_ligands <- suppressMessages(reshape2::melt(avg_expression_ligands))
    avg_expression_ligands$origin.population <- nm
    names(avg_expression_ligands)[1] <- "ECM.component"
  }else{
    avg_expression_ligands <- data.frame(ECM.component=best_upstream_ligands,
                                         value=avg_expression_ligands,
                                         origin.population=nm)
    names(avg_expression_ligands)[2] <- "value"
  }
  sends <- distinct(merge(sends,avg_expression_ligands,by=c("origin.population","ECM.component")))
  sends <- subset(sends,sends$value>expr.filter)

  if(nrow(sends)<1){
    cat(crayon::red("\nno population expresses any ligand above 1, execution stops \n"))
    cat(crayon::red("\ntry changing the min.pct or the expr.filter parameter \n"))
    stop()
  }else{
    names(sends)[3] <- "ECM.value"}

  lndf <- bind_rows(lndf)
  k <- list()
  for(i in names(ligand_activities)){
    z <- ligand_activities[i][[1]]
    r <- subset(lndf,lndf$target.cell==i)
    r <- subset(r,r$from%in%unique(z$test_ligand))
    if(nrow(r)<1){next}else{
      avg_expression_r <- suppressWarnings(as.data.frame(t(unlist(AverageExpression(obj, features = unique(r$to))[[1]]))))
      avg_expression_r <- avg_expression_r[rownames(avg_expression_r)%in%i,]
      if(length(avg_expression_r)>1){
        nm <- rownames(avg_expression_r)
        avg_expression_r <- suppressMessages(reshape2::melt(avg_expression_r))
        avg_expression_r$receiving.population <- nm
        names(avg_expression_r)[1] <- "to"
      }else{
        avg_expression_r <- data.frame(to=r$to,
                                       value=avg_expression_r,
                                       receiving.population=i)

      }
      r <- r[,c(1,2,5)]
      names(r)[3] <- "receiving.population"
      r <- distinct(merge(r,avg_expression_r,by=c("to","receiving.population")))
      r <- subset(r,r$value>expr.filter)
      if(nrow(r)<1){next}else{
        r <- r[,c(3,1,2,4)]
        k[[i]] <- r}
    }}
  k <- bind_rows(k)
  names(k)[1:2] <- c("ECM.component","ECM.receptor")
  names(k)[4] <- "receptor.value"

  if(nrow(k)<1){
    cat(crayon::red("\nno population expresses any receptor above 1, execution stops \n"))
    cat(crayon::red("\ntry changing the min.pct or the expr.filter parameter \n"))
    stop()
  }

  k2 <- list()
  for(i in names(ligand_activities)){
    z <- ls_of_interest[i][[1]]
    avg_expression_t <- suppressWarnings(as.data.frame(t(unlist(AverageExpression(obj, features = rownames(z))[[1]]))))
    avg_expression_t <- avg_expression_t[rownames(avg_expression_t)%in%i,]
    if(length(avg_expression_t)>1){
      avg_expression_t <- suppressMessages(reshape2::melt(avg_expression_t))
      avg_expression_t$receiving.population <- i
      names(avg_expression_t)[1] <- "target.gene"
    }else{
      avg_expression_t <- data.frame(target.gene=rownames(z),
                                     value=avg_expression_t,
                                     receiving.population=i)

    }
    df <- subset(avg_expression_t,avg_expression_t$value>expr.filter)
    if(nrow(df)<1){next}else{
      df <- df[order(-df$value),]
      if(nrow(df)>max.targets){df<-df[1:max.targets,]}else{df<-df}
      k2[[i]] <- df}
  }
  k2 <- bind_rows(k2)
  names(k2)[2] <- "target.value"

  if(nrow(k2)<1){
    cat(crayon::red("\nno population expresses any target above 1, execution stops \n"))
    cat(crayon::red("\ntry changing the min.pct or the expr.filter parameter \n"))
    stop()
  }

  fin <- distinct(merge(sends,k,by="ECM.component"))
  fin <- fin[,c(2,1,3,5,4,6)]
  fin <- distinct(merge(fin,k2,by="receiving.population"))
  fin <- fin[,c(2:4,1,5,6:8)]
  if(isTRUE(verbose)){cat(crayon::green("done \n"))}

  #explainable paths
  ### quick differential correlation - VERSION 2

  if(!isFALSE(explainable)){
    if(isTRUE(verbose)){cat(crayon::white(paste0("(8/9) finding explainable paths... ")))}

    om <- as.matrix(obj@assays$RNA@counts)
    om[is.na(om)] <- 0
    rands <- list()
    for(i in 1:10){
      rm <- om[sample(rownames(om),100,replace = T),
               sample(colnames(om),100,replace = T)]
      cf <- mean(rowMeans(suppressWarnings(cor(rm,method = "spearman"))))
      cf <- ifelse(is.na(cf),0,cf)
      rands[[i]] <- cf
    }

    rands <- mean(unlist(rands))

    rp <- as.character(unique(fin$receiving.population))
    l1 <- list()
    for(i in rp){
      z <- subset(fin,fin$receiving.population==i)
      g <- as.character(unique(z$target.gene))
      rs <- as.character(unique(z$ECM.receptor))
      s <- obj@meta.data[obj@meta.data[,group.column]==i,]
      s <- unique(rownames(s))
      if(length(s)<1){next}else{
        #s2 <- obj@meta.data[obj@meta.data[,group.column]!=i,]
        #s2 <- unique(rownames(s2))
        tarrec.m <- om[rownames(om)%in%unique(c(g,rs)),
                       colnames(om)%in%s]
        mat1 <- suppressWarnings(cor(t(tarrec.m),method="spearman"))
        mat1[upper.tri(mat1)] <- NA
        rw <- rownames(mat1)
        cl <- colnames(mat1)
        mat1 <- as.data.frame(cbind(which(!is.na(mat1),arr.ind = TRUE),na.omit(as.vector(mat1))))
        mat1 <- mat1[mat1$V3>rands,]
        mat1 <- mat1[mat1$row!=mat1$col,]
        if(nrow(mat1)<1){next}else{
          mat1$row <- rw[mat1$row]
          mat1$col <- cl[mat1$col]
          names(mat1) <- c("ECM.receptor","target.gene")
          tf <- new_gr[new_gr$to%in%mat1$target.gene,]
          rs <- as.character(unique(tf$from))
          #s <- obj@meta.data[obj@meta.data[,group.column]==i,]
          #s <- unique(rownames(s))
          #s2 <- obj@meta.data[obj@meta.data[,group.column]!=i,]
          #s2 <- unique(rownames(s2))
          tartf.m <- om[rownames(om)%in%unique(c(mat1$target.gene,rs)),
                        colnames(om)%in%s]
          mat2 <- suppressWarnings(cor(t(tartf.m),method="spearman"))
          mat2[upper.tri(mat2)] <- NA
          rw <- rownames(mat2)
          cl <- colnames(mat2)
          mat2 <- as.data.frame(cbind(which(!is.na(mat2),arr.ind = TRUE),na.omit(as.vector(mat2))))
          mat2 <- mat2[mat2$V3>rands,]
          mat2 <- mat2[mat2$row!=mat2$col,]
          if(nrow(mat2)<1){next}else{
            mat2$row <- rw[mat2$row]
            mat2$col <- cl[mat2$col]
            names(mat2) <- c("TF","target.gene")
            fifin <- distinct(merge(mat1,mat2,by="target.gene"))
            if(nrow(fifin)<1){next}else{
              nobj <- om[rownames(om)%in%unique(fifin$TF),
                         colnames(om)%in%s]
              df <- data.frame(TF=rownames(nobj),TF.value=rowMeans(nobj))
              df$receiving.population <- i
              df <- merge(fifin,df,by=c("TF"))
              df <- df[,c(1,2,3,6,7)]
              l1[[i]] <- df
            }
          }
        }
      }}
    l1 <- bind_rows(l1)
    rownames(l1) <- NULL

    new_fin <- distinct(merge(fin,l1,by=c("receiving.population","target.gene","ECM.receptor")))
    new_fin <- new_fin[,c(4:6,1,3,7,2,8,9,10)]
  }else{
    if(isTRUE(verbose)){cat(crayon::white("(8/9) not searching for explainable paths \n"))}
    new_fin <- fin
    new_fin$TF <- "not searched"
    new_fin$TF.value <- 1
  }

  if(isTRUE(add.RF)){
    rlst <- list()
    for(i in unique(new_fin$receiving.population)){
      z <- subset(new_fin,new_fin$receiving.population==i)
      z <- unique(c(z$ECM.receptor,z$TF))
      z <- z[z!="not searched"]
      if(length(z)<1){next}else{
        mk <- obj@assays$RNA@counts[rownames(obj@assays$RNA@counts)%in%z,]
        mk <- as.matrix(mk)
        Y <- rownames(obj@meta.data[obj@meta.data[,group.column]==i,])
        Y <- ifelse(colnames(mk)%in%Y,"first","second")
        Y <- factor(Y,levels=c("first","second"))
        mk[is.na(mk)] <- 0
        mk <- t(mk)
        mk <- as.data.frame(apply(mk,2,scale))
        mk$identity <- Y

        set.seed(12345)
        sds <- vector(mode = "list", length = 5)
        inTrain <- createDataPartition(mk$identity, p = 0.7, list = FALSE)
        trainData <- mk[inTrain,]
        testData <- mk[-inTrain,]

        suppressMessages(suppressWarnings(
          my_control <- trainControl(
            method="boost",
            number=5,
            savePredictions="final",
            classProbs=TRUE,
            index=createResample(trainData$identity,5),
            summaryFunction=twoClassSummary
          )
        ))

        suppressMessages(suppressWarnings(
          rf_default <- train(identity~.,
                              data=mk,
                              method='rf',
                              metric='Accuracy',
                              trControl=my_control
          )
        ))

        vs <- vs <- as.data.frame(varImp(rf_default)$importance)
        vs$gene <- rownames(vs)
        vs <- vs[order(-vs$Overall),]
        if(nrow(vs)>max.ligands){
          vs <- vs[1:max.ligands,]
        }else{vs <- vs}

        fincut <- subset(new_fin,new_fin$receiving.population==i)
        fincut <- fincut[fincut$ECM.receptor%in%vs$gene & fincut$TF%in%vs$gene,]

        if(nrow(fincut)<1){
          fincut <- subset(new_fin,new_fin$receiving.population==i)
          fincut <- fincut[fincut$ECM.receptor%in%vs$gene | fincut$TF%in%vs$gene,]

          if(nrow(fincut)<1){

        }else{
          rlst[[i]] <- fincut}
        }

      }
    }
    rlst <- bind_rows(rlst)
    new_fin <- rlst
  }

  if(nrow(new_fin)<1){
    cat(crayon::red("\nno ligand path was found, execution stops \n"))
    cat(crayon::red("\ntry changing the min.pct or the expr.filter parameter or to set add.RF to FALSE\n"))
    stop()
  }
  if(isTRUE(verbose)){cat(crayon::green("done \n"))}

  #export the results

  if(isTRUE(scale.results)){
    if(isTRUE(verbose)){cat(crayon::white("(9/9) scaling the results to [0,1] and annotating... "))}
    library(scales)
    df_scaled <- new_fin %>% group_by(origin.population) %>% mutate(ECM.value=round(rescale(ECM.value,to=c(0.01,1)),3)) %>%
      group_by(receiving.population) %>% mutate(receptor.value=round(rescale(receptor.value,to=c(0.01,1)),3)) %>%
      group_by(receiving.population) %>% mutate(target.value=round(rescale(target.value,to=c(0.01,1)),3)) %>%
      group_by(receiving.population) %>% mutate(TF.value=round(rescale(TF.value,to=c(0.01,1)),3))
    fin <- data.frame(df_scaled)

    annofin <- list()
    for(i in unique(fin$receiving.population)){
      z <- fin[fin$receiving.population==i,]
      z <- unique(z$ECM.receptor)
      l <- list()
      for(w in ncpl$complex){
        s <- subset(ncpl,ncpl$complex==w)
        if(isFALSE(use.ensembl)){
          s <- unlist(strsplit(s$genes,split=";"))}else{
            s <- unlist(strsplit(s$genes,split=";"))
            s <- suppressMessages(mapIds(org.Hs.eg.db, keys=s, column="ENSEMBL", keytype="SYMBOL", multiVals="first"))
          }
        y <- length(intersect(s,z))
        if(y<2){next}else{
          y <- intersect(s,z)
          l[[w]] <- data.frame(receiving.population=i,ECM.receptor=y,annotation=paste0(w," is enriched in the population"))}
      }
      l <- bind_rows(l)
      annofin[[i]] <- bind_rows(l)
    }
    annofin <- bind_rows(annofin)
    if(nrow(annofin)>0){
      fin <- distinct(merge(fin,annofin,by=c("receiving.population","ECM.receptor"),all=T))
      fin[is.na(fin)] <- ""
      fin <- fin[,c(3:5,1,2,6:ncol(fin))]}else{
        fin$annotation <- ""
      }

    attributes(fin)$mode <- "matricom.scaled"}else{
      if(isTRUE(verbose)){cat(crayon::white("(9/9) annotating results... "))}
      fin<-new_fin

      annofin <- list()
      for(i in unique(fin$receiving.population)){
        z <- fin[fin$receiving.population==i,]
        z <- unique(z$ECM.receptor)
        l <- list()
        for(w in ncpl$complex){
          s <- subset(ncpl,ncpl$complex==w)
          s <- unlist(strsplit(s$genes,split=";"))
          y <- length(intersect(s,z))
          if(y<2){next}else{
            y <- intersect(s,z)
            l[[w]] <- data.frame(receiving.population=i,ECM.receptor=y,annotation=paste0(w," is enriched in the receiver population"))}
        }
        l <- bind_rows(l)
        annofin[[i]] <- bind_rows(l)
      }
      annofin <- bind_rows(annofin)
      if(nrow(annofin)>0){
        fin <- distinct(merge(fin,annofin,by=c("receiving.population","ECM.receptor"),all=T))
        fin[is.na(fin)] <- ""
        fin <- fin[,c(3:5,1,2,6:ncol(fin))]}else{
          fin$annotation <- ""
        }
      attributes(fin)$mode <- "matricom.raw"
    }


  # did we have Ensembl IDs in the output? Convert to gene names.
  if(isTRUE(use.ensembl)){
    fin$ECM.component <- as.character(mapIds(org.Hs.eg.db, keys=fin$ECM.component, column="SYMBOL", keytype="ENSEMBL"))
    fin$ECM.receptor <- as.character(mapIds(org.Hs.eg.db, keys=fin$ECM.receptor, column="SYMBOL", keytype="ENSEMBL"))
    fin$target.gene <- as.character(fin$target.gene)
    fin$target.gene <- as.character(mapIds(org.Hs.eg.db, keys=fin$target.gene, column="SYMBOL", keytype="ENSEMBL"))
    fin$TF <- as.character(mapIds(org.Hs.eg.db, keys=fin$TF, column="SYMBOL", keytype="ENSEMBL"))
  }

  if(isTRUE(verbose)){cat(crayon::green("done \n "))}

  return(fin)
}

#' flowgraph: representing cell-ECM communications with alluvials

#' @param result.obj character. The results from a previous call to matricom.
#' @param targets logical. Whether to add target genes (and their connections) to the graph. Since it might clutter the graph, default is FALSE.
#' @param color.by characters. Color the ribbons (alluvials) by the sender or the receiver population, or by the ECM component. Default is "sender".
#' @param simplify logical. Whether to cut all triplets where a component's expression (scaled) is below 0.1. Default is TRUE, to de-clutter the graph. A failsafe mechanism is triggered if all values fall below 0.1 and this option is TRUE.

#' @return an alluvial plot made with ggplot2
#' @export

#' @examples flowgraph(results)

flowgraph <- function(result.obj=NULL,
                      targets=TRUE,
                      color.by=c("sender", "ecm", "receiver"),
                      simplify=TRUE){

  #error handling
  if(is.null(result.obj)){
    cat(crayon::red("no data provided, execution stops \n"))
    stop()
  }
  if(class(result.obj)!="data.frame"){
    cat(crayon::red("data are not in a valid format (data.frame), execution stops \n"))
    stop()
  }

  if(attributes(result.obj)$mode=="matricom.scaled"){
    mode <- "scaled"}else{mode <- "raw"}

  obj <- result.obj
  library(ggalluvial)

  if(isFALSE(targets)){
    obj$target.gene <- NULL
    obj$target.value <- NULL
    obj <- distinct(obj)

    if(mode=="scaled"){obj <- obj}else{
      obj <- obj %>% group_by(origin.population) %>% mutate(ECM.value=round(rescale(ECM.value,to=c(0.01,1)),3)) %>%
        group_by(receiving.population) %>% mutate(receptor.value=round(rescale(receptor.value,to=c(0.01,1)),3)) %>%
        group_by(receiving.population) %>% mutate(target.value=round(rescale(target.value,to=c(0.01,1)),3)) %>%
        group_by(receiving.population) %>% mutate(TF.value=round(rescale(TF.value,to=c(0.01,1)),3))
      }

    if(isTRUE(simplify)){
      p <- obj[obj$ECM.value>=0.1 & obj$receptor.value>=0.1,]
      if(nrow(p)>=1){
        obj<-p
        tit <- "simplified graph (values are scaled by population and only those >= 0.1 are reported)"
      }else{
        obj <- obj[obj$ECM.value>=0.01 & obj$receptor.value>=0.01,]
        tit <- "simplified graph (values are scaled by population and only those >= 0.01 are reported) \n(failsafe mechanism triggered))"
      }
      }else{
      tit <- "full graph (values are scaled by population and all are reported)"
    }
    p <- NULL
    g <- as.data.frame(table(obj$origin.population,obj$ECM.component,obj$receiving.population))
    g <- g[g$Freq>0,]

    cl <- match.arg(NULL,color.by)
    s <- cl
    if(cl=="sender"){cl<-g$Var1}else{
      if(cl=="ecm"){cl<-g$Var2}else{
        cl<-g$Var3
      }
    }

    p <- ggplot(g,
                aes(y = Freq, axis1 = Var1, axis2 = Var2, axis3 = Var3)) +
      geom_alluvium(aes(fill = cl), color="grey" ,width = 1/12) +
      geom_stratum(width = 1/12, fill = "grey80", color = "black") +
      geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
      scale_x_discrete(limits = c("original.population","ECM.component","receiver.population"),
                       expand = c(.05, .05)) +
      theme_bw() + theme(axis.ticks.y = element_blank(),
                         axis.text.y = element_blank()) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      ylab("") + ggtitle(tit)

    if(s=="sender"){p<-p+scale_fill_discrete(name = "population originating \nthe ECM component")}else{
      if(s=="ecm"){p<-p+scale_fill_discrete(name = "ECM component")}else{
        p<-p+scale_fill_discrete(name = "population receiving \nthe ECM component")
      }}

    suppressWarnings(print(p))}else{

      if(mode=="scaled"){obj <- obj}else{
        obj <- obj %>% group_by(origin.population) %>% mutate(ECM.value=round(rescale(ECM.value,to=c(0.01,1)),3)) %>%
          group_by(receiving.population) %>% mutate(receptor.value=round(rescale(receptor.value,to=c(0.01,1)),3)) %>%
          group_by(receiving.population) %>% mutate(target.value=round(rescale(target.value,to=c(0.01,1)),3)) %>%
          group_by(receiving.population) %>% mutate(TF.value=round(rescale(TF.value,to=c(0.01,1)),3))
        }

      if(isTRUE(simplify)){
        p <- obj[obj$ECM.value>=0.1 & obj$receptor.value>=0.1 & obj$TF.value>=0.1,]
        if(nrow(p)>=1){
          obj<-p
          tit <- "simplified graph (values are scaled by population and only those >= 0.1 are reported)"
        }else{
          obj <- obj[obj$ECM.value>=0.01 & obj$receptor.value>=0.01 & obj$TF.value>=0.01,]
          tit <- "simplified graph (values are scaled by population and only those >= 0.01 are reported) \n(failsafe mechanism triggered))"
        }
        }else{
        tit <- "full graph (values are scaled by population and all are reported)"
        }
      p <- NULL
      g <- as.data.frame(table(obj$origin.population,obj$ECM.component,obj$target.gene,obj$receiving.population))
      g <- g[g$Freq>0,]

      cl <- match.arg(NULL,color.by)
      s <- cl
      if(cl=="sender"){cl<-g$Var1}else{
        if(cl=="ecm"){cl<-g$Var2}else{
          cl<-g$Var4
        }
      }

      p <- ggplot(g,
                  aes(y = Freq, axis1 = Var1, axis2 = Var2, axis3 = Var4, axis4 = Var3)) +
        geom_alluvium(aes(fill = as.factor(cl)), color="grey" ,width = 1/12) +
        geom_stratum(width = 1/12, fill = "grey80", color = "black") +
        geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
        scale_x_discrete(limits = c("original.population","ECM.component","receiver.population","target"),
                         expand = c(.05, .05)) +
        theme_bw() + theme(axis.ticks.y = element_blank(),
                           axis.text.y = element_blank()) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        ylab("") + ggtitle(tit)

      if(s=="sender"){p<-p+scale_fill_discrete(name = "population originating \nthe ECM component")}else{
        if(s=="ecm"){p<-p+scale_fill_discrete(name = "ECM component")}else{
          p<-p+scale_fill_discrete(name = "population receiving \nthe ECM component")
        }}

      suppressWarnings(print(p))
    }
}

#' annograph: matrisome annotations with nested pie charts

#' @param result.obj character. The results from a previous call to matricom.
#' @param select character. Annotate matrisome families and categories for the sender or the receiver population. Default is "sender".

#' @return a PieDonut (a nested pie chart) made with ggplot2
#' @export

#' @examples annograph(results)

annograph <- function(result.obj=NULL,
                      select=c("sender","receiver")){

  #error handling
  if(is.null(result.obj)){
    cat(crayon::red("no data provided, execution stops \n"))
    stop()
  }
  if(class(result.obj)!="data.frame"){
    cat(crayon::red("data are not in a valid format (data.frame), execution stops \n"))
    stop()
  }

  obj <- result.obj
  data(mat, package = "matRicom")

  cl <- match.arg(NULL,select)

  if(cl=="sender"){
    obj <- distinct(obj[,c(1,2)])
    g <- as.data.frame(table(obj$origin.population,obj$ECM.component))
    g <- g[g$Freq>0,]

    g1 <- distinct(merge(g,mat,by.x="Var2",by.y="gene",all.x=T))
    g1[is.na(g1)] <- "not.matrisome"
    g1 <- aggregate(g1$Freq,list(g1$Var1,g1$family),sum)
    names(g1) <- c("sender.population","ecm.component","value")

    suppressWarnings(PieDonut(g1,aes(sender.population,ecm.component,count=value)))
  }else{
    obj <- distinct(obj[,c(4,2)])
    g <- as.data.frame(table(obj$receiving.population,obj$ECM.component))
    g <- g[g$Freq>0,]

    g2 <- distinct(merge(g,mat,by.x="Var2",by.y="gene",all.x=T))
    g2[is.na(g2)] <- "not.matrisome"
    g2 <- aggregate(g2$Freq,list(g2$Var1,g2$family),sum)
    names(g2) <- c("receiver.population","ecm.component","value")

    suppressWarnings(PieDonut(g2,aes(receiver.population,ecm.component,count=value)))

  }
}

#' compgraph: representing differential cell-ECM communications with alluvials

#' @param result.obj character. A list with strictly two results, each scaled, from previous calls to matricom.
#' @param targets logical. Whether to add target genes (and their connections) to the graph. Since it might clutter the graph, default is FALSE.
#' @param only.unique logical. Whether to plot only the communications (alluvials) that differ between conditions (present in one and absent in another). Default is TRUE.

#' @return an alluvial plot made with ggplot2
#' @export

#' @examples compgraph(list(a=results_1,b=results_2))

compgraph <- function(result.obj=NULL,
                      targets=FALSE,
                      only.unique=TRUE){

  #error handling
  if(is.null(result.obj)){
    cat(crayon::red("no data provided, execution stops \n"))
    stop()
  }
  if(class(result.obj)!="list"){
    cat(crayon::red("data are not in a valid format (list), execution stops \n"))
    stop()
  }
  if(length(result.obj)>2){
    cat(crayon::red("only two result files can be compared at once, execution stops \n"))
    stop()
  }
  if(attributes(result.obj[[1]])$mode!="matricom.scaled"){
    cat(crayon::red("only scaled results are accepted, execution stops \n"))
    stop()
  }
  if(attributes(result.obj[[2]])$mode!="matricom.scaled"){
    cat(crayon::red("only scaled results are accepted, execution stops \n"))
    stop()
  }
  if(length(names(result.obj))<2){
    cat(crayon::red("data are not a named list, execution stops \n"))
    stop()
  }

  obj <- result.obj
  library(ggalluvial)

  if(isFALSE(targets)){
    g <- lapply(obj, function(x){
      x$target.gene <- NULL
      x$target.value <- NULL
      x <- distinct(x)
      x <- x[x$ECM.value>0.01 & x$receptor.value>0.01,]
      x <- as.data.frame(table(x$ECM.component,x$ECM.receptor))
      x <- x[x$Freq>0,]
    })

    for(i in names(g)){
      g[i][[1]]$condition <- i
      g[i][[1]]$edge <- paste(g[i][[1]]$Var1,g[i][[1]]$Var2,sep="_")
    }

    ints <- unique(c(intersect(g[[1]]$edge,g[[2]]$edge),
                     intersect(g[[2]]$edge,g[[1]]$edge)))

    for(i in names(g)){
      g[i][[1]]$evaluation <- ifelse(g[i][[1]]$edge%in%ints,"common","unique")
    }

    g <- bind_rows(g)

    if(isTRUE(only.unique)){

      g <- subset(g,g$evaluation=="unique")
      g$evaluation <- paste0(g$evaluation,"_",g$condition)
      p <- ggplot(g,
                  aes(y = Freq, axis1 = Var1, axis2 = Var2)) +
        geom_alluvium(aes(fill = evaluation), color="grey" ,width = 1/12) +
        geom_stratum(width = 1/12, fill = "grey80", color = "black") +
        geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
        scale_x_discrete(limits = c("ECM.component","ECM.receptor"),
                         expand = c(.05, .05)) +
        theme_bw() + theme(axis.ticks.y = element_blank(),
                           axis.text.y = element_blank()) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        ylab("") + ggtitle("unique signals")

      suppressWarnings(print(p))}else{

        g$evaluation <- ifelse(g$evaluation=="common","common",paste0(g$evaluation,"_",g$condition))
        p <- ggplot(g,
                    aes(y = Freq, axis1 = Var1, axis2 = Var2)) +
          geom_alluvium(aes(fill = evaluation), color="grey" ,width = 1/12) +
          geom_stratum(width = 1/12, fill = "grey80", color = "black") +
          geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
          scale_x_discrete(limits = c("ECM.component","ECM.receptor"),
                           expand = c(.05, .05)) +
          theme_bw() + theme(axis.ticks.y = element_blank(),
                             axis.text.y = element_blank()) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          ylab("") + ggtitle("unique and common signals")

        suppressWarnings(print(p))}
  }else{
    g <- lapply(obj, function(x){
      x <- x[x$ECM.value>0.01 & x$receptor.value>0.01 & x$target.value>0.01,]
      x <- as.data.frame(table(x$ECM.component,x$ECM.receptor,x$target.gene))
      x <- x[x$Freq>0,]
    })

    for(i in names(g)){
      g[i][[1]]$condition <- i
      g[i][[1]]$edge <- paste(g[i][[1]]$Var1,g[i][[1]]$Var2,g[i][[1]]$Var3,sep="_")
    }

    ints <- unique(c(intersect(g[[1]]$edge,g[[2]]$edge),
                     intersect(g[[2]]$edge,g[[1]]$edge)))

    for(i in names(g)){
      g[i][[1]]$evaluation <- ifelse(g[i][[1]]$edge%in%ints,"common","unique")
    }

    g <- bind_rows(g)

    if(isTRUE(only.unique)){

      g <- subset(g,g$evaluation=="unique")
      g$evaluation <- paste0(g$evaluation,"_",g$condition)
      p <- ggplot(g,
                  aes(y = Freq, axis1 = Var1, axis2 = Var2, axis3 = Var3)) +
        geom_alluvium(aes(fill = evaluation), color="grey" ,width = 1/12) +
        geom_stratum(width = 1/12, fill = "grey80", color = "black") +
        geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
        scale_x_discrete(limits = c("ECM.component","ECM.receptor","target.gene"),
                         expand = c(.05, .05)) +
        theme_bw() + theme(axis.ticks.y = element_blank(),
                           axis.text.y = element_blank()) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        ylab("") + ggtitle("unique signals")

      suppressWarnings(print(p))}else{

        g$evaluation <- ifelse(g$evaluation=="common","common",paste0(g$evaluation,"_",g$condition))
        p <- ggplot(g,
                    aes(y = Freq, axis1 = Var1, axis2 = Var2, axis3 = Var3)) +
          geom_alluvium(aes(fill = evaluation), color="grey" ,width = 1/12) +
          geom_stratum(width = 1/12, fill = "grey80", color = "black") +
          geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
          scale_x_discrete(limits = c("ECM.component","ECM.receptor","target.gene"),
                           expand = c(.05, .05)) +
          theme_bw() + theme(axis.ticks.y = element_blank(),
                             axis.text.y = element_blank()) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          ylab("") + ggtitle("unique and common signals")

        suppressWarnings(print(p))}
  }
}

#' netgraph: representing cell-ECM communications with networks, focusing on intracellular paths

#' @param result.obj character. The results from a previous call to matricom.
#' @param receiving.population character. The name of one receiving population in the results. If NULL (default), a random population is picked.
#' @param ECM.component character. The name of one ECM component impacting on a receiving population in the results. If NULL (default), a random ECM component is picked.
#' @param label.size numeric. Size of the vertex labels in the network plot. If NULL (default), no labels are added to the nodes of the plot.
#' @param simplified logical. Whether to remove ECM-loops (an ECM component might also work as a receptor and/or be a target of its own stimulation). Since it helps de-cluttering the graph, default is TRUE.
#' @param interactive logical. Whether to use the VisNetwork package to plut an interactive form of the graph. Default is FALSE.

#' @return a network plot made with igraph
#' @export

#' @examples netgraph(results)

netgraph <- function(result.obj=NULL,
                     receiving.population=NULL,
                     ECM.component=NULL,
                     label.size=NULL,
                     simplified=TRUE,
                     interactive=FALSE){

  #error handling
  if(is.null(result.obj)){
    cat(crayon::red("no data provided, execution stops \n"))
    stop()
  }
  if(is.null(receiving.population)){
    cat(crayon::white("no population chosen, picking a random one \n"))
  }
  if(is.null(ECM.component)){
    cat(crayon::white("no ECM component chosen, picking a random one \n"))
  }
  require(igraph)
  require(visNetwork)

  obj <- result.obj

  if(is.null(receiving.population)){
    receiving.population <- sample(unique(obj$receiving.population),1)
  }else{
    receiving.population <- receiving.population
  }

  z <- subset(obj,obj$receiving.population==receiving.population)

  if(is.null(ECM.component)){
    ECM.component <- sample(unique(z$ECM.component),1)
  }else{
    ECM.component <-ECM.component
  }

  z <- z[z$ECM.component%in%ECM.component,]
  z <- distinct(z[,c(1,2,4,5,7,9)])

  b <- z[,c(2,3)]
  n2 <- graph.data.frame(b)
  n2 <- simplify(n2)

  d <- z[,c(3,4)]
  n3 <- graph.data.frame(d)
  n3 <- simplify(n3)

  e <- z[,c(4,6)]
  if(isTRUE(simplified)){
  e <- e[!(e$ECM.receptor %in% unique(b$ECM.component)),]}
  n4 <- graph.data.frame(e)
  n4 <- simplify(n4)

  g <- z[,c(6,5)]
  if(isTRUE(simplified)){
  g <- g[!(g$target.gene %in% unique(b$ECM.component)),]}
  n5 <- graph.data.frame(g)
  n5 <- simplify(n5)

  fin <- n2 %u% n3 %u% n4 %u% n5
  fin <- simplify(fin)

  V(fin)$color <- ifelse(
    names(V(fin))%in%unique(b$ECM.component),"red",
    ifelse(names(V(fin))%in%unique(d$ECM.receptor),"orange",
    ifelse(names(V(fin))%in%unique(b$receiving.population),"yellow",
    ifelse(names(V(fin))%in%unique(g$TF),"pink","green"))))

  if(isFALSE(interactive)){

  if(is.null(label.size)){
  plot.igraph(fin,
              layout=layout_with_sugiyama(fin,
                                          attributes="all",
                                          hgap=100, vgap=100),
              edge.arrow.size=0,edge.curved=F,
              edge.width=2,
              vertex.label=NA)}else{
  lb <- label.size
  plot.igraph(fin,
              layout=layout_with_sugiyama(fin,
                                          attributes="all",
                                          hgap=100, vgap=100),
              edge.arrow.size=0,edge.curved=F,
              edge.width=2,
              vertex.label.cex=lb)
              }
  title(paste0("paths of ",ECM.component, " in ",receiving.population),cex.main=1,col.main="black")
  }else{
    data <- toVisNetworkData(fin)
    visNetwork(data$nodes, data$edges,
               main=paste0("paths of ",ECM.component, " in ",receiving.population)) %>%
      visEdges(arrows = "to") %>%
      visNodes(size = 10) %>%
      visOptions(highlightNearest = list(enabled = T, hover = T),
                 nodesIdSelection = T) %>%
      visInteraction(navigationButtons = TRUE)
  }
}

#' matricom.spatial: cell-ECM communications in Spatial RNAseq data

#' @param seurat.obj character. A spatial Seurat object.
#' @param group.column character. Name of the column in metadata with cell identities to be used. Strongly suggested to use cluster identities.
#' @param min.pct numeric. Minimum % of cells with non-zero expression of a ligand/receptor/target within a cluster. All ligands/receptors/targets less present will be ignored. Only for ligands, a total of min.pct is also allowed from the total, on top of cell-specific quantities. Default is 0.5 (50%).
#' @param expr.filter numeric. Mean expression of a ligand/receptor/target within a cluster. All ligands/receptors/targets less expressed will be ignored. Default is 1.
#' @param target character. A receiver cell type to focus on. If selected, only communications affecting that cell type are considered. Default is NULL, all cell types are considered.
#' @param geneset character. A response geneset (a vector of gene names matching those of the counts in the object) to test ligands against. If NULL (default), a recursive cell type-vs.-all DEG search is launched using Seurat functions.
#' @param max.ligands numeric. Maximum number of ligands to test for response. Default is 10, in decreasing Pearson order.
#' @param max.targets numeric. Maximum number of targets to continue with after having tested the response. Default is 10, in decreasing gene expression order.
#' @param explainable logical. Whether to apply differential correlation analysis to targets' regulators. In the analysis, all known regulators of each target gene per receiving population are correlated to their targets with Spearman correlation. In parallel, the same tests are repeated across the rest of sample (non-receiving population). Results are confronted and the regulator with the highest differential (positive) correlation in the receiving population is chosen as the master regulator and reported. Default is TRUE.
#' @param add.RF logical. Whether to add a random forest run to the results from the correlation analysis. Features are ordered by their VarImp and as many as max.ligands (up to twice that value, once for receptors and once for TF) are chosen. Default is TRUE.
#' @param scale.results logical. Whether to rescale the gene expression values of ligands/receptors/targets to the [0.01-1] interval by cell type to ease cross-comparisons. Default is TRUE.
#' @param verbose logical. Progress messages are printed to the console. Default is TRUE.

#' @return a dataframe with ligand/receptor/TF/target quadruplets, the sender and receiver cells for each, and the mean expression values of each triplet.
#' @export

#' @examples matricom.spatial(obj,"group.column")

matricom.spatial <- function(seurat.obj=NULL,
                             group.column=NULL,
                             min.pct=0.50,
                             expr.filter=1,
                             target=NULL,
                             geneset=NULL,
                             max.ligands=10,
                             max.targets=10,
                             explainable=TRUE,
                             add.RF=TRUE,
                             scale.results=TRUE,
                             verbose=T){

  #error handling
  if(is.null(seurat.obj)){
    cat(crayon::red("no data provided, execution stops \n"))
    stop()
  }
  if(class(seurat.obj)!="Seurat"){
    cat(crayon::red("data are not a valid Seurat object, execution stops \n"))
    stop()
  }
  if(nrow(seurat.obj@meta.data)<1){
    cat(crayon::red("data are not a valid Seurat object, execution stops \n"))
    stop()
  }
  #different scaling techniques need to be taken into consideration!
  if(nrow(seurat.obj@assays$SCT@scale.data)<1){
    cat(crayon::red("data must be scaled, execution stops \n"))
    stop()
  }
  if(length(intersect(names(seurat.obj@assays),"SCT"))<1){
    cat(crayon::red("spatial data must be SCT-scaled, execution stops \n"))
    stop()
  }
  if(is.null(group.column)){
    cat(crayon::red("no group assignment column provided, execution stops \n"))
    stop()
  }
  if(length(intersect(colnames(seurat.obj@meta.data),group.column))<1){
    cat(crayon::red("group assignment column missing in object metadata slot, execution stops \n"))
    stop()
  }

  #read in the matricom objects
  if(isTRUE(verbose)){cat(crayon::white("(1/9) loading matRicom objects... "))}
  data(matricom_obj, package = "matRicom")
  data(ncpl, package = "matRicom")
  new_lr <- matricom_obj$new_lr
  new_gr <- matricom_obj$new_gr
  new_ligand_target_matrix <- matricom_obj$new_ligand_target_matrix
  if(isTRUE(verbose)){cat(crayon::green("done \n"))}

  #isolation of valid ECM genes (should be expressed at least by 30% of a population but can be changed via the min.pct parameter)
  if(isTRUE(verbose)){cat(crayon::white("(2/9) evaluating active matrisome genes in cells... "))}

  obj <- seurat.obj

  Idents(obj) <- obj@meta.data[,group.column]

  if(length(intersect(colnames(new_ligand_target_matrix),
                      rownames(obj@assays$SCT@counts)))<1){
    cat(crayon::red("\nno gene in data match with possible ligands, execution stops \n"))
    cat(crayon::red("\ntry checking the format of genes in your data \n"))
    stop()
  }

  sender_celltypes <- unique(obj@meta.data[,group.column])
  list_expressed_genes_sender <- list()
  for(i in unique(sender_celltypes)){
    k <- seurat.obj@meta.data[seurat.obj@meta.data[,group.column]==i,]
    k <- rownames(k)
    m <- seurat.obj@assays$SCT@counts[,colnames(seurat.obj@assays$SCT@counts)%in%k]
    v1 <- apply(m,1,function(x){length(x[x>0])})
    min.pct.cells <- round(length(k)*min.pct,0)
    v1 <- v1[v1>=min.pct.cells]
    if(length(v1)<1){next}else{
      v1 <- names(v1)
      m <- seurat.obj@assays$SCT@counts[rownames(seurat.obj@assays$SCT@counts)%in%v1,
                                        colnames(seurat.obj@assays$SCT@counts)%in%k]
      k2 <- rowMeans(m)
      k2 <- k2[k2>=expr.filter]
      if(length(k2)<1){next}else{
        list_expressed_genes_sender[[i]] <- names(k2)
      }
    }
  }

  expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

  if(length(expressed_genes_sender)<1){
    cat(crayon::red("\nno gene passes the threshold, execution stops \n"))
    cat(crayon::red("\ntry changing the min.pct or the expr.filter parameter \n"))
    stop()
  }

  ligands <- new_lr %>% pull(from) %>% unique()
  expressed_ligands <- intersect(ligands,expressed_genes_sender)

  if(length(expressed_ligands)<1){
    cat(crayon::red("\nno gene passes the threshold, execution stops \n"))
    cat(crayon::red("\ntry changing the min.pct or the expr.filter parameter \n"))
    stop()
  }
  if(isTRUE(verbose)){cat(crayon::green("done \n"))}

  #isolation of valid ECM receptors (should be expressed at least by 30% of a population but can be changed via the min.pct parameter)
  if(isTRUE(verbose)){cat(crayon::white("(3/9) evaluating matrisome receptor genes in target population(s)... "))}

  if(!is.null(target)){receiver_celltypes <- target}else{
    receiver_celltypes <- unique(obj@meta.data[,group.column])}

  list_expressed_genes_receiver <- list()
  for(i in unique(receiver_celltypes)){
    k <- seurat.obj@meta.data[seurat.obj@meta.data[,group.column]==i,]
    k <- rownames(k)
    m <- seurat.obj@assays$SCT@counts[,colnames(seurat.obj@assays$SCT@counts)%in%k]
    v1 <- apply(m,1,function(x){length(x[x>0])})
    min.pct.cells <- round(length(k)*min.pct,0)
    v1 <- v1[v1>=min.pct.cells]
    if(length(v1)<1){next}else{
      v1 <- names(v1)
      m <- seurat.obj@assays$SCT@counts[rownames(seurat.obj@assays$SCT@counts)%in%v1,
                                        colnames(seurat.obj@assays$SCT@counts)%in%k]
      k2 <- rowMeans(m)
      k2 <- k2[k2>=expr.filter]
      if(length(k2)<1){next}else{
        list_expressed_genes_receiver[[i]] <- names(k2)
      }
    }
  }
  background_expressed_genes <- lapply(list_expressed_genes_receiver,function(x){
    x[x%in%rownames(new_ligand_target_matrix)]
  })
  nnn <- unique(unlist(background_expressed_genes))

  if(length(nnn)<1){
    cat(crayon::red("\nno gene passes the threshold, execution stops \n"))
    cat(crayon::red("\ntry changing the min.pct or the expr.filter parameter \n"))
    stop()
  }

  receptors <- new_lr %>% pull(to) %>% unique()
  expressed_receptors <- lapply(background_expressed_genes,function(x){
    intersect(receptors,x)})
  nnn <- unique(unlist(expressed_receptors))

  if(length(nnn)<1){
    cat(crayon::red("\nno gene passes the threshold, execution stops \n"))
    cat(crayon::red("\ntry changing the min.pct or the expr.filter parameter \n"))
    stop()
  }
  if(isTRUE(verbose)){cat(crayon::green("done \n"))}

  #filtering potentially active ligands
  if(isTRUE(verbose)){cat(crayon::white("(4/9) extracting potential matrisome ligands based on receptors... "))}
  potential_ligands <- new_lr %>% filter(from %in% expressed_ligands & to %in% nnn) %>% pull(from) %>% unique()

  if(length(potential_ligands)<1){
    cat(crayon::red("\nno cell-specific pair was found, execution stops \n"))
    cat(crayon::red("\ntry changing the min.pct or the expr.filter parameter \n"))
    stop()
  }

  lndf <- list()
  for(i in names(expressed_receptors)){
    y <- expressed_receptors[i][[1]]
    nr <- new_lr[new_lr$from%in%potential_ligands &
                   new_lr$to%in%y,]
    if(nrow(nr)<1){next}else{
      nr$target.cell <- i
      lndf[[i]] <- nr
    }
  }
  nnn <- bind_rows(lndf)

  if(nrow(nnn)<1){
    cat(crayon::red("\nno path explains ligand activities in cells, execution stops \n"))
    cat(crayon::red("\ntry changing the min.pct or the expr.filter parameter \n"))
    stop()
  }else{
    potential_ligands <- unique(nnn$from)
    expressed_receptors <- expressed_receptors[unique(nnn$target.cell)]
  }
  if(isTRUE(verbose)){cat(crayon::green("done \n"))}

  #extracting response signatures
  if(!is.null(geneset)){
    if(isTRUE(verbose)){cat(crayon::white("(5/9) applying the response geneset... "))}
    ls_of_interest <- list()
    ls_of_interest[[1]] <- geneset
  }else{
    if(isTRUE(verbose)){cat(crayon::white("(5/9) creating population-specific response geneset(s)... "))}

    ls_of_interest <- list()
    for(i in names(expressed_receptors)){
      z <- expressed_receptors[i][[1]]
      if(length(z)<1){next}else{
        g <- background_expressed_genes[i][[1]]
        g <- g[!g%in%z]
        ln <- table(obj@meta.data[,group.column])
        ln <- ln[names(ln)%in%i]
        if(ln<10){next}else{
          seurat_obj_receiver <- SetIdent(obj,value=ifelse(Idents(obj)%in%i,"1","0"))
          seurat_obj_receiver <- subset(seurat_obj_receiver,features=g)
          DE_table_receiver <- FindMarkers(object = seurat_obj_receiver, ident.1 = "1", ident.2 = "0", min.pct = min.pct, max.cells.per.ident = ln,verbose = F)
          DE_table_receiver$gene <- rownames(DE_table_receiver)
          geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & avg_log2FC >= 0.25) %>% pull(gene)
          geneset_oi <- geneset_oi %>% .[. %in% rownames(new_ligand_target_matrix)]
          lll <- list_expressed_genes_receiver[i][[1]]
          geneset_oi <- geneset_oi[geneset_oi%in%lll]
          ls_of_interest[[i]] <- geneset_oi
        }
      }
    }
  }
  nnn <- unique(unlist(ls_of_interest))

  if(length(nnn)<1){
    cat(crayon::red("no gene for geneset found, execution stops \n"))
    cat(crayon::red("try changing the min.pct or the expr.filter parameter \n"))
    stop()
  }
  if(isTRUE(verbose)){cat(crayon::green("done \n"))}

  #filtering potentially active ligands
  if(isTRUE(verbose)){cat(crayon::white("(6/9) testing potential ligand(s) against response genesets (clipping of ligands might occurr if > max.ligands)... "))}

  ligand_activities <- list()
  for(i in names(ls_of_interest)){
    act <- suppressWarnings(predict_ligand_activities(geneset = ls_of_interest[i][[1]], background_expressed_genes = background_expressed_genes[i][[1]], ligand_target_matrix = new_ligand_target_matrix, potential_ligands = potential_ligands[potential_ligands%in%colnames(new_ligand_target_matrix)]))
    act <- act %>% arrange(-pearson)
    act <- act[act$pearson>0,]
    if(nrow(act)<1){next}else{
      if(nrow(act)>max.ligands){
        act<-act[1:max.ligands,]
      }else{act<-act}
      act$target.cell <- i
      ligand_activities[[i]] <- act}
  }
  nnn <- bind_rows(ligand_activities)
  nnn <- unique(nnn$test_ligand)
  nnn <- as.character(na.omit(nnn))

  if(length(nnn)<1){
    cat(crayon::red("\nno ligand explains geneset(s), execution stops \n"))
    cat(crayon::red("\ntry changing the min.pct or the expr.filter parameter \n"))
    stop()
  }
  if(isTRUE(verbose)){cat(crayon::green("done \n"))}

  #combining the results
  if(isTRUE(verbose)){cat(crayon::white("(7/9) forming triplets (clipping of targets might occurr if > max.targets)... "))}

  best_upstream_ligands <- as.character(na.omit(unique(bind_rows(ligand_activities)$test_ligand)))

  sends <- list()
  for(i in names(list_expressed_genes_sender)){
    z <- list_expressed_genes_sender[i][[1]]
    z <- z[z%in%best_upstream_ligands]
    if(length(z)<1){next}else{
      sends[[i]] <- data.frame(origin.population=i,ECM.component=z)
    }
  }
  sends <- bind_rows(sends)
  avg_expression_ligands <- suppressWarnings(as.data.frame(t(unlist(AverageExpression(obj, assays = "SCT" ,features = best_upstream_ligands)[[1]]))))
  nm <- rownames(avg_expression_ligands)
  if(length(avg_expression_ligands)>1){
    avg_expression_ligands <- suppressMessages(reshape2::melt(avg_expression_ligands))
    avg_expression_ligands$origin.population <- nm
    names(avg_expression_ligands)[1] <- "ECM.component"
  }else{
    avg_expression_ligands <- data.frame(ECM.component=best_upstream_ligands,
                                         value=avg_expression_ligands,
                                         origin.population=nm)
    names(avg_expression_ligands)[2] <- "value"
  }
  sends <- distinct(merge(sends,avg_expression_ligands,by=c("origin.population","ECM.component")))
  sends <- subset(sends,sends$value>expr.filter)

  if(nrow(sends)<1){
    cat(crayon::red("\nno population expresses any ligand above 1, execution stops \n"))
    cat(crayon::red("\ntry changing the min.pct or the expr.filter parameter \n"))
    stop()
  }else{
    names(sends)[3] <- "ECM.value"}

  lndf <- bind_rows(lndf)
  k <- list()
  for(i in names(ligand_activities)){
    z <- ligand_activities[i][[1]]
    r <- subset(lndf,lndf$target.cell==i)
    r <- subset(r,r$from%in%unique(z$test_ligand))
    if(nrow(r)<1){next}else{
      avg_expression_r <- suppressWarnings(as.data.frame(t(unlist(AverageExpression(obj,assays="SCT",features = unique(r$to))[[1]]))))
      avg_expression_r <- avg_expression_r[rownames(avg_expression_r)%in%i,]
      if(length(avg_expression_r)>1){
        nm <- rownames(avg_expression_r)
        avg_expression_r <- suppressMessages(reshape2::melt(avg_expression_r))
        avg_expression_r$receiving.population <- nm
        names(avg_expression_r)[1] <- "to"
      }else{
        avg_expression_r <- data.frame(to=r$to,
                                       value=avg_expression_r,
                                       receiving.population=i)

      }
      r <- r[,c(1,2,5)]
      names(r)[3] <- "receiving.population"
      r <- distinct(merge(r,avg_expression_r,by=c("to","receiving.population")))
      r <- subset(r,r$value>expr.filter)
      if(nrow(r)<1){next}else{
        r <- r[,c(3,1,2,4)]
        k[[i]] <- r}
    }}
  k <- bind_rows(k)
  names(k)[1:2] <- c("ECM.component","ECM.receptor")
  names(k)[4] <- "receptor.value"

  if(nrow(k)<1){
    cat(crayon::red("\nno population expresses any receptor above 1, execution stops \n"))
    cat(crayon::red("\ntry changing the min.pct or the expr.filter parameter \n"))
    stop()
  }

  k2 <- list()
  for(i in names(ligand_activities)){
    z <- ls_of_interest[i][[1]]
    avg_expression_t <- suppressWarnings(as.data.frame(t(unlist(AverageExpression(obj,assays = "SCT" ,features = rownames(z))[[1]]))))
    avg_expression_t <- avg_expression_t[rownames(avg_expression_t)%in%i,]
    if(length(avg_expression_t)>1){
      avg_expression_t <- suppressMessages(reshape2::melt(avg_expression_t))
      avg_expression_t$receiving.population <- i
      names(avg_expression_t)[1] <- "target.gene"
    }else{
      avg_expression_t <- data.frame(target.gene=rownames(z),
                                     value=avg_expression_t,
                                     receiving.population=i)

    }
    df <- subset(avg_expression_t,avg_expression_t$value>expr.filter)
    if(nrow(df)<1){next}else{
      df <- df[order(-df$value),]
      if(nrow(df)>max.targets){df<-df[1:max.targets,]}else{df<-df}
      k2[[i]] <- df}
  }
  k2 <- bind_rows(k2)
  names(k2)[2] <- "target.value"

  if(nrow(k2)<1){
    cat(crayon::red("\nno population expresses any target above 1, execution stops \n"))
    cat(crayon::red("\ntry changing the min.pct or the expr.filter parameter \n"))
    stop()
  }

  fin <- distinct(merge(sends,k,by="ECM.component"))
  fin <- fin[,c(2,1,3,5,4,6)]
  fin <- distinct(merge(fin,k2,by="receiving.population"))
  fin <- fin[,c(2:4,1,5,6:8)]

  #cut to only in-cluster communications
  fin <- fin[fin$origin.population==fin$receiving.population,]
  if(nrow(fin)<1){
    cat(crayon::red("\nno within-cluster ECM-cell communication (trans-community communication is not yet supported), execution stops \n"))
    cat(crayon::red("\ntry changing the min.pct or the expr.filter parameter \n"))
    stop()
  }
  if(isTRUE(verbose)){cat(crayon::green("done \n"))}

  #explainable paths
  ### quick differential correlation - VERSION 2

  if(!isFALSE(explainable)){
    if(isTRUE(verbose)){cat(crayon::white(paste0("(8/9) finding explainable paths... ")))}

    om <- as.matrix(obj@assays$SCT@counts)
    om[is.na(om)] <- 0
    rands <- list()
    for(i in 1:10){
      rm <- om[sample(rownames(om),100,replace = T),
               sample(colnames(om),100,replace = T)]
      cf <- mean(rowMeans(suppressWarnings(cor(rm,method = "spearman"))))
      cf <- ifelse(is.na(cf),0,cf)
      rands[[i]] <- cf
    }

    rands <- mean(unlist(rands))

    rp <- as.character(unique(fin$receiving.population))
    l1 <- list()
    for(i in rp){
      z <- subset(fin,fin$receiving.population==i)
      g <- as.character(unique(z$target.gene))
      rs <- as.character(unique(z$ECM.receptor))
      s <- obj@meta.data[obj@meta.data[,group.column]==i,]
      s <- unique(rownames(s))
      if(length(s)<1){next}else{
        #s2 <- obj@meta.data[obj@meta.data[,group.column]!=i,]
        #s2 <- unique(rownames(s2))
        tarrec.m <- om[rownames(om)%in%unique(c(g,rs)),
                       colnames(om)%in%s]
        mat1 <- suppressWarnings(cor(t(tarrec.m),method="spearman"))
        mat1[upper.tri(mat1)] <- NA
        rw <- rownames(mat1)
        cl <- colnames(mat1)
        mat1 <- as.data.frame(cbind(which(!is.na(mat1),arr.ind = TRUE),na.omit(as.vector(mat1))))
        mat1 <- mat1[mat1$V3>rands,]
        mat1 <- mat1[mat1$row!=mat1$col,]
        if(nrow(mat1)<1){next}else{
          mat1$row <- rw[mat1$row]
          mat1$col <- cl[mat1$col]
          names(mat1) <- c("ECM.receptor","target.gene")
          tf <- new_gr[new_gr$to%in%mat1$target.gene,]
          rs <- as.character(unique(tf$from))
          #s <- obj@meta.data[obj@meta.data[,group.column]==i,]
          #s <- unique(rownames(s))
          #s2 <- obj@meta.data[obj@meta.data[,group.column]!=i,]
          #s2 <- unique(rownames(s2))
          tartf.m <- om[rownames(om)%in%unique(c(mat1$target.gene,rs)),
                        colnames(om)%in%s]
          mat2 <- suppressWarnings(cor(t(tartf.m),method="spearman"))
          mat2[upper.tri(mat2)] <- NA
          rw <- rownames(mat2)
          cl <- colnames(mat2)
          mat2 <- as.data.frame(cbind(which(!is.na(mat2),arr.ind = TRUE),na.omit(as.vector(mat2))))
          mat2 <- mat2[mat2$V3>rands,]
          mat2 <- mat2[mat2$row!=mat2$col,]
          if(nrow(mat2)<1){next}else{
            mat2$row <- rw[mat2$row]
            mat2$col <- cl[mat2$col]
            names(mat2) <- c("TF","target.gene")
            fifin <- distinct(merge(mat1,mat2,by="target.gene"))
            if(nrow(fifin)<1){next}else{
              nobj <- om[rownames(om)%in%unique(fifin$TF),
                         colnames(om)%in%s]
              df <- data.frame(TF=rownames(nobj),TF.value=rowMeans(nobj))
              df$receiving.population <- i
              df <- merge(fifin,df,by=c("TF"))
              df <- df[,c(1,2,3,6,7)]
              l1[[i]] <- df
            }
          }
        }
      }}
    l1 <- bind_rows(l1)
    rownames(l1) <- NULL

    new_fin <- distinct(merge(fin,l1,by=c("receiving.population","target.gene","ECM.receptor")))
    new_fin <- new_fin[,c(4:6,1,3,7,2,8,9,10)]
  }else{
    if(isTRUE(verbose)){cat(crayon::white("(8/9) not searching for explainable paths \n"))}
    new_fin <- fin
    new_fin$TF <- "not searched"
    new_fin$TF.value <- 1
  }

  if(isTRUE(add.RF)){
    rlst <- list()
    for(i in unique(new_fin$receiving.population)){
      z <- subset(new_fin,new_fin$receiving.population==i)
      z <- unique(c(z$ECM.receptor,z$TF))
      z <- z[z!="not searched"]
      if(length(z)<1){next}else{
        mk <- obj@assays$SCT@counts[rownames(obj@assays$SCT@counts)%in%z,]
        Y <- rownames(obj@meta.data[obj@meta.data[,group.column]==i,])
        Y <- ifelse(colnames(mk)%in%Y,"first","second")
        Y <- factor(Y,levels=c("first","second"))
        mk[is.na(mk)] <- 0
        mk <- t(mk)
        mk <- as.data.frame(apply(mk,2,scale))
        mk$identity <- Y

        set.seed(12345)
        sds <- vector(mode = "list", length = 5)
        inTrain <- createDataPartition(mk$identity, p = 0.7, list = FALSE)
        trainData <- mk[inTrain,]
        testData <- mk[-inTrain,]

        suppressMessages(suppressWarnings(
          my_control <- trainControl(
            method="boost",
            number=5,
            savePredictions="final",
            classProbs=TRUE,
            index=createResample(trainData$identity,5),
            summaryFunction=twoClassSummary
          )
        ))

        suppressMessages(suppressWarnings(
          rf_default <- train(identity~.,
                              data=mk,
                              method='rf',
                              metric='Accuracy',
                              trControl=my_control
          )
        ))

        vs <- vs <- as.data.frame(varImp(rf_default)$importance)
        vs$gene <- rownames(vs)
        vs <- vs[order(-vs$Overall),]
        if(nrow(vs)>max.ligands){
          vs <- vs[1:max.ligands,]
        }else{vs <- vs}

        fincut <- subset(new_fin,new_fin$receiving.population==i)
        fincut <- fincut[fincut$ECM.receptor%in%vs$gene & fincut$TF%in%vs$gene,]

        if(nrow(fincut)<1){
          fincut <- subset(new_fin,new_fin$receiving.population==i)
          fincut <- fincut[fincut$ECM.receptor%in%vs$gene | fincut$TF%in%vs$gene,]

          if(nrow(fincut)<1){

          }else{
            rlst[[i]] <- fincut}
        }

      }
    }
    rlst <- bind_rows(rlst)
    new_fin <- rlst
  }

  if(nrow(new_fin)<1){
    cat(crayon::red("\nno ligand path was found, execution stops \n"))
    cat(crayon::red("\ntry changing the min.pct or the expr.filter parameter or to set add.RF to FALSE\n"))
    stop()
  }
  if(isTRUE(verbose)){cat(crayon::green("done \n"))}

  #export the results

  if(isTRUE(scale.results)){
    if(isTRUE(verbose)){cat(crayon::white("(9/9) scaling the results to [0,1] and annotating... "))}
    df_scaled <- new_fin %>% group_by(origin.population) %>% mutate(ECM.value=round(rescale(ECM.value,to=c(0.01,1)),3)) %>%
      group_by(receiving.population) %>% mutate(receptor.value=round(rescale(receptor.value,to=c(0.01,1)),3)) %>%
      group_by(receiving.population) %>% mutate(target.value=round(rescale(target.value,to=c(0.01,1)),3)) %>%
      group_by(receiving.population) %>% mutate(TF.value=round(rescale(TF.value,to=c(0.01,1)),3))
    fin <- data.frame(df_scaled)

    annofin <- list()
    for(i in unique(fin$receiving.population)){
      z <- fin[fin$receiving.population==i,]
      z <- unique(z$ECM.receptor)
      l <- list()
      for(w in ncpl$complex){
        s <- subset(ncpl,ncpl$complex==w)
        s <- unlist(strsplit(s$genes,split=";"))
        y <- length(intersect(s,z))
        if(y<2){next}else{
          y <- intersect(s,z)
          l[[w]] <- data.frame(receiving.population=i,ECM.receptor=y,annotation=paste0(w," is enriched in the population"))}
      }
      l <- bind_rows(l)
      annofin[[i]] <- bind_rows(l)
    }
    annofin <- bind_rows(annofin)
    fin <- distinct(merge(fin,annofin,by=c("receiving.population","ECM.receptor"),all=T))
    fin[is.na(fin)] <- ""
    fin <- fin[,c(3:5,1,2,6:ncol(fin))]

    attributes(fin)$mode <- "matricom.scaled"}else{
      if(isTRUE(verbose)){cat(crayon::white("(9/9) annotating results... "))}
      fin<-new_fin

      annofin <- list()
      for(i in unique(fin$receiving.population)){
        z <- fin[fin$receiving.population==i,]
        z <- unique(z$ECM.receptor)
        l <- list()
        for(w in ncpl$complex){
          s <- subset(ncpl,ncpl$complex==w)
          s <- unlist(strsplit(s$genes,split=";"))
          y <- length(intersect(s,z))
          if(y<2){next}else{
            y <- intersect(s,z)
            l[[w]] <- data.frame(receiving.population=i,ECM.receptor=y,annotation=paste0(w," is enriched in the receiver population"))}
        }
        l <- bind_rows(l)
        annofin[[i]] <- bind_rows(l)
      }
      annofin <- bind_rows(annofin)

      if(nrow(annofin)<1){
        cat(crayon::red("no ligand path was found, execution stops \n"))
        cat(crayon::red("try changing the min.pct or the expr.filter parameter \n"))
        stop()
      }

      if(nrow(fin)<1){
        cat(crayon::red("no ligand path was found, execution stops \n"))
        cat(crayon::red("try changing the min.pct or the expr.filter parameter \n"))
        stop()
      }

      fin <- distinct(merge(fin,annofin,by=c("receiving.population","ECM.receptor"),all=T))
      fin[is.na(fin)] <- ""
      fin <- fin[,c(3:5,1,2,6:ncol(fin))]

      attributes(fin)$mode <- "matricom.raw"
    }

  if(isTRUE(verbose)){cat(crayon::green("done \n"))}
  return(fin)
}

#' Spatplot: representing cell-ECM communication in spatial RNAseq data

#' @param result.obj character. The results from a previous call to matricom.spatial.
#' @param seurat.obj character. A spatial Seurat object.
#' @param group.column character. Name of the column in metadata with cell identities to be used.
#' @param feature character. A quadruplet to be plotted, collapsed by underscores (e.g.,"gene_receptor_target_TF"). If NULL (default), the function will pick a random quadruplet.
#' @param singles logical. Whether to plot also spatial plots of the ECM/receptor/target genes. Default is FALSE

#' @return a multiplot with at least a spatial cluster plot, a spatial plot with clusters showing the chosen signal, and a spatial plot with average ECM/receptor/target intensities per dot. Might also include the spatial plots of the single elements.
#' @export

#' @examples spatplot(results,obj,"group.column")

spatplot <- function(result.obj=NULL,
                     seurat.obj=NULL,
                     group.column=NULL,
                     feature=NULL,
                     singles=FALSE){
  #error handling
  if(is.null(result.obj)){
    cat(crayon::red("no data provided, execution stops \n"))
    stop()
  }
  if(class(result.obj)!="data.frame"){
    cat(crayon::red("data are not in a valid format (data.frame), execution stops \n"))
    stop()
  }
  if(is.null(seurat.obj)){
    cat(crayon::red("no data provided, execution stops \n"))
    stop()
  }
  if(class(seurat.obj)!="Seurat"){
    cat(crayon::red("data are not a valid Seurat object, execution stops \n"))
    stop()
  }
  if(nrow(seurat.obj@meta.data)<1){
    cat(crayon::red("data are not a valid Seurat object, execution stops \n"))
    stop()
  }
  if(is.null(group.column)){
    cat(crayon::red("no group assignment column provided, execution stops \n"))
    stop()
  }

  require(Matrix)
  res <- seurat.obj

  result.obj$feat <- paste(res$ECM.component,
                           res$ECM.receptor,
                           res$target.gene,
                           res$TF,
                           sep="_")
  result.obj$feat.value <- rowMeans(result.obj[,c(3,6,8,10)])

  if(is.null(feature)){
    cat(crayon::white("no feature chosen, picking a random one for the graph \n"))
    sel <- sample(result.obj$feat,1)
  }else{sel <- feature}

  m0 <- as.matrix(seurat.obj@assays$SCT@counts)

  k <- data.frame(sample=rownames(seurat.obj@meta.data),
                  cl=seurat.obj@meta.data[,group.column])

  z <- distinct(result.obj[result.obj$feat==sel,])
  z <- z[,c(1,12,13)]
  k1 <- distinct(merge(z,k,by.x="origin.population",by.y="cl",all.y=T))
  k1$feat.value[is.na(k1$feat.value)] <- 0
  v <- as.data.frame(t(k1$feat.value))
  names(v) <- k1$sample
  tsel <- gsub("_not searched","",sel)
  rownames(v) <- tsel
  v <- v[,colnames(m0)]

  m0 <- rbind(m0,v)
  m0 <- as.matrix(m0)
  m0 <- as(m0,"sparseMatrix")

  seurat.obj@assays$SCT@data <- m0

  p1 <- suppressWarnings(SpatialDimPlot(seurat.obj, label = TRUE, label.size = 3)) +
    ggtitle("clusters") + ggplot2::theme(legend.position = "none")

  seurat.obj@meta.data$newclust <- ifelse(Idents(seurat.obj)%in%unique(z$origin.population),"1","0")
  Idents(seurat.obj) <- seurat.obj@meta.data$newclust
  p2 <- suppressWarnings(SpatialDimPlot(seurat.obj, label = F)+
                           scale_fill_manual(values = c("blue","red"))) +
    ggtitle("clusters with signal") + ggplot2::theme(legend.position = "none")

  p3 <- SpatialFeaturePlot(seurat.obj,
                           slot = "data",
                           features = tsel) +
    ggtitle("matRicom signal")

  tsel <- unlist(strsplit(tsel,"_"))

  if(isFALSE(singles)){
    suppressWarnings(plot(p1+p2+p3))
  }else{

    l <- list()
    for(i in 1:length(tsel)){
      stp <- i+3
      l[[stp]] <- SpatialFeaturePlot(seurat.obj, features = tsel[i])
    }
    cat(crayon::yellow("NOTE: gene data (counts) are not scaled like the matRicom results \n"))
    suppressWarnings(plot(p1+p2+p3+l))
  }

}

#' filter.finder: graphical utility to explore filters (expression and min.pct) for matRicom

#' @param seurat.obj character. A Seurat, Seurat-like object, or a spatial Seurat object.
#' @param group.column character. Name of the column in metadata with cell identities to be used.

#' @return a faceted plot showing mean gene expression and % expression per cell identity together with results from 100 random samples and the defulat thresholds.
#' @export

#' @examples filter.finder(seurat.obj,group.column)

filter.finder <- function(seurat.obj=NULL,
                          group.column=NULL){

  library(ggplot2)

  #error handling
  if(is.null(seurat.obj)){
    cat(crayon::red("no data provided, execution stops \n"))
    stop()
  }
  if(class(seurat.obj)!="Seurat"){
    cat(crayon::red("data are not a valid Seurat object, execution stops \n"))
    stop()
  }
  if(nrow(seurat.obj@meta.data)<1){
    cat(crayon::red("data are not a valid Seurat object, execution stops \n"))
    stop()
  }
  if(is.null(group.column)){
    cat(crayon::red("no group assignment column provided, execution stops \n"))
    stop()
  }

  #read in the matricom objects
  data(matricom_obj, package = "matRicom")
  new_ligand_target_matrix <- matricom_obj$new_ligand_target_matrix

  #expression levels and % comparison
  obj <- seurat.obj
  Idents(obj) <- obj@meta.data[,group.column]

  ls <- list()
  for(i in 1:100){

    if(isTRUE("Spatial" %in% names(obj@assays))){
      k <- sample(rownames(obj@assays$SCT@counts),100)
      g <- sample(colnames(obj@assays$SCT@counts),100)
      v <- obj@assays$SCT@counts[rownames(obj@assays$SCT@counts)%in%k,
                                 colnames(obj@assays$SCT@counts)%in%g]
    }else{
      k <- sample(rownames(obj@assays$RNA@counts),100)
      g <- sample(colnames(obj@assays$RNA@counts),100)
      v <- obj@assays$RNA@counts[rownames(obj@assays$RNA@counts)%in%k,
                                 colnames(obj@assays$RNA@counts)%in%g]
    }

    a <- apply(v,1,function(x){
      length(x[x>0])
    })
    a <- length(a[a>0])
    a <- a/100

    b <- mean(rowMeans(v))

    ls[[i]] <- data.frame(a=a,b=b)
  }
  ls <- bind_rows(ls)
  a <- mean(ls$a)
  b <- mean(ls$b)

  if(isTRUE("Spatial" %in% names(obj@assays))){
    m <- obj@assays$SCT@counts[rownames(obj@assays$SCT@counts)%in%colnames(new_ligand_target_matrix),]
  }else{
    m <- obj@assays$RNA@counts[rownames(obj@assays$RNA@counts)%in%colnames(new_ligand_target_matrix),]
  }

  ls2 <- list()
  for(i in unique(Idents(obj))){
    n <- as.character(unlist(obj@meta.data[group.column]))
    names(n) <- rownames(obj@meta.data)
    n <- n[n%in%i]
    n <- names(n)

    if(isTRUE("Spatial" %in% names(obj@assays))){
      v <- m[,colnames(obj@assays$SCT@counts)%in%n]
    }else{
      v <- m[,colnames(obj@assays$RNA@counts)%in%n]
    }

    if(nrow(v)<2){next}else{
      if(ncol(v)<2){next}else{

        mp <- apply(v,1,function(x){
          length(x[x>0])
        })
        mp <- mp/length(n)
        mg <- rowMeans(v)

        ls2[[i]] <- data.frame(mean.percentile=mp,mean.expression=mg,cell.type=i)
      }
    }
  }
  ls2 <- bind_rows(ls2)

  if(isFALSE("Spatial" %in% names(obj@assays))){
    p <- ggplot(ls2,aes(x=mean.percentile,y=mean.expression)) +
      geom_point() +
      theme_bw() + xlab("mean % positivity") + ylab("mean gene expression value") +
      facet_wrap(~unique(ls2$cell.type)) +
      geom_vline(xintercept = rep(b,length(unique(ls2$cell.type))), color="black", linetype=2) +
      geom_hline(yintercept = rep(a,length(unique(ls2$cell.type))), color="black", linetype=2) +
      geom_vline(xintercept = rep(0.3,length(unique(ls2$cell.type))), color="red", linetype=2) +
      geom_hline(yintercept = rep(1,length(unique(ls2$cell.type))), color="red", linetype=2) +
      xlim(c(0,1)) +
      ggtitle("gene expression and % positivity vs. \nBLACK: random 100x100 tables (mean of 100 draws) and \nRED: default matRicom thresholds")
    suppressWarnings(print(p))}else{
      p <- ggplot(ls2,aes(x=mean.percentile,y=mean.expression)) +
        geom_point() +
        theme_bw() + xlab("mean % positivity") + ylab("mean gene expression value") +
        facet_wrap(~unique(ls2$cell.type)) +
        geom_vline(xintercept = rep(b,length(unique(ls2$cell.type))), color="black", linetype=2) +
        geom_hline(yintercept = rep(a,length(unique(ls2$cell.type))), color="black", linetype=2) +
        geom_vline(xintercept = rep(0.5,length(unique(ls2$cell.type))), color="red", linetype=2) +
        geom_hline(yintercept = rep(1,length(unique(ls2$cell.type))), color="red", linetype=2) +
        xlim(c(0,1)) +
        ggtitle("gene expression and % positivity vs. \nBLACK: random 100x100 tables (mean of 100 draws) and \nRED: default matRicom thresholds")
      suppressWarnings(print(p))
    }
}

#' make.object: a wrapper around the Seurat workflow to prepare 10X data for matRicom

#' @param data.name character. Internal data name.
#' @param data.dir character. The full path to the directory containing the objects needed for the Load10X function in Seurat
#' @param project.name character. How the project will be called. If NULL, data.name are used.
#' @param data.multiple logical. IWhether data contain multiple types. If TRUE, "Gene Expression" is used. Default is FALSE.
#' @param dirty.trick logical. Whether to turn gene names to uppercase in case of non-human data (matRicom would otherwise crash). Default is FALSE.
#' @param min.cells numeric. The minimum number of cells expressing a given feature.
#' @param min.features numeric. The minimum number of features per cell.
#' @param percent.mt numeric. Cutoff for mitochondrial RNA content.
#'
#' @param gene.column.rx numeric. Option gene.column of Read10X (defaults to 2)
#' @param cell.column.rx numeric. Option cell.column of Read10X (defaults to 1)
#' 
#' @return a Seurat object
#' @export

#' @examples make.object(data.name,data.dir)

make.object <- function(data.name = NULL,
                        data.dir = NULL,
                        project.name = NULL,
                        gene.column.rx = 2,
                        cell.column.rx = 1,
                        data.multiple = FALSE,
                        dirty.trick = FALSE,
                        min.cells = 3,
                        min.features = 200,
                        percent.mt = 5,
                        verbose = TRUE){

  #error handling
  if(is.null(data.name)){
    cat(crayon::red("no data provided, execution stops \n"))
    stop()
  }

  if(is.null(data.dir)){
    cat(crayon::red("no data provided, execution stops \n"))
    stop()
  }

  # set defaults if no values provided
  if(is.null(min.cells)){
    mincells <- 3
    if(isTRUE(verbose)){
    cat(crayon::white("setting min.cells = 3 (default) \n"))}
  }else{
    mincells <- min.cells
    if(isTRUE(verbose)){
    cat(crayon::white("setting min.cells = ",mincells, "\n"))}
  }

  if(is.null(min.features)){
    minfeatures <- 200
    if(isTRUE(verbose)){
    cat(crayon::white("setting min.features = 200 (default) \n"))}
  }else{
    minfeatures <- min.features
    if(isTRUE(verbose)){
    cat(crayon::white("setting min.features = ",minfeatures, "\n"))}
  }

  if(is.null(percent.mt)){
    percentmt <- 5
    if(isTRUE(verbose)){
    cat(crayon::white("setting percent.mt = 5 (default) \n"))}
  }else{
    percentmt <- percent.mt
    if(isTRUE(verbose)){
    cat(crayon::white("setting percent.mt = ",percentmt, "\n"))}
  }

  if(is.null(project.name)){
    project.name <- data.name
  }else{
    project.name <- project.name
  }

  #RDS object from 10xGenomics
  library(dplyr)
  library(Seurat)
  library(patchwork)

  if(isTRUE(verbose)){
  cat(crayon::white("1/3) loading data... "))}

  ddd <- data.dir

  tengen.data <- Read10X(ddd, gene.column = gene.column.rx, cell.column = cell.column.rx)

  # are we using data with more than one type?
  if(isTRUE(data.multiple)){
    counts.data <- tengen.data$`Gene Expression`
  }else{
    counts.data <- tengen.data
  }

  if(isTRUE(dirty.trick)){
    if(isTRUE(verbose)){
    cat(crayon::italic(" applying dirty trick for non-human data "))}
    rownames(counts.data) <- casefold(rownames(counts.data), upper = T)
  }
  if(isTRUE(verbose)){
  cat(crayon::green("done \n"))}

  tengen <- suppressMessages(CreateSeuratObject(counts = counts.data, project.name, min.cells = mincells, min.features = minfeatures))
  tengen[["percent.mt"]] <- suppressMessages(PercentageFeatureSet(tengen, pattern = "^MT-"))
  if(isTRUE(dirty.trick)){
    if(isTRUE(verbose)){
      cat(crayon::white("2/3) preparing the Seurat object, mitochondrial removal for non-human data might fail... "))}
    rownames(counts.data) <- casefold(rownames(counts.data), upper = T)
  }
  if(isTRUE(verbose)){
    cat(crayon::white("2/3) preparing the Seurat object... "))}
  tengen <- suppressMessages(subset(tengen, subset = nFeature_RNA > minfeatures & nFeature_RNA < 2500 & percent.mt < percentmt))
  tengen <- suppressMessages(NormalizeData(tengen, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F))
  #tengen <- suppressMessages(NormalizeData(tengen)) #done it already just above
  tengen <- suppressMessages(FindVariableFeatures(tengen, selection.method = "vst", nfeatures = 2000, verbose = F))
  all.genes <- rownames(tengen)
  tengen <- suppressMessages(ScaleData(tengen, features = all.genes,verbose = F))
  tengen <- suppressMessages(RunPCA(tengen, features = VariableFeatures(object = tengen),verbose = F))
  tengen <- suppressMessages(FindNeighbors(tengen, dims = 1:10,verbose = F))
  tengen <- suppressMessages(FindClusters(tengen, resolution = 0.5,verbose = F))
  tengen <- suppressMessages(RunUMAP(tengen, dims = 1:10,verbose = F))
  if(isTRUE(verbose)){
  cat(crayon::green("done \n"))}

  if(isTRUE(verbose)){
  cat(crayon::white("3/3) exporting the Seurat object... "))}
  if(isTRUE(verbose)){
  cat(crayon::green("done \n"))}
  return(tengen)
}

#' make.object.spatial: a wrapper around the Seurat workflow to prepare spatial 10X data for matRicom

#' @param data.name character. Internal data name.
#' @param data.dir character. The full path to the directory containing the objects needed for the Load10X_spatial function in Seurat (including h5 file). The folder must contain the h5 file plus a subfolder with the other requested data.
#' @param project.name character. How the project will be called. If NULL, data.name are used.
#' @param h5.filename character. Name (only the name, without full path) of the h5 file. Must contain the file extension.
#' @param dirty.trick logical. Whether to turn gene names to uppercase in case of non-human data (matRicom would othewise crash)
#' @param min.cells numeric. The minimum number of cells expressing a given feature.
#' @param min.features numeric. The minimum number of features per cell.
#' @param percent.mt numeric. Cutoff for mitochondrial RNA content.

#' @return a Seurat spatial object
#' @export

#' @examples make.object(data.dir,h5.filename)

make.object.spatial <- function(data.name = NULL,
                                data.dir = NULL,
                                project.name = NULL,
                                h5.filename = NULL,
                                slice.number = "slice1",
                                dirty.trick = FALSE,
                                percent.mt = 5,
                                verbose = TRUE){

  # RDS object from 10xGenomics
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(rhdf5) #needed for Linux users
  library(hdf5r) #needed for Linux users
  #hdf5 is also needed for Linux users (outside R)

  # error handling
  if(is.null(data.name)){
    cat(crayon::red("no data provided, execution stops \n"))
    stop()
  }
  if(is.null(data.dir)){
    cat(crayon::red("no data provided, execution stops \n"))
    stop()
  }
  if(is.null(h5.filename)){
    cat(crayon::red("no h5 file name provided, execution stops \n"))
    stop()
  }

  if(is.null(percent.mt)){
    percentmt <- 5
    if(isTRUE(verbose)){
      cat(crayon::white("setting percent.mt = 5 (default) \n"))}
  }else{
    percentmt <- percent.mt
    if(isTRUE(verbose)){
      cat(crayon::white("setting percent.mt = ",percentmt, "\n"))}
  }

  if(is.null(project.name)){
    project.name <- data.name
  }else{
    project.name <- project.name
  }

  cat(crayon::white("1/3) loading data... "))

  ddd <- data.dir
  hff <- h5.filename

  visium <- Load10X_Spatial(data.dir = ddd, filename = hff, assay = "Spatial", slice = slice.number, filter.matrix = TRUE, to.upper = dirty.trick)

  # set data and project names
  visium$orig.ident <- data.name
  visium@project.name <- project.name

  if(isTRUE(verbose)){
    cat(crayon::green("done \n"))}

  if(isTRUE(dirty.trick)){
    if(isTRUE(verbose)){
      cat(crayon::white("2/3) preparing the spatial Seurat object, mitochondrial removal for non-human data might fail... "))}
  }else{
    if(isTRUE(verbose)){
      cat(crayon::white("2/3) preparing the spatial Seurat object... "))}

  }

  visium[["percent.mt"]] <- PercentageFeatureSet(visium, pattern = "^MT-")
  visium <- subset(visium, subset = percent.mt < percentmt)

  # normalizes the data, detects high-variance features, and stores the data in the SCT assay
  visium <- suppressMessages(SCTransform(visium, assay = "Spatial", verbose = FALSE))

  # dimensionality reduction and clustering on the RNA expression
  visium <- suppressMessages(RunPCA(visium, assay = "SCT", verbose = FALSE))
  visium <- suppressMessages(FindNeighbors(visium, reduction = "pca", dims = 1:30))
  visium <- suppressMessages(FindClusters(visium, verbose = FALSE))
  visium <- suppressMessages(RunUMAP(visium, reduction = "pca", dims = 1:30, verbose=FALSE))
  if(isTRUE(verbose)){
    cat(crayon::green("done \n"))}
  if(isTRUE(verbose)){
    cat(crayon::white("3/3) exporting the Seurat spatial object... "))}
  if(isTRUE(verbose)){
    cat(crayon::green("done \n"))}
  return(visium)
}
