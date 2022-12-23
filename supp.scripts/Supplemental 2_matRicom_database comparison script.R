### Source evaluation script ###

library(dplyr)
library(data.table)
library(nichenetr)
library(igraph)
library(eulerr)
library(matRicom)

df1 <- as.data.frame(nichenetr::lr_network)
matricom_obj <- matRicom::matricom_obj
new_lr <- matricom_obj$new_lr
mat <- matRicom::mat
mat <- subset(mat,mat$category=="CORE MATRISOME")

length(unique(df1$from[df1$from%in%unique(mat$gene)])) #43 ECM components as ligands
length(unique(new_lr$from[new_lr$from%in%unique(mat$gene)])) #97 ECM components as ligands (2.25 times more)
a <- c(1:43)
b <- c(1:97)

eu <- euler(list(original=a,
                 updated=b))
plot(eu,
     quantities=T,
     shape="ellipse",
     fills=c("lavenderblush2","lightblue2"),
     labels=c("ECM components as ligands in NicheNet","ECM components as ligands in matRicom"))

a <- df1[df1$from%in%unique(mat$gene),]
a <- paste0(a$from,"_",a$to)
a <- unique(a) #347 unique ECM-receptor interactions (from ECM) 
b <- new_lr[new_lr$from%in%unique(mat$gene),]
b <- paste0(b$from,"_",b$to)
b <- unique(b) #1374 unique ECM-receptor interactions (from ECM) (3.95 times more)
eu <- euler(list(original=a,
                 updated=b))
plot(eu,
     quantities=T,
     shape="ellipse",
     fills=c("lavenderblush2","lightblue2"),
     labels=c("ECM interactions in NicheNet","ECM interactions in matRicom"))

### Graphs ###

library(ggplot2)

x <- data.frame(source=c("NicheNet","matRicom"),value=c(43,97))
x$source <- factor(x$source,levels=c("NicheNet","matRicom"))
ggplot(x,aes(source,value,fill=source)) +
  geom_bar(stat="Identity") +
  scale_fill_manual(values=c("lavenderblush2","lightblue2")) +
  theme_bw() + xlab("") + ylab("ECM components \nas ligands")

x <- data.frame(source=c("NicheNet","matRicom"),value=c(347,1374))
x$source <- factor(x$source,levels=c("NicheNet","matRicom"))
ggplot(x,aes(source,value,fill=source)) +
  geom_bar(stat="Identity") +
  scale_fill_manual(values=c("lavenderblush2","lightblue2")) +
  theme_bw() + xlab("") + ylab("ECM-receptor \ncomplexes")
