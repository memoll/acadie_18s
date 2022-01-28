###################################################################
# Network analysis (nematode-bacteria co-occurrence networks)     #
# Null models and network metrics                                 #
# Data: Miseq-18S & 16S - L'Acadie (ACA)                          #
# Mona Parizadeh & Mathieu landry - 2020-2021                     #
###################################################################

#Load libraries ####
library(rgr); packageVersion("rgr") #‘1.1.15’
library(phyloseq); packageVersion("phyloseq") #‘1.34.0’
library(ggplot2); packageVersion("ggplot2") #‘3.3.3’
library(picante); packageVersion("picante") #‘1.8.2’
library(ape); packageVersion("ape") #‘5.4.1’
library(SpiecEasi); packageVersion("SpiecEasi") #1.1.1
library(dplyr); packageVersion("dplyr") #‘1.0.4’
library(Matrix); packageVersion("Matrix") #‘1.3.2’
library(igraph); packageVersion("igraph") #‘1.2.6’
library(ggnetwork); packageVersion("ggnetwork") #‘0.5.10’

#load inputs ####
#inputs: phyloseq object and sparcc bootstrap
ps.ctl = readRDS("~/Documents/article2/ps.ctl5.rds"); ps.ctl
sparcc.boot.ctl1000 = readRDS("~/Documents/article2/sparcc.boot.ctl1000.rds")
ps.neo = readRDS("~/Documents/article2/ps.neo5.rds"); ps.neo
sparcc.boot.neo1000 = readRDS("~/Documents/article2/sparcc.boot.neo1000.rds")

#Replace NAs ####
#%Control####
comm.ctl = otu_table(ps.ctl)
taxo.ctl = tax_table(ps.ctl)
meta.ctl = sample_data(ps.ctl)
#replace unassigned by the lowest  assigned taxonomic level
for(i in c(1:length(colnames(tax_table(taxo.ctl))))){
  for(j in rownames(tax_table(taxo.ctl))){
    if(is.na(tax_table(taxo.ctl)[j,i]) == TRUE | tax_table(taxo.ctl)[j,i] == "Unknown_Family"){
      if(substr(tax_table(taxo.ctl)[j,i-1], start=1,stop=2)=="UA"){
        tax_table(taxo.ctl)[j,i] <- tax_table(taxo.ctl)[j,i-1]
      }
      else{
        tax_table(taxo.ctl)[j,i] <- paste0("UA_",tax_table(taxo.ctl)[j,i-1])
      }
    }
  }
}
#put them back together
ps.ctl1 = phyloseq(comm.ctl, taxo.ctl, meta.ctl)
ps.ctl1
#%Neonic ####
comm.neo = otu_table(ps.neo)
taxo.neo = tax_table(ps.neo)
meta.neo = sample_data(ps.neo)
#replace unassigned by the lowest  assigned taxonomic level
for(i in c(1:length(colnames(tax_table(taxo.neo))))){
  for(j in rownames(tax_table(taxo.neo))){
    if(is.na(tax_table(taxo.neo)[j,i]) == TRUE | tax_table(taxo.neo)[j,i] == "Unknown_Family"){
      if(substr(tax_table(taxo.neo)[j,i-1], start=1,stop=2)=="UA"){
        tax_table(taxo.neo)[j,i] <- tax_table(taxo.neo)[j,i-1]
      }
      else{
        tax_table(taxo.neo)[j,i] <- paste0("UA_",tax_table(taxo.neo)[j,i-1])
      }
    }
  }
}
#put them back together
ps.neo1 = phyloseq(comm.neo, taxo.neo, meta.neo)

#Build networks ####
eco_network_2<-function(sparcc.boot,phylo, threshold, corThreshold=0.3, sign="both"){
  
  #do a first run
  sparcc <- sparcc(otu_table(phylo))
  #calculate p-values
  sparcc.pvals <- pval.sparccboot(sparcc.boot)
  pvals.dat <- data.frame(sparcc.pvals$cors, sparcc.pvals$pvals) 
  #reorganize the p-value matrix
  cors <- sparcc.pvals$cors
  pvals <- sparcc.pvals$pvals
  sparCCpcors <- diag(0.5, nrow = dim(sparcc$Cor)[1], ncol = dim(sparcc$Cor)[1])
  sparCCpcors[upper.tri(sparCCpcors, diag=FALSE)] <- cors
  sparCCpcors <- sparCCpcors + t(sparCCpcors)
  
  sparCCpval <- diag(0.5, nrow = dim(sparcc$Cor)[1], ncol = dim(sparcc$Cor)[1])
  sparCCpval[upper.tri(sparCCpval, diag=FALSE)] <- pvals
  sparCCpval <- sparCCpval + t(sparCCpval)
  
  rownames(sparCCpcors) <- colnames(otu_table(phylo))
  colnames(sparCCpcors) <- colnames(otu_table(phylo))
  rownames(sparCCpval) <- colnames(otu_table(phylo))
  colnames(sparCCpval) <- colnames(otu_table(phylo))
  
  sparCCpcors[1:5, 1:5]
  sparCCpval[1:5, 1:5]
  
  # Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  
  reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }
  
  reorder_cor_and_p <- function(cormat, pmat){
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
    pmat <- pmat[hc$order, hc$order]
    list(r = cormat, p = pmat)
  }
  
  reordered_all_sparcc <- reorder_cor_and_p(sparCCpcors, sparCCpval)
  reordered_sparccCor <- reordered_all_sparcc$r
  reordered_sparccP<- reordered_all_sparcc$p
  
  
  
  sparccCor_processed <- reordered_sparccCor  %>% get_upper_tri() %>% reshape2::melt() %>% na.omit() %>% 
    dplyr::rename(cor = value)
  sparccP_processed <- reordered_sparccP  %>% get_upper_tri() %>% reshape2::melt() %>% na.omit() %>% 
    dplyr::rename(p = value)
  
  SparccP <- left_join(sparccCor_processed, sparccP_processed, by = c("Var1", "Var2")) %>%
    
    mutate(fdr = p.adjust(p, method = "BH"))
  
  #remove NAs
  SparccP$p[is.na(SparccP$p)] <- 1
  SparccP$fdr[is.na(SparccP$fdr)] <- 1
  
  
  sparcc.high.corr <- SparccP%>% filter(abs(cor) > corThreshold)
  sparccOkP <- sparcc.high.corr%>% filter(fdr < threshold)
  #keep only positive interactions
  sparcc.pos <- sparcc.high.corr%>% filter(cor > 0.00)
  sparcc.neg <- sparcc.high.corr%>% filter(cor < 0.00)
  
  print(paste(dim(sparccOkP)[1]," correlations with adjusted p-values under your", " X ", "threshold were found!", sep=""))
  
  #generate a matrix of fdr-adjusted p-values for igraph
  if (sign == "both"){
    mat.sparccP <- reshape2::acast(sparcc.high.corr, Var1~Var2 )
    sparcc.pos$sign <- "pos"
    sparcc.neg$sign <- "neg"
    list <- rbind(sparcc.pos,sparcc.neg)
  }
  else if (sign == "positive"){
    mat.sparccP <- reshape2::acast(sparcc.pos, Var1~Var2 )
    list <- sparcc.pos
  }
  else if (sign == "negative"){
    mat.sparccP <- reshape2::acast(sparcc.neg, Var1~Var2 )
    list <- sparcc.neg
  }
  else {
    return('Wrong sign argument, please choose between "positive", "negative" and "both"')
  }
  
  mat.sparccP[lower.tri(mat.sparccP)] <- t(mat.sparccP)[lower.tri(mat.sparccP)]
  mat.sparccP[is.na(mat.sparccP)] <- 1
  
  SparccP_plot <- sparcc.high.corr %>% ggplot(aes(x = Var2, y = Var1, fill = cor)) + geom_tile() + scale_fill_gradient2() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_point(data = sparccOkP, shape = 1)
  
  sparcc.graph <- abs(mat.sparccP) < threshold
  print(paste("So a matrix with ",table(sparcc.graph > 0)[2]/2," non-zero relationships was made",sep=""))
  diag(sparcc.graph) <- 0
  #keep a vector of which taxa have more than one relationship over the threshhold
  posit.vector <- apply(sparcc.graph,2,sum) > 0
  table(posit.vector)
  
  print(paste("Out of the ",length(posit.vector)," taxa in the analysis, ",length(posit.vector)-table(posit.vector)["TRUE"]," did not have any significant relationship with another taxa",sep=""))
  #keep a sparcc.graph version of non-zero relationships
  sparcc.graph.filt <- sparcc.graph[posit.vector,posit.vector]
  #edit our phyloseq object to keep ASVs with a least one significant relationship
  vect.filt <- rownames(tax_table(phylo)) %in% rownames(sparcc.graph.filt)
  
  #make a new phyloseq object
  OTU <- otu_table(phylo)[,vect.filt]
  TAX <- tax_table(phylo)[vect.filt,]
  DATA <- sample_data(phylo)
  phylo.filt <- phyloseq(OTU,TAX,DATA)
  print(paste("So ",length(rownames(tax_table(phylo.filt)))," out of the ",length(rownames(tax_table(phylo)))," taxa were kept in the final phyloseq object",sep=""))
  sparcc.graph.filt <- Matrix(sparcc.graph.filt, sparse=TRUE)
  
  return(list(SparccP_plot, sparcc.graph.filt,phylo.filt,list))
  
}

#BOTH####
eco.net.ctl.both = eco_network_2(sparcc.boot.ctl1000, ps.ctl1, threshold = 0.01, corThreshold = 0.3, sign="both")
eco.net.neo.both = eco_network_2(sparcc.boot.neo1000, ps.neo1, threshold = 0.01, corThreshold = 0.3, sign="both")
#Correlations ####
eco.net.ctl.both[1] 
eco.net.neo.both[1]

#%Control ####
ps.filt.ctl <- eco.net.ctl.both[3][[1]]
filt.tax.ctl <- as.data.frame(tax_table(ps.filt.ctl))
levels(factor(filt.tax.ctl$Phylum)) #not all of them pass the fdr < 0.01 threshold
filt.tax.ctl.bac = filt.tax.ctl[which(filt.tax.ctl$Phylum != "Nematoda"),]
levels(factor(filt.tax.ctl.bac$Family))
filt.tax.ctl.nem = filt.tax.ctl[which(filt.tax.ctl$Phylum == "Nematoda"),]
levels(factor(filt.tax.ctl.nem$Family))
list.ctl = eco.net.ctl.both[4][[1]]
length(which(list.ctl$fdr < 0.01))
length(which(list.ctl$fdr < 0.01 & list.ctl$sign == "pos"))
length(which(list.ctl$fdr < 0.01 & list.ctl$sign == "neg"))
#take out network results
table_cor.ctl <- list.ctl 
#choose desired taxonomic rank
taxo.ctl1 = tax_table(ps.ctl1)
table_cor.ctl$taxa1 <- taxo.ctl1[,5][match(table_cor.ctl$Var1, rownames(taxo.ctl1),nomatch = NA)]
table_cor.ctl$taxa2 <- taxo.ctl1[,5][match(table_cor.ctl$Var2, rownames(taxo.ctl1),nomatch = NA)]
table_cor.ctl

table_cor.ctl01 = table_cor.ctl[which(table_cor.ctl$fdr < 0.01),]
table_cor.ctl01.pos = table_cor.ctl01[table_cor.ctl01$sign == "pos",]
table_cor.ctl01.pos_18s = table_cor.ctl01.pos[c(grep("18S_",table_cor.ctl01.pos$Var1),grep("18S_",table_cor.ctl01.pos$Var2)),]
#Table 2-1####
table_cor.ctl01.pos_18s

table_cor.ctl01.neg = table_cor.ctl01[table_cor.ctl01$sign == "neg",]
table_cor.ctl01.neg_18s = table_cor.ctl01.neg[c(grep("18S_",table_cor.ctl01.neg$Var1),grep("18S_",table_cor.ctl01.neg$Var2)),]
#Table 2-2####
table_cor.ctl01.neg_18s

#%Neonic ####
ps.filt.neo <- eco.net.neo.both[3][[1]]
filt.tax.neo <- as.data.frame(tax_table(ps.filt.neo))
levels(factor(filt.tax.neo$Phylum)) #not all of them pass the fdr < 0.01 threshold
levels(factor(filt.tax.neo[which(filt.tax.neo$Phylum == "Nematoda"),][,5]))
list.neo = eco.net.neo.both[4][[1]]
length(which(list.neo$fdr < 0.01))
length(which(list.neo$fdr < 0.01 & list.neo$sign == "pos"))
length(which(list.neo$fdr < 0.01 & list.neo$sign == "neg"))
#take out network results
table_cor.neo <- list.neo 
#choose desired taxonomic rank
taxo.neo1 = tax_table(ps.neo1)
table_cor.neo$taxa1 <- taxo.neo1[,5][match(table_cor.neo$Var1, rownames(taxo.neo1),nomatch = NA)]
table_cor.neo$taxa2 <- taxo.neo1[,5][match(table_cor.neo$Var2, rownames(taxo.neo1),nomatch = NA)]
table_cor.neo
#taxo.neo1[which(rownames(taxo.neo1) == "16S_ASV14628")] #verification

table_cor.neo01 = table_cor.neo[which(table_cor.neo$fdr < 0.01),]
table_cor.neo01.pos = table_cor.neo01[table_cor.neo01$sign == "pos",]
table_cor.neo01.pos_18s = table_cor.neo01.pos[c(grep("18S_",table_cor.neo01.pos$Var1),grep("18S_",table_cor.neo01.pos$Var2)),]
#Table 2-3####
table_cor.neo01.pos_18s

table_cor.neo01.neg = table_cor.neo01[table_cor.neo01$sign == "neg",]
table_cor.neo01.neg_18s = table_cor.neo01.neg[c(grep("18S_",table_cor.neo01.neg$Var1),grep("18S_",table_cor.neo01.neg$Var2)),]
#Table 2-4####
table_cor.neo01.neg_18s

#Network plots ####
#%---Control ####
sparcc.graph.filt.ctl = eco.net.ctl.both[2][[1]]
#igraph 
ig.sparcc.ctl = adj2igraph(sparcc.graph.filt.ctl, vertex.attr=list(name=colnames(sparcc.graph.filt.ctl)))
plot.igraph(ig.sparcc.ctl) #check the nodes (compare randomly with View(table_cor.ctl01))
vsize.ctl = colMeans(clr(otu_table(ps.filt.ctl)[,colnames(sparcc.graph.filt.ctl)],1)) #for ggnetwork (should give ASVs)
#attributes #### 
length(E(ig.sparcc.ctl)) #number of edges
length(V(ig.sparcc.ctl)) #number of nodes
vertex_attr(ig.sparcc.ctl, "Kingdom" ) <- as.character(tax_table(ps.filt.ctl)[,1])
vertex_attr(ig.sparcc.ctl, "Phylum" ) <- as.character(tax_table(ps.filt.ctl)[,2])
vertex_attr(ig.sparcc.ctl, "Family" ) <- as.character(tax_table(ps.filt.ctl)[,5])
vertex_attr(ig.sparcc.ctl, "vsize" ) <- vsize.ctl
#edge_list
edge_list.ctl <- data.frame(as_edgelist(ig.sparcc.ctl)) 
#put the correlation table in the same order as the list of  edges
rownames(eco.net.ctl.both[[4]]) <- paste(eco.net.ctl.both[[4]]$Var1,eco.net.ctl.both[[4]]$Var2)
rownames(edge_list.ctl) <- paste(edge_list.ctl$X1,edge_list.ctl$X2)
list.ctl <- eco.net.ctl.both[[4]][rownames(edge_list.ctl),]
#add the sign as an edge attribute 
edge_attr(ig.sparcc.ctl, "sign") <- list.ctl$sign
levels(factor(as.character(tax_table(ps.filt.ctl)[,2]))); length(levels(factor(as.character(tax_table(ps.filt.ctl)[,2]))))
levels(factor(as.character(tax_table(ps.filt.ctl)[,5]))); length(levels(factor(as.character(tax_table(ps.filt.ctl)[,5]))))
length(grep("18S_",rownames(tax_table(ps.filt.ctl)))) #nematode families
length(grep("16S_",rownames(tax_table(ps.filt.ctl)))) #bacteria families
#data frame
ig.sparcc.ctl.df=as.data.frame(vertex_attr(ig.sparcc.ctl))
#bacteria phyla
levels(factor(ig.sparcc.ctl.df[grep("Bacteria",ig.sparcc.ctl.df$Kingdom),][,3]))
#nematodes families
levels(factor(ig.sparcc.ctl.df[grep("Nematoda",ig.sparcc.ctl.df$Phylum),][,4]))

#%plot_network ####
plot_network(ig.sparcc.ctl, ps.filt.ctl, point_size = vsize.ctl, alpha = 0.1,
             type='taxa', color="Phylum", shape = "Kingdom", line_color = "Phylum", line_alpha = 0.5)

#convert to a network object
#Table 3-1 ###
library(intergraph);packageVersion("intergraph") #2.0.2
ne.sparcc.ctl <- asNetwork(x = ig.sparcc.ctl); ne.sparcc.ctl
ggnet.ctl <- ggnetwork(ne.sparcc.ctl, layout = "fruchtermanreingold")
#FIG. 3A ####
#replace Solibacteraceae_(Subgroup_3) w/ Solibacteraceae
ggnet.ctl$Family[which(ggnet.ctl$Family == "Solibacteraceae_(Subgroup_3)")] = "Solibacteraceae"
plot.ctl.both <- ggplot(ggnet.ctl, aes(x = x, y = y, xend = xend, yend = yend)) + 
  theme_blank() + 
  geom_edges(aes(linetype = sign), color = "antiquewhite4",size=0.3,alpha = 0.5) + 
  scale_linetype_manual(labels=c("Negative","Positive"), values=c(2,1)) + 
  labs(linetype="Type of correlation") + 
  geom_nodes(aes(color = Phylum, size = vsize,shape = Kingdom),alpha = 0.8) + 
  scale_size_continuous(range=c(2,12)) +
  scale_color_manual("Phylum",values = c("yellowgreen","tomato3","darkolivegreen4","indianred1","aquamarine2",
                                         "darkgoldenrod1","cyan3","red2","darkorange2","chartreuse3",
                                         "mediumvioletred","seagreen","antiquewhite4","cornflowerblue","tan4")) +
  geom_nodetext_repel(aes(label = Family),max.overlaps = Inf,cex = 2,color="black",fontface = "bold",
                      segment.size = 0.02, segment.colour = "gray") + #, box.padding = unit(1, "lines")
  scale_shape_manual("Microbial Community", labels = c("Bacteria"="Bacteria", "Metazoa"="Nematodes"), values = c(19,17)) +
  labs(tag = "A)") + theme(plot.tag = element_text(size = 18, face = "bold")) +
  theme(legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14)) + #title = element_text(size=22,face="bold")
  guides(colour = guide_legend(override.aes = list(size=6)),
         shape = guide_legend(override.aes = list(size=6)),
         size = FALSE) 
plot.ctl.both

#%---Neonic ####
sparcc.graph.filt.neo = eco.net.neo.both[2][[1]]
#igraph 
ig.sparcc.neo = adj2igraph(sparcc.graph.filt.neo, vertex.attr=list(name=colnames(sparcc.graph.filt.neo)))
plot.igraph(ig.sparcc.neo) #check the nodes (compare randomly with View(table_cor.neo01))
vsize.neo = colMeans(clr(otu_table(ps.filt.neo)[,colnames(sparcc.graph.filt.neo)],1)) #for ggnetwork (should give ASVs)

tax_new = tax_table(ps.filt.neo)[match(colnames(sparcc.graph.filt.neo), rownames(tax_table(ps.filt.neo)),nomatch = NA)]

#attributes #### 
length(E(ig.sparcc.neo)) #number of edges
length(V(ig.sparcc.neo)) #number of nodes
vertex_attr(ig.sparcc.neo, "Kingdom" ) <- as.character(tax_new[,1])
vertex_attr(ig.sparcc.neo, "Phylum" ) <- as.character(tax_new[,2])
vertex_attr(ig.sparcc.neo, "Family" ) <- as.character(tax_new[,5])
vertex_attr(ig.sparcc.neo, "vsize" ) <- vsize.neo
#edge_list
edge_list.neo <- data.frame(as_edgelist(ig.sparcc.neo)) 
#put the correlation table in the same order as the list of  edges
rownames(eco.net.neo.both[[4]]) <- paste(eco.net.neo.both[[4]]$Var1,eco.net.neo.both[[4]]$Var2)
rownames(edge_list.neo) <- paste(edge_list.neo$X1,edge_list.neo$X2)
list.neo <- eco.net.neo.both[[4]][rownames(edge_list.neo),]
#add the sign as a edge attribute 
edge_attr(ig.sparcc.neo, "sign") <- list.neo$sign
levels(factor(as.character(tax_table(ps.filt.neo)[,2]))); length(levels(factor(as.character(tax_table(ps.filt.neo)[,2]))))
levels(factor(as.character(tax_table(ps.filt.neo)[,5]))); length(levels(factor(as.character(tax_table(ps.filt.neo)[,5]))))
length(grep("18S_",rownames(tax_table(ps.filt.neo)))) #nematode families
length(grep("16S_",rownames(tax_table(ps.filt.neo)))) #bacteria families
#data frame
ig.sparcc.neo.df=as.data.frame(vertex_attr(ig.sparcc.neo))
#bacteria phyla
levels(factor(ig.sparcc.neo.df[grep("Bacteria",ig.sparcc.neo.df$Kingdom),][,3]))
#nematodes families
levels(factor(ig.sparcc.neo.df[grep("Nematoda",ig.sparcc.neo.df$Phylum),][,4]))

#%plot_network ####
plot_network(ig.sparcc.neo, ps.filt.neo, point_size = (vsize.neo*3), alpha = 0.8,
             type='taxa', label = "Family",color="Phylum", shape = "Kingdom", line_color = "Phylum")

#convert to a network object
#Table 3-2 ###
ne.sparcc.neo <- asNetwork(x = ig.sparcc.neo); ne.sparcc.neo
ggnet.neo <- ggnetwork(ne.sparcc.neo, layout = "fruchtermanreingold")

#FIG. 3B ####
#replace Solibacteraceae_(Subgroup_3) w/ Solibacteraceae
ggnet.neo$Family[which(ggnet.neo$Family == "Solibacteraceae_(Subgroup_3)")] = "Solibacteraceae"
plot.neo.both <- ggplot(ggnet.neo, aes(x = x, y = y, xend = xend, yend = yend)) + 
  theme_blank() + 
  geom_edges(aes(linetype = sign), color = "antiquewhite4",size=0.3,alpha = 0.5) + 
  scale_linetype_manual(labels=c("Negative","Positive"), values=c(2,1)) + 
  labs(linetype="Type of correlation") + 
  geom_nodes(aes(color = Phylum, size = vsize,shape = Kingdom),alpha = 0.8) + 
  scale_size_continuous(range=c(2,12)) +
  scale_color_manual("Phylum",values = c("yellowgreen","tomato3","darkolivegreen4","indianred1","darkgoldenrod1",
                                         "cyan2","red2","darkorange2","chartreuse3","mediumvioletred",
                                         "seagreen","antiquewhite4","cornflowerblue","tan4")) +
  geom_nodetext_repel(aes(label = Family),max.overlaps = Inf,cex = 2,color="black",fontface = "bold",
                      segment.size = 0.02, segment.colour = "gray") + #, box.padding = unit(1, "lines")
  scale_shape_manual("Microbial Community", labels = c("Bacteria"="Bacteria", "Metazoa"="Nematodes"), values = c(19,17)) +
  labs(tag = "B)") + theme(plot.tag = element_text(size = 18, face = "bold")) +
  theme(legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14)) + #title = element_text(size=22,face="bold")
  guides(colour = guide_legend(override.aes = list(size=6)),
         shape = guide_legend(override.aes = list(size=6)),
         size = FALSE) 
plot.neo.both

#NULL MODELS ####
#Networks ####
net_null <- function (network, runs = 999)
{
  # actually randomize edges
  network.adj <- as_adjacency_matrix(network)
  network.adj.d <- as.dist(network.adj)
  network.adj.m <- as.matrix(network.adj.d)
  
  ###list of metrics###
  metrics.list <- c("diameter","motifs","edge_connectivity","vertex_connectivity","girth","radius","centralized_betweenness","centralized_degree","centralized_eigen_centrality","modularity")
  
  # diameter
  diameter.obs <- diameter(network)
  diameter.rnd <- vector(length=runs)
  
  # number of motifs
  motifs.obs <- count_motifs(network)
  motifs.rnd <- vector(length=runs)
  
  # edge_connectivity
  edge_connectivity.obs <- edge_connectivity(network)
  edge_connectivity.rnd <- vector(length=runs)
  
  # vertex_connectivity
  vertex_connectivity.obs <- vertex_connectivity(network)
  vertex_connectivity.rnd <- vector(length=runs)
  
  # girth
  girth.obs <- girth(network)[[1]]
  girth.rnd <- vector(length=runs)
  
  # radius
  radius.obs <- radius(network)
  radius.rnd <- vector(length=runs)
  
  # centralized_betweenness
  centralized_betweenness.obs <- centr_betw(network)$centralization
  centralized_betweenness.rnd <- vector(length=runs)
  
  # # centralized_closeness
  # centralized_closeness.obs <- centr_clo(network)$centralization
  # centralized_closeness.rnd <- vector(length=runs)
  
  # centralized_degree
  centralized_degree.obs <- centr_degree(network)$centralization
  centralized_degree.rnd <- vector(length=runs)
  
  # centralized_eigen_centrality
  centralized_eigen_centrality.obs <- centr_eigen(network)$centralization
  centralized_eigen_centrality.rnd <- vector(length=runs)
  
  # # number of edges
  # edge_number.obs <- length(E(network)) 
  # edge_number.rnd <- vector(length=runs)
  
  # modularity
  wtc.obs <- cluster_walktrap(network)
  modularity.obs <- modularity(network, membership(wtc.obs))
  modularity.rnd <- vector(length=runs)
  
  for (i in 1:runs) 
  {
    network.adj.m.shuf <- network.adj.m
    network.adj.m.shuf[lower.tri(network.adj.m.shuf)]  <- sample(network.adj.d)
    network.adj.d.shuf <- as.dist(network.adj.m.shuf)
    network.shuf <- graph_from_adjacency_matrix(network.adj.d.shuf, mode="undirected")
    
    # to change metric, change this line
    #metric.rnd[i] <- diameter(network.shuf)
    
    # diameter
    diameter.rnd[i] <- diameter(network.shuf)
    
    # number of motifs
    motifs.rnd[i] <- count_motifs(network.shuf)
    
    # edge_connectivity
    edge_connectivity.rnd[i] <- edge_connectivity(network.shuf)
    
    # vertex_connectivity
    vertex_connectivity.rnd[i] <- vertex_connectivity(network.shuf)
    
    # girth
    girth.rnd[i] <- girth(network.shuf)[[1]]
    
    # radius
    radius.rnd[i] <- radius(network.shuf)
    
    # centralized_betweenness
    centralized_betweenness.rnd[i] <- centr_betw(network.shuf)$centralization
    
    # # centralized_closeness
    # centralized_closeness.rnd[i] <- centr_clo(network.shuf)$centralization
    
    # centralized_degree
    centralized_degree.rnd[i] <- centr_degree(network.shuf)$centralization
    
    # centralized_eigen_centrality
    centralized_eigen_centrality.rnd[i] <- centr_eigen(network.shuf)$centralization
    
    # # edge_number
    # edge_number.rnd[i] <- length(E(network))
    
    #modularity
    modularity.rnd[i] <- modularity(network.shuf,membership(cluster_walktrap(network.shuf)))
  }
  
  #iterate over each metric and calculate statistics
  results <- data.frame()
  distribution <- data.frame(matrix(nrow=runs,ncol=length(metrics.list)))
  rownames(distribution) <- c(1:runs)
  colnames(distribution) <- metrics.list
  
  
  for (metric in metrics.list) {
    
    metric.obs <- get(paste0(metric,".obs"))
    
    metric.rnd <- get(paste0(metric,".rnd"))
    distribution[,metric] <- metric.rnd
    metric.rand.mean <- mean(metric.rnd)
    metric.rand.sd <- sd(metric.rnd)
    metric.obs.z <- (metric.obs - metric.rand.mean) / metric.rand.sd
    metric.obs.rank <- rank(c(metric.obs,metric.rnd))[1]
    metric.obs.p <- metric.obs.rank / (runs + 1)
    
    stats <- c(metric.obs, metric.rand.mean, metric.rand.sd, metric.obs.z, metric.obs.rank, metric.obs.p)
    
    results <- rbind(results,stats)
    
  }
  rownames(results) <- metrics.list
  colnames(results) <- c("Observed metric", "Null model mean", "Null model standard deviation", "Observed metric z score", "Observed metric rank", "Observed metric p-value")
  
  return(list(results,distribution))
  
}

#%---Control ####
net.null_ctl <- net_null(ig.sparcc.ctl); net.null_ctl
#extract the variables
net.null.stats_ctl <- net.null_ctl[[1]]
net.null.dist_ctl <- data.frame(net.null_ctl[[2]][,c("motifs","modularity")])
net.null.dist.melt_ctl <- reshape2::melt(net.null.dist_ctl)
net.null.dist.melt_ctl$variable <- as.character(net.null.dist.melt_ctl$variable)
for (row in rownames(net.null.dist.melt_ctl)){
  metric <- net.null.dist.melt_ctl[row,"variable"]
  net.null.dist.melt_ctl[row,"Obs"] <- net.null.stats_ctl[metric,"Observed metric"]
}
net.null.dist.melt_ctl$variable <- as.factor(net.null.dist.melt_ctl$variable)
#Table 3-3####
rbind(net.null.stats_ctl[1:2,],net.null.stats_ctl[7,],net.null.stats_ctl[9:10,])

#%---Neonic ####
net.null_neo <- net_null(ig.sparcc.neo); net.null_neo
#extract the variables
net.null.stats_neo <- net.null_neo[[1]]
net.null.dist_neo <- data.frame(net.null_neo[[2]][,c("motifs","modularity")])
net.null.dist.melt_neo <- reshape2::melt(net.null.dist_neo)
net.null.dist.melt_neo$variable <- as.character(net.null.dist.melt_neo$variable)
for (row in rownames(net.null.dist.melt_neo)){
  metric <- net.null.dist.melt_neo[row,"variable"]
  net.null.dist.melt_neo[row,"Obs"] <- net.null.stats_neo[metric,"Observed metric"]
}
net.null.dist.melt_neo$variable <- as.factor(net.null.dist.melt_neo$variable)
#Table 3-4####
rbind(net.null.stats_neo[1:2,],net.null.stats_neo[7,],net.null.stats_neo[9:10,])

#Compare to networks ####
#combine two treatment
net.null.dist.melt_ctl$trt <- "Control"
net.null.dist.melt_neo$trt <- "Neonicotinoid-treated"
net.null.dist.melt_ctl$trtx <- 1
net.null.dist.melt_neo$trtx <- 2
net.null.dist.melt_trt <- rbind(net.null.dist.melt_ctl,net.null.dist.melt_neo)

#violin graphs
ctl.pos.v <- ggplot(net.null.dist.melt_ctl, aes(x=1, y=value)) + 
  geom_violin() + facet_wrap(net.null.dist.melt_ctl$variable,scales = "free_y") + 
  geom_hline(aes(yintercept=Obs, color=variable), size=3) + theme(legend.position = "none") 
ctl.pos.v
#library(ggbeeswarm)
neo.pos.v <- ggplot(net.null.dist.melt_neo, aes(x=1, y=value)) + 
  geom_violin() + facet_wrap(net.null.dist.melt_neo$variable,scales = "free_y") + 
  geom_hline(aes(yintercept=Obs, color=variable), size=3) + theme(legend.position = "none") 
neo.pos.v 
#Fig. S3 ####
trt.pos.v <- ggplot(net.null.dist.melt_trt, aes(x=trt, y=value)) + 
  theme_bw() +
  geom_violin() + 
  facet_wrap(~ variable,scales = "free_y",
             labeller = labeller(variable = c(motifs = "Motifs", modularity = "Modularity"))) + 
  geom_segment(aes(x=trtx-0.5,xend=trtx+0.5,y=Obs,yend=Obs,color=trt),size=2) +
  scale_color_manual(values = c("cornflowerblue","mediumvioletred")) +
  xlab("\nTreatments") +
  ylab("Values") +
  theme(legend.position = "none", 
        axis.title = element_text(size=16,face="bold"),
        strip.text = element_text(size=16,face="bold"),
        axis.text.x = element_text(size = 14,face="bold"),
        axis.text.y = element_text(size = 14),
        strip.background = element_rect(color = "white", size = 1))
trt.pos.v

#Network vertex ####
net_vertex_null <- function (network, runs = 999, model="unrestricted")
{
  # actually randomize edges
  network.adj <- as_adjacency_matrix(network)
  network.adj.d <- as.dist(network.adj)
  network.adj.m <- as.matrix(network.adj.d)
  
  ###extract vertex names###
  if (!is.null(network.adj@Dimnames[[1]])){
    vertex.names <- network.adj@Dimnames[[1]]
  }
  else{
    vertex.names <- c(1:dim(network.adj)[1])
  }
  
  
  ###list of metrics###
  metrics.list <- c("authority_score","closeness","coreness","degree","estimated_betweenness","average_nearest_neighbor_degree","strength","betweenness_centrality","closeness_centrality","degree_centrality","eigen_centrality")
  
  #authority_score
  authority_score.obs <- authority_score(network)$vector 
  authority_score.rnd <- data.frame(matrix(0,nrow = length(vertex.names)))
  
  
  #closeness
  closeness.obs <- closeness(network)
  closeness.rnd <- data.frame(matrix(0,nrow = length(vertex.names)))
  
  
  #coreness
  coreness.obs <- coreness(network)
  coreness.rnd <- data.frame(matrix(0,nrow = length(vertex.names)))
  
  #degree
  degree.obs <- degree(network) 
  degree.rnd <- data.frame(matrix(0,nrow = length(vertex.names)))
  
  
  
  #estimated_betweenness
  estimated_betweenness.obs <- estimate_betweenness(network,cutoff = 0) 
  estimated_betweenness.rnd <- data.frame(matrix(0,nrow = length(vertex.names)))
  
  
  #average_nearest_neighbor_degree
  average_nearest_neighbor_degree.obs <- knn(network)$knn 
  average_nearest_neighbor_degree.rnd <- data.frame(matrix(0,nrow = length(vertex.names)))
  
  
  #strength
  strength.obs <- strength(network) 
  strength.rnd <- data.frame(matrix(0,nrow = length(vertex.names)))
  
  
  #betweenness_centrality
  betweenness_centrality.obs <- centr_betw(network)[[1]]
  betweenness_centrality.rnd <- data.frame(matrix(0,nrow = length(vertex.names)))
  names(betweenness_centrality.obs) <- vertex.names
  
  
  #closeness_centrality
  closeness_centrality.obs <- centr_clo(network)[[1]]
  closeness_centrality.rnd <- data.frame(matrix(0,nrow = length(vertex.names)))
  names(closeness_centrality.obs) <- vertex.names
  
  
  #degree_centrality
  degree_centrality.obs <- centr_degree(network)[[1]]
  degree_centrality.rnd <- data.frame(matrix(0,nrow = length(vertex.names)))
  names(degree_centrality.obs) <- vertex.names
  
  
  #eigen_centrality
  eigen_centrality.obs <- centr_eigen(network)[[1]]
  eigen_centrality.rnd <- data.frame(matrix(0,nrow = length(vertex.names)))
  names(eigen_centrality.obs) <- vertex.names
  
  for (i in 1:runs) 
  {
    if (model=="unrestricted"){
      
      network.adj.m.shuf <- network.adj.m
      network.adj.m.shuf[lower.tri(network.adj.m.shuf)]  <- sample(network.adj.d)
      network.adj.d.shuf <- as.dist(network.adj.m.shuf)
      network.shuf <- graph_from_adjacency_matrix(network.adj.d.shuf, mode="undirected")
      vertex_attr(network.shuf)$name <- vertex_attr(network)$name
      
    }
    
    else if (model=="restricted"){
      
      network.dim <- dim(network.adj.m)[1]
      network.label.shuf <- sample(network.dim)
      network.adj.m.shuf <- network.adj.m[network.label.shuf, network.label.shuf]
      network.shuf <- graph_from_adjacency_matrix(network.adj.m.shuf, mode="undirected")
      vertex_attr(network.shuf)$name <- vertex_attr(network)$name
    }
    
    
    
    #authority_score
    
    authority_score.rnd[,i] <- authority_score(network.shuf)$vector 
    
    #closeness
    closeness.rnd[,i] <- closeness(network.shuf)
    
    #coreness
    coreness.rnd[,i] <- igraph::coreness(network.shuf)
    
    #degree
    degree.rnd[,i] <- igraph::degree(network.shuf)
    
    #estimated_betweenness
    estimated_betweenness.rnd[,i] <- estimate_betweenness(network.shuf,cutoff = 0)
    
    #average_nearest_neighbor_degree
    average_nearest_neighbor_degree.rnd[,i] <- knn(network.shuf)$knn 
    
    #strength
    strength.rnd[,i] <- strength(network.shuf) 
    
    #betweenness_centrality
    rnd.cb <- igraph::centr_betw(network.shuf)[[1]]
    names(rnd.cb) <- vertex.names
    betweenness_centrality.rnd[,i] <- rnd.cb
    
    #closeness_centrality
    rnd.cc <- centr_clo(network.shuf)[[1]]
    names(rnd.cc) <- vertex.names
    closeness_centrality.rnd[,i] <- rnd.cc
    
    #degree_centrality
    rnd.cd <- centr_degree(network.shuf)[[1]]
    names(rnd.cd) <- vertex.names
    degree_centrality.rnd[,i] <- rnd.cd
    
    #eigen_centrality
    rnd.ce <- centr_eigen(network.shuf)[[1]]
    
    names(rnd.ce) <- vertex.names
    
    eigen_centrality.rnd[,i] <- rnd.ce
    
  }
  
  
  results.list <- vector("list", 11)
  #iterate over each metric and calculate statistics
  for (n in c(1:length(metrics.list))) {
    metric <- metrics.list[n]
    results <- data.frame()
    
    metric.dat <- as.data.frame(get(paste0(metric,".rnd")))
    rownames(metric.dat) <- vertex.names
    
    
    
    for (vertex in vertex.names){
      
      metric.obs <- get(paste0(metric,".obs"))[vertex]
      metric.rnd <- as.matrix(metric.dat[vertex,])
      
      
      metric.rand.mean <- mean(metric.rnd)
      
      metric.rand.sd <- sd(metric.rnd)
      
      metric.obs.z <- (metric.obs - metric.rand.mean) / metric.rand.sd
      
      metric.obs.rank <- rank(c(metric.obs,metric.rnd))[1]
      metric.obs.p <- metric.obs.rank / (runs + 1)
      
      stats <- c(metric.obs, metric.rand.mean, metric.rand.sd, metric.obs.z, metric.obs.rank, metric.obs.p)
      results <- rbind(results,stats)
      
      
    }
    
    #put the vertex names
    rownames(results) <- vertex.names
    colnames(results) <- c("Observed metric", "Null model mean", "Null model standard deviation", "Observed metric z score", "Observed metric rank", "Observed metric p-value")
    
    results.list[n] <- list(results)
    
    
    
    
    
  }
  
  names(results.list) <- metrics.list
  
  return(results.list)
  
  
}

#%---Unrestricted graphs ####
vert.null.unrest.ctl <- net_vertex_null(ig.sparcc.ctl) 
vert.null.unrest.neo <- net_vertex_null(ig.sparcc.neo) 
#keep only interesting metrics
vert.null.unrest.ctl1 = vert.null.unrest.ctl[c("betweenness_centrality","coreness")]
vert.null.unrest.ctl.df = as.data.frame(vert.null.unrest.ctl[c("betweenness_centrality","coreness")])
vert.null.unrest.neo1 = vert.null.unrest.neo[c("betweenness_centrality","coreness")]
vert.null.unrest.neo.df = as.data.frame(vert.null.unrest.neo[c("betweenness_centrality","coreness")])

#Table S1 ####
#Based on ranking of Tables S4A and S4B
#Betweenness Centrality 
View(vert.null.unrest.ctl1$betweenness_centrality[order(vert.null.unrest.ctl1$betweenness_centrality$`Observed metric`, decreasing = TRUE),])
View(vert.null.unrest.neo1$betweenness_centrality[order(vert.null.unrest.neo1$betweenness_centrality$`Observed metric`, decreasing = TRUE),])
#Coreness
View(vert.null.unrest.ctl1$coreness[order(vert.null.unrest.ctl1$coreness$`Observed metric`, decreasing = TRUE),])
View(vert.null.unrest.neo1$coreness[order(vert.null.unrest.neo1$coreness$`Observed metric`, decreasing = TRUE),])


