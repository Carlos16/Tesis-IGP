
#Read data
bodysizes<-read.delim('bodysizes_2008.txt')
#explore structure
str(bodysizes)

#subset
benguela <- subset(bodysizes,Geographic.location == "Africa, Benguela ecosystem")

write.csv(benguela,'benguela.csv',row.names = FALSE)

library(igraph)

#igraph

# reordering benguela , stackoverflow

# Put the columns in icols to the start of the dataframe
ReorderColumns<-function(DF,icols){
  cols <- c(icols,names(DF)[-which(names(DF) %in% icols)])
  DF <- DF[cols]
  return(DF)
}

# adding foraging mode

SetForagingMode <- function(Prey_Pred, MovData){
    Prey <-  Prey_Pred[1]
    Pred <- Prey_Pred[2]
    M.Pred<- MovData[MovData[1] == Pred , 2]
    M.Prey <- MovData[MovData[1] == Prey , 2]
    return(ifelse(M.Pred == "A", ifelse(M.Prey == "A", "Ac", "Gr"), ifelse(M.Prey == "A", "Sw", NA)))
}

Mov.data <- read.csv("MovdataBenguela.csv")
icols <- c("Common.name.s..resource","Common.name.s..consumer")
benguela<-read.csv('benguela.csv',as.is = T)
benguela<-ReorderColumns(benguela,icols)
FM <- apply(benguela[,1:2], 1, SetForagingMode , Mov.data)
benguela$FM <- FM

G <- graph.data.frame(benguela)



E(G)$color <- ifelse(E(G)$FM == "Ac", "darkred", "darkgreen")
V(G)$size <- degree(G)/2.5

layout <- layout.fruchterman.reingold(G, circular = T)
plot(G, layout = layout ,  vertex.label.cex = 0.7 , vertex.label.color = "black" , edge.arrow.size = 0.3, vertex.color = "lightblue")
pdf("Benguela.pdf")
dev.off()
library(ggplot2)
ggsave("Benguela.pdf",limitsize=FALSE, height = 8, width = 12)  
G2 <-  network(benguela, directed =  TRUE)

ggnet2(G2, size = 6, color = "Metabolic.category.resource")
##Function for finding IGP modules

findIGPs<-function(G){
  FocusVertices <- V(G)[degree(G,mode="out") >=2]
  IGPs <- findModules(FocusVertices,G)
  return(IGPs)
}

findModules<-function(V,G){
  Modules <- c()
  for(resource in V){
    Modules<-append(Modules,IGPmod(G,resource))
  }
  return(Modules)
}  
IGPmod<-function(G,resource){
  Predators <- neighborhood(G, 1, nodes = resource, mode ="out")
  IGPs <- c()
  n <- length(Predators[[1]])
  for(i in 2:n){
    IGPs<-append(IGPs,findIGPrey(resource, i, Predators, G))
  }
  return(IGPs)
}

findIGPrey<-function(resource, predindex,Predators,G){
  pIGPs<-c()
  NeighPred <- neighborhood(G, 1, nodes = Predators[[1]][predindex], mode="out")
  for(p in Predators[[1]][-c(1,predindex)]){
    if(p %in% NeighPred[[1]]){
      pIGPs<-append(pIGPs,list(c(resource,Predators[[1]][predindex],p)))
    }
  }
  return(pIGPs)
}


IGPs<-findIGPs(G)
##Finding size ratios for the IGP modules

findIGPsizeRatiosFM <- function(IGPs,G){
  K_RC<-c()
  K_CP<-c()
  m_P <-c()
  fm1 <- c()
  fm2 <- c()
  fm3 <- c()
  for(i in 1:length(IGPs)){
      SRFM <- getsizeRatiosFM(G,IGPs[i][[1]])
      SR <- SRFM[[1]]
      FM <- SRFM[[2]]
      K_RC[i] = SR[1]
      K_CP[i] = SR[2]
      m_P[i] = SR[3]
      fm1[i] = FM[1]
      fm2[i] = FM[3]
      fm3[i] = FM[2]
  }
  
  return(data.frame(K_RC, K_CP, m_P, fm1, fm2, fm3))
}

getsizeRatiosFM<-function(G,igpL){
  SR<-c()
  fm <- c()
  igpEdges = getIGPedges(G,igpL)
  for(edge in igpEdges[1:2]){
    E <- E(G)[edge]
    SR<-append(SR,c(1/E$Consumer.resource.body.mass.ratio))
    fm <- append(fm,c(E$FM))
  }
  fm <- c(fm,E(G)[igpEdges[3]]$FM)
  m_P =E(G)[igpEdges[2]]$Mean.mass..g..consumer
  SR[3] <- m_P
  return(list(SR,fm))
}

getIGPedges<-function(G,igpL){
  edgeInd<-c()
  for(i in 1:2){
    edgeInd<-append(edgeInd,c(igpL[i],igpL[i+1]))
  }
  edgeInd <-c(edgeInd,c(igpL[1], igpL[3]))
  edges<-get.edge.ids(G,edgeInd)
  return(edges)
}

SR <- findIGPsizeRatiosFM(IGPs,G)
Res <- names(sapply(IGPs, function(x) V(G)[x][1]))
IGPreys <- names(sapply(IGPs, function(x) V(G)[x][2]))
IGPreds <- names(sapply(IGPs, function(x) V(G)[x][3]))
SR$Res <- Res
SR$IGPreys <- IGPreys
SR$IGPreds <- IGPreds
write.csv(SR,"BenguellaIGPS.csv",row.names = FALSE)
#Exploratory plots

#equation function

lm_eqn = function(m) {
  
  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits = 3));
  
  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
  }
  
  as.character(as.expression(eq));                 
}


library(ggplot2)
reg1<-lm(K_RC~K_CP , data = SR)

plot <- ggplot(SR,aes(K_CP,K_RC)) + geom_point(aes(color=m_pInt),size=3,shape=16)+labs(x = expression("log "* K[CP]), y = expression("log "* K[RC])) + geom_smooth(method="lm",col="red") +theme_bw(base_size = 10)

plot + annotate("text", x = -6 , y = -5., label = lm_eqn(reg1), colour="black", size = 5, parse=TRUE)


bodysizes<-read.csv('bodysizes.csv')
