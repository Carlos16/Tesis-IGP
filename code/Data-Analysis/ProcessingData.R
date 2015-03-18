setwd("C:\\Users\\Carlos\\Documents\\Thesis\\Tesis\\DataAnalysis")

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

icols <- c("Common.name.s..resource","Common.name.s..consumer")
benguela<-read.csv('benguela.csv')
benguela<-ReorderColumns(benguela,icols)

G <- graph.data.frame(benguela)
plot(G,label= V(G))

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
  Predators <- neighborhood(G,1,nodes = resource,mode ="out")
  IGPs <- c()
  n <- length(Predators[[1]])
  for(i in 2:n){
    IGPs<-append(IGPs,findIGPrey(resource,i,Predators,G))
  }
  return(IGPs)
}

findIGPrey<-function(resource,predindex,Predators,G){
  pIGPs<-c()
  NeighPred <- neighborhood(G,1, nodes = Predators[[1]][predindex],mode="out")
  for(p in Predators[[1]][-c(1,predindex)]){
    if(p %in% NeighPred[[1]]){
      pIGPs<-append(pIGPs,list(c(resource,Predators[[1]][predindex],p)))
    }
  }
  return(pIGPs)
}


IGPs<-findIGPs(G)
##Finding size ratios for the IGP modules

findIGPsizeRatios <- function(IGPs,G){
  K_RC<-c()
  K_CP<-c()
  m_P <-c()
  for(i in 1:length(IGPs)){
    SR <- getsizeRatios(G,IGPs[i][[1]])
    K_RC[i] = log10(SR[1])
    K_CP[i] = log10(SR[2])
    m_P[i] = log10(SR[3])
  }
  
  return(data.frame(K_RC,K_CP,m_P))
}

getsizeRatios<-function(G,igpL){
  SR<-c()
  igpEdges = getIGPedges(G,igpL)
  for(edge in igpEdges){
    E <- E(G)[edge]
    
    SR<-append(SR,c(1/E$Consumer.resource.body.mass.ratio))
  }
  
  m_P =E(G)[igpEdges[2]]$Mean.mass..g..consumer
  SR[3] <- m_P
  return(SR)
  
}

getIGPedges<-function(G,igpL){
  edgeInd<-c()
  for(i in 1:2){
    edgeInd<-append(edgeInd,c(igpL[i],igpL[i+1]))
  }
  edges<-get.edge.ids(G,edgeInd)
  return(edges)
}

SR <- findIGPsizeRatios(IGPs,G)

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
