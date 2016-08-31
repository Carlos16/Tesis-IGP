Data <- read.csv("BenguellaIGPS.csv")
Data$m_P <-  1e-3*Data$m_P
Data$fc <-  sapply(seq(1,nrow(Data)), function(x) paste(Data$fm1[x],Data$fm2[x],Data$fm3[x],sep = "-"))
Data$fc <- as.factor(Data$fc)
levels(Data$fc ) <-  c(3,1)
library(ggplot2)
library(latex2exp)

plot <- ggplot(Data,aes(log10(K_CP),log10(K_RC))) + geom_point(aes(color=log10(m_P)),size=3,shape = 16)+labs(x = expression("log "* K[CP]), y = expression("log "* K[RC])) + theme_bw(base_size = 13) + (scale_colour_gradientn(colours = rainbow(10))) + facet_grid(.~fc)
plot


#Influence of parameters
source("~/Tesis-IGP/code/Theory/Analisis/InvasibilityFunctions.r")
source("~/Tesis-IGP/code/Theory/Analisis/BasicParameters.r")

n <- 20000
limin <- -10
limsup <- 2
range <- limsup - limin
k0 <-  10**seq(limin,limsup, range/n)
D <- 3
B <- c()
B2 <- c()
B3 <- c()
phi11 = 2.0
phi21 = 2.0
phi12 = 0.2
phi22 = 0.2
phi13 = 0.02
phi23 = 0.02
for(k in k0){
    B <- c(B,GetProp(Data,D,k, e1, e2, e3, f, phi11, phi21))
    B2 <- c(B2,GetProp(Data, D, k, e1, e2, e3, f, phi12, phi22))
    B3 <- c(B3,GetProp(Data, D, k, e1, e2, e3, f, phi13, phi23))
}

nk0 <- as.numeric(sapply(k0, rep, 8))
T <- as.factor(rep(seq(1:8),length(k0)))
levels(T) <- c("InvCR","InvPR","mu3", "mu4", "InvPCR", "InvCPR", "Coex", "MutInv")
phi <- as.numeric(sapply(c(2.0, 0.2, 0.02), rep, length(B2)))
W <- data.frame(Propor = c(B,B2,B3), class = c(T,T,T), k0 = c(nk0,nk0,nk0), phi = as.factor(phi))
W$class <- as.factor(W$class)
levels(W$class) <- levels(T)


p1 <- ggplot(W,aes(x = log10(k0), y = Propor)) + geom_line(aes(color = as.factor(class)),size = 0.55) +  ylab("Proporción(p)")  + xlab(expression(log[10](kappa[0]))) + facet_grid(.~phi) + theme_bw() + scale_color_discrete( name = "Condición satisfecha", labels = list(TeX('$\\mu_1$'),TeX('$\\mu_2$'),TeX('$\\mu_3$'), TeX('$\\mu_4$'), TeX('$\\mu_1 \\vee \\mu_3$'), TeX('$\\mu_2 \\vee \\mu_4$'),TeX('$\\mu_3 \\vee \\mu_4$'), TeX('$\\mu_1 \\vee \\mu_2 \\vee \\mu_3 \\vee \\mu_4$')))
p1


ggsave("DataAna.pdf",limitsize=FALSE, height = 8, width = 12)  
