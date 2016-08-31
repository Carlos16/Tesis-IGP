source("BasicParameters.r")
source("InvasibilityFunctions.r")
source("PlottingFunc.r")


## Parameters especifications
e = 0.3
n = 700

D <- c(2,3)
phi <- c(0.1, 0.5, 1.0)
kLims <- c(-5,8)
kRange <- kLims[2]-kLims[1]
krc <- 10**seq(-5,8,kRange/n)
fm <- c("Ac", "Gr", "Sw")
# set parameter set
params <- list(D = D, phi = phi, krc = krc, fm = fm)
params <- ConstructParamSet(params)
k0 <- ifelse(params$D == 2 , 0.1, 30)
params$k0 = k0
# get mc
mc <-  InvasibilityRC(params$krc,params$D,params$k0,e,f,params$fm,params$phi)
mr <- krc*mc
# set main data frame for plotting
M <- data.frame(MassC =  mc , MassR = mr, krc = params$krc , Dim = factor(params$D), Prod = factor(params$k0), phi = factor(params$phi), ForagingMode = factor(params$fm))
levels(M$Dim) <- c("2D","3D")
levels(M$ForagingMode) <-  c("Captura Activa", "Pastoreo", "Captura Pasiva")

# plot
library(ggplot2)

Plot <-  ggplot(M,aes(x= log10(krc), y = log10(MassC),color = phi))  + geom_line() + facet_grid(Dim~ForagingMode) + xlab(expression(log[10](k[RC]))) + ylab(expression(log[10](m[C]))) + labs(colour = expression(phi)) + annotate("text", x = 0 , y = 6 , label= "frac(dC,dt) > 0",color = "black", fontface = "italic" , parse = TRUE) + theme_bw()

Plot
ggsave("R-CInv.pdf",limitsize=FALSE, height = 8, width = 8)  
