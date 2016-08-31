source("BasicParameters.r")
source("InvasibilityFunctions.r")
source("MTP.r")
source("PlottingFunc.r")


breaks <- c(-10,10)
lengths <- c(700)

krc <- 10**createSequence(breaks, lengths)
fc <- c(1)

# set parameter set
params <- list(D = D, phi = phi, krc = krc, kcp = krc, fc = fc)
params <- ConstructParamSet(params)
k0 <- ifelse(params$D == 2 , 0.1, 30)
params$k0 = k0
N <- MTPCond(params$krc, params$kcp, params$D, params$k0, e1, e2, e3, f,params$fc, params$phi)
#ConCPR <-  CondCPR(params$krc, params$kcp, params$D, params$k0, e1, e2, e3, f,params$fc, params$phi)
N <- ifelse(N < 0 , NA, N)
# set main data frame for plotting
M <- data.frame(kcp = params$kcp, krc = params$krc, Dim = factor(params$D), Prod = factor(params$k0), phi = factor(params$phi), ForagingMode = factor(params$fc), N = N)
levels(M$Dim) <- c("2D","3D")
levels(M$ForagingMode) <-  c("Gr-Gr-Ac","Ac-Sw-Sw","Ac-Ac-Ac")

### plot
library(ggplot2)

p <- ggplot(M[M$ForagingMode == "Gr-Gr-Ac",],aes(x = log10(kcp), y = log10(krc), fill = log10(N))) +facet_grid(Dim~phi) + geom_raster() + xlab(expression(log[10](k[CP]))) + ylab(expression(log[10](k[RC]))) + labs(fill = expression(log[10](omega))) + coord_cartesian(xlim = c(-10,10))  + (scale_fill_gradientn(colours = rainbow(13), na.value = "transparent")) + theme_bw()

p
ggsave("MTPvar.pdf",limitsize=FALSE, height = 8, width = 8)  
