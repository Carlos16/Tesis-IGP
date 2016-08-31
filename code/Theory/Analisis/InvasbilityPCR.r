source("BasicParameters.r")
source("InvasibilityFunctions.r")
source("PlottingFunc.r")



breaks <- c(-10,10)
lengths <- c(500)


krc <- 10**createSequence(breaks, lengths)
fc <- c(1)

# set parameter set
params <- list(D = D, phi = phi, krc = krc, kcp = krc, fc = fc)
params <- ConstructParamSet(params)
k0 <- ifelse(params$D == 2 , 0.1, 30)
params$k0 = k0
N <- InvPR(params$krc, params$kcp, params$D, params$k0, e1, e2, e3, f,params$fc, params$phi)
N2 <- InvCPR(params$krc, params$kcp, params$D, params$k0, e1, e2, e3, f,params$fc, params$phi)
N3 <- SufEstabCoex(params$krc, params$kcp, params$D, params$k0, e1, e2, e3, f,params$fc, params$phi)
#ConCPR <-  CondCPR(params$krc, params$kcp, params$D, params$k0, e1, e2, e3, f,params$fc, params$phi)

# set main data frame for plotting
M <- data.frame(kcp = params$kcp, krc = params$krc, Dim = factor(params$D), Prod = factor(params$k0), phi = factor(params$phi), ForagingMode = factor(params$fc), N = N2)
levels(M$Dim) <- c("2D","3D")
levels(M$ForagingMode) <-  c("Gr-Gr-Ac","Ac-Sw-Sw","Ac-Ac-Ac")

M$N2 = N2
M$N3 = N3
M$Cond <- ConCPR
m5 <-  getZone(1e5, M$N, M$N2, M$Cond)
M$m <-  m5
# plot
library(ggplot2)

p <- ggplot(M[M$ForagingMode == "Gr-Gr-Ac",],aes(x = log10(kcp), y = log10(krc), fill = log10(N))) +facet_grid(Dim~phi) + geom_raster() + xlab(expression(log[10](k[CP]))) + ylab(expression(log[10](k[RC]))) + labs(fill = expression(log[10](m[p]))) + coord_cartesian(xlim = c(-10,10))  + (scale_fill_gradientn(colours = rainbow(13), na.value = "transparent")) + theme_bw()

p
ggsave("InvCPRUp.pdf",limitsize=FALSE, height = 8, width = 8)  






### variable approach

M <- M[complete.cases(M),]
M$N <-  log10(M$N)

library(akima)
library(reshape)
fld <- with(M[M$phi == 2 & M$D == "2D",], interp(x = log10(kcp), y = log10(krc), z = N, duplicate = "mean"))
df <-  melt(fld$z, na.rm =TRUE)
names(df) <- c("kcp", "krc", "M")
df$kcp <- fld$x[df$kcp]
df$krc <- fld$y[df$krc]

p <- ggplot(data = df, aes( x = kcp, y = krc, z = M)) + stat_contour(aes(colour = ..level..), breaks = c(-10, -5, -3, 0, 5, 8))
p


