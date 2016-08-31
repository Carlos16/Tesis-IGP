source("BasicParameters.r")
source("InvasibilityFunctions.r")
source("PlottingFunc.r")

    
n = 100
kRange <- kLims[2]-kLims[1]
krc <- 10**seq(kLims[1],kLims[2],kRange/n)
fc <- c(1,2,3)


# set parameter set
params <- list(D = D, phi = phi, krc = krc, kcp = krc, fc = fc)
params <- ConstructParamSet(params)
k0 <- ifelse(params$D == 2 , 0.1, 30)
params$k0 = k0
N <- NecInvCPR(params$krc, params$kcp, params$D, params$k0, e1, e2, e3, f,params$fc, params$phi)
# set main data frame for plotting
M <- data.frame(kcp = params$kcp, krc = params$krc , Dim = factor(params$D), Prod = factor(params$k0), phi = factor(params$phi), ForagingMode = factor(params$fc), N = N)
levels(M$Dim) <- c("2D","3D")
levels(M$ForagingMode) <-  c("Gr-Gr-Ac","Ac-Sw-Sw","Ac-Ac-Ac")

# plot
library(ggplot2)

# data set for adding text
Dim_t = factor(c(2,2,2,3,3,3))
levels(Dim_t) = c("2D","3D")
ForagingMode_t = factor(c(1,2,3,1,2,3))
levels(ForagingMode_t) = c("Gr-Gr-Ac","Ac-Sw-Sw","Ac-Ac-Ac")

ann_text <- data.frame(kcp = 10**c(-2,-1,-2,-4,-4,-4) , krc = 10**c(-6,-6,-6,-4,-4,-4) , lab = rep(c("chi[3] + chi[4] - q[0] > 0"),6), Dim = Dim_t , ForagingMode = ForagingMode_t,phi = rep(c(0.1),6), N = rep(c(0),6))


p <- ggplot(M,aes(x = log10(kcp), y = log10(krc), z =N)) + facet_grid(Dim~ForagingMode)  + stat_contour(aes(colour = phi), breaks = c(0)) + xlab(expression(log[10](k[CP]))) + ylab(expression(log[10](k[RC]))) + labs(colour = expression(phi)) + theme_bw() + coord_cartesian(xlim = c(-13,10)) + geom_text(data = ann_text , aes(label = lab), parse = TRUE, size = 3)

p
    
ggsave("NecCPR.pdf",limitsize=FALSE, height = 8, width = 8)


#Necessity C to P-R

ann_text2 <- data.frame(kcp = 10**c(6, 6, -10 ,6 ,6,-10) , krc = 10**c(-6,-6,8,-4,-4,8) , lab = rep(c("chi[3] > q[0]"),6) , Dim = Dim_t , ForagingMode = ForagingMode_t, phi = rep(c(0.015, 0.15, 2.0 ),2), N = rep(c(0),6))

p <- ggplot(M,aes(x = log10(kcp), y = log10(krc), z = N)) + facet_grid(Dim~phi)  + stat_contour(aes(colour = ForagingMode), breaks = c(0)) + xlab(expression(log[10](k[CP]))) + ylab(expression(log[10](k[RC]))) + labs(colour = expression(Fm)) + theme_bw() + coord_cartesian(xlim = c(-13,10))  + geom_text(data = ann_text2 , aes(label = lab), parse = TRUE, size = 3)

p


