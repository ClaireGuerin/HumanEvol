rm(list = ls())

library(TeachingDemos)
library(car)
library(extrafont)
loadfonts()

x = seq(0,1,length.out = 10000)
mateAccess = (1 - x ^ 2)
reproSuccess = mateAccess*x

imgPath = "C:/Users/Claire/Dropbox/MEME/Montpellier/AdaptiveDynamicsStratification/choosinessFunction.pdf"
#pdf(imgPath, family="CM Roman Greek", width=7, height=7)
pdf(imgPath, width=7, height=7)

X11()
plot(x,reproSuccess, type = 'l', lwd = 2, col = 'grey', xlab = expression(paste('Male choosiness level ', varphi)), ylab = 'Male reproductive success',asp=1)
par(mar = rep(2, 4))
subplot(plot(x,mateAccess, type = 'l', col = 'grey', asp=1, cex.axis=0.5),"topright",type='fig')


dev.off()