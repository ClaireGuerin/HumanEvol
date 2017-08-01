rm(list = ls())

library(car)
library(extrafont)
loadfonts()

x = seq(0,1,length.out = 10000)
y = (1 - x ^ 2)*x

imgPath = "C:/Users/Claire/Dropbox/MEME/Montpellier/AdaptiveDynamicsStratification/choosinessFunction.pdf"
#pdf(imgPath, family="CM Roman Greek", width=7, height=7)
pdf(imgPath, width=7, height=7)

plot(x,y, type = 'l', lwd = 2, col = 'grey', xlab = expression(paste('Male choosiness level ', varphi)), ylab = 'Male reproductive success')

dev.off()