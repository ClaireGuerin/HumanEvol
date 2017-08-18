#---- Initialize ----

rm(list = ls())

#Sys.setenv(R_GSCMD = normalizePath("C:/Program Files/gs/gs9.21/bin/gswin64c.exe"))

#library(lattice)
#library(XML)
#library(grImport) # import vector images into plot
#library(Hmisc)
library(jpeg)
library(latex2exp)
library(grDevices)
library(TeachingDemos) # insert a subplot
#library(extrafont)
#loadfonts()

imgPath = "C:/Users/Claire/Dropbox/MEME/Montpellier/Figures/"
setwd(imgPath)

#---- Stratified society ----

viability = function(xVal, muVal, sigmaVal, maxVal = 1) {
  v = maxVal * exp(-(xVal - muVal) ^ 2 / (2 * sigmaVal ^ 2))
  return(v)
}

muV = c(.3,.8)
sigmaV = c(.2,.1)

x = seq(-1,1,length.out = 10000)
vLow = viability(x,muV[1],sigmaV[1],.8)
vHigh = viability(x,muV[2],sigmaV[2])


socialLadder = readJPEG("socialStrata.jpg")


dev.off()
graphics.off() 
pdf('viabilityClass.pdf', width=14, height=7)
par(mfrow = c(1,2))
#layout(mat = matrix(C(1,2)), widths = , heights = )
plot(1:10,type='n', axes=F,xlab='',ylab='',asp=1)
rasterImage(socialLadder,2,2,10,10)

strataCols = rgb(c(0,204),c(204,204),c(203,254), maxColorValue = 255)

plot(c(-.5,1),c(0:1), type='n', xlab=TeX('Diet requirement $X$'), ylab='', cex.lab = 2, cex.axis = 1.5)
mtext(TeX('Viability $V_i(X)$'), 2, cex = 2, adj = .5, line = 2)
lines(vHigh~x, lwd = 3, col = strataCols[1])
lines(vLow~x, lwd = 3, col = strataCols[2] )
legend("topleft",legend = c(TeX('Upper-class'),TeX('Lower-class')),lty = c(1,1),lwd = c(3,3), col = strataCols, bty = 'n', seg.len = .5, x.intersp = .5, cex = 1.5)
dev.off()

#---- Male mate Choice ----

x = seq(0,1,length.out = 10000)
mateAccess = (1 - x ^ 2)
reproSuccess = mateAccess*x


dev.off()
graphics.off() 
pdf('choosinessFunction.pdf', width=7, height=7)

plot(reproSuccess~x, type = 'l', lwd = 2, col = 'grey', xlab = TeX('Male choosiness level $\\phi$'), ylab = '', asp=1, xlim = c(0,1), ylim = c(0,.6), cex.lab = 1.5)
mtext(TeX('Male reproductive success $W(\\varphi)$'),2, cex = 1.5, adj = .5, line = 2)

subplot(plot(x,mateAccess, type = 'l', col = 'darkred', asp=1, cex.axis=0.5),'topright', size = c(1.8,1.8), type='fig', pars=list(oma=c(0,0,0,0), mar=c(0,0,0,0), tcl=-0.1, mgp=c(0,0,0)))


text(0.65,0.55,labels = TeX('Female pool portion $P(\\varphi)$'), srt=90, cex =.9)
dev.off()
