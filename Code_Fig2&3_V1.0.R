#-----------------------
#Code to produce Fig2&3
#----------------------


#install.packages("deSolve")
library("deSolve")


setwd("C:/Users/Sébastien/Desktop/To do/CodeIsodynNEW/GitHub")#Working directory


rm(list=ls())# clear the current environment
graphics.off()# clear the current plots

#####"
#Function
######

source("Functions_IsoDynV1.0.R")# File hosting all needed functions to run the code


############
##General properties of IsoDyn model
##Fig 2
###########

time<-seq(0,600,0.1)# in days
Ei<-2#Isotopic discrimination on the flux of body mass gains 
Eo<--2#Isotopic discrimination on the flux of body mass losses

W0=0.1# initial body mass of all virtual species
r_ref=0.2# a reference rate
ri=c(r_ref,r_ref/2,r_ref*3)#Rate of body mass gains (or inputs)
ro=c(r_ref/4,r_ref/8,r_ref)#Rate of body mass losses (or outputs) 

#(ri/ro)^3#Max body mass of the three species respectively

namesp<-c("Sp1","Sp2","Sp3")# names of the species

#------------------
#body mass dynamics
#Fig2 A
#------------------

marge<-c(4, 5, 3, 1)#c(bottom, left, top, right)
s<-1.5
par(mar = marge, col='black',col.axis = 'black',col.lab='black')
par(cex.axis=s,cex.lab=s,las = 1,tck=-0.02)
Color<-c("slategray2","slategray3","slategray4")

Xname1<-c(expression(paste("Time"," (d)")))
Yname1<-c(expression(paste("Body mass"," (g)")))
plot(time,IGMW(W0=W0,ri=ri[1],ro=ro[1],beta=2/3,time),type="l",xlab=Xname1,ylab=Yname1,lwd=5,cex=s,col=Color[1])
lines(time,IGMW(W0=W0,ri=ri[2],ro=ro[2],beta=2/3,time),lwd=5,cex=s,col=Color[2])
lines(time,IGMW(W0=W0,ri=ri[3],ro=ro[3],beta=2/3,time),lwd=5,cex=s,col=Color[3])
axis(side=2, lwd = 2, col.ticks = 'black')
axis(side = 1, lwd = 2,col.ticks = 'black')
box(lwd = 2)
mtext("A.", side=3, line=1, adj=0.0, cex=1.4, font=2)

dev.copy(png,'Fig2fakeanimalW.png')
dev.off()

#------------------
#Lambda dynamics
#Fig2 B
#------------------


l1<-c(namesp[1],namesp[2],namesp[3])
Xname3<-c(expression(paste("Ln (Body mass)")))
Yname3<-c(expression(paste("Ln(", lambda,")")))

plot(NULL,xlab=Xname3,ylab=Yname3,lwd=5,cex=s,ylim=c(-4,1),xlim=c(-3,5))

for (i in 1:3){
    lines(log(IGMW(W0=W0,ri=ri[i],ro=ro[i],beta=2/3,time)),
        log(ri[i]*(IGMW(W0=W0,ri=ri[i],ro=ro[i],beta=2/3,time))^(-1/3)),
        lwd=5,cex=s,col=Color[i])
}
axis(side=2, lwd = 2, col.ticks = 'black')
axis(side = 1, lwd = 2,col.ticks = 'black')
box(lwd = 2)
legend("topright", inset = .02, legend = l1, col=c("slategray2","slategray3","slategray4"),
       lty = 1, lwd=2, bg=NA, cex = 1.2, bty="n")
mtext("B.", side=3, line=1, adj=0.0, cex=1.4, font=2)

dev.copy(png,'Fig2fakeanimalLambda.png')
dev.off()

#------------------
#TDF vs kg
#Fig2 B
#------------------


Xname3<-c(expression(paste("Specific growth rate",", k"[g]," (1/d)")))
Yname3<-c(expression(paste(delta^15,"N"[b],"-",delta^15,"N"[d])))

plot(NULL,xlab=Xname3,ylab=Yname3,lwd=5,cex=s,ylim=c(0,-Eo+Ei+1),xlim=c(0,1))

for (i in 1:3){
  kg<-ri[i]*(IGMW(W0=W0,ri=ri[i],ro=ro[i],beta=2/3,time))^(-1/3)-ro[i]
  lines(kg,runisodyn(time=time,state0=c(SI=0),
                     param=c(W0=W0,SIdiet=0,ri=ri[i],ro=ro[i],beta=2/3,deltai=Ei,deltao=Eo)),
        lwd=5,cex=s,col=Color[i])
}
axis(side=2, lwd = 2, col.ticks = 'black')
axis(side = 1, lwd = 2,col.ticks = 'black')
box(lwd = 2)
mtext("C.", side=3, line=1, adj=0.0, cex=1.4, font=2)

dev.copy(png,'Fig2fakeanimalTEFvskg.png')
dev.off()


#######################
#Diet Switch experiment 
#Fig 3
#######################

time2<-seq(0,100,0.1)

Color<-c("slategray3","slategray4")
leg<-c("juvenile","adult")
namepng<-c('Fig3DSE_sp1.png','Fig3DSE_sp3.png','Fig3DSE_sp3.png')

s<-1.5
Xname3<-c(expression(paste("Time (d)")))
Yname3<-c(expression(paste(delta^15,"N"[b])))


for (j in 1:3){
Wdse<-c(W0,(ri[j]/ro[j])^(3))

par(mar = marge, col='black',col.axis = 'black',col.lab='black')
par(cex.axis=s,cex.lab=s,las = 1,tck=-0.02)#c(bottom, left, top, right)

plot(NULL,xlab=Xname3,ylab=Yname3,lwd=5,cex=s,ylim=c(0,-Eo+Ei+3),xlim=c(0,time2[length(time2)]))
abline(a=2,b=0, lty="longdash",lwd=5,col="slategray2")

for (i in 1:2){

lines(time2,runisodyn(time=time2,state0=c(SI=0),
                     param=c(W0=Wdse[i],SIdiet=2,ri=ri[j],ro=ro[j],beta=2/3,deltai=Ei,deltao=Eo)),
      lwd=5,cex=s,col=Color[i])
}
box(lwd = 2)
legend("bottomright", inset = .01, legend = leg , col=Color,
       lty = 1, lwd=2, bg=NA, cex = s, bty="n")
mtext(namesp[j], side=3, line=1, adj=0.0, cex=1.4, font=2)
dev.copy(png,namepng[j])
dev.off()
}






