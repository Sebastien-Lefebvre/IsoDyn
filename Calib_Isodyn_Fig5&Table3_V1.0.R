############################
#Isodyn with Von Bertalanffy
############################


#install.packages("deSolve")
library("deSolve")


setwd("C:/Users/Sébastien/Desktop/To do/CodeIsodynNEW/GitHub")#your working directory

rm(list=ls())# clear the current environment
graphics.off()# clear the current plots


#--------------------------------------------------------
# functions dedicated to link IsoDyn model with DEB Theory
#--------------------------------------------------------

source(file = 'Functions_IsodynV1.0.R')


r_o<-function(f,g,k_M,Tcor){3*1/3/(1+f/g)*k_M*Tcor}
W_m<-function(f,Lm,sm,w){(f*Lm*sm)^3*(1+f*w)}
Temp_cor<-function(Te,T_A,T_ref){exp(T_A/T_ref-T_A/(Te+273.15))}

IGM_DEB<-function(W0,f,g,k_M,Tcor,Lm,sm,w,t){
  ro<-r_o(f,g,k_M,Tcor)
  Wm<-W_m(f,Lm,sm,w)
  W<-(Wm^(1/3)-(Wm^(1/3)-W0^(1/3))*exp(-ro/3*t))^3
return(W)}

#---------------------------------------------------------
#Function to be minimized by the Optim function
#--------------------------------------------------------

min.func_f_DEB <- function(data, par) {
  with(data, 
       sbLOSS(W,IGM_DEB(W0,f=par[1],g,k_M,Tcor,Lm,sm,w,t))#loss func for weight
  )}


min.func_SI4 <- function(data, par) {#loss func for iso
  with(data, 
       sbLOSS(c(0,TEF[1]),runisodyn(time=c(0,t[1]),state0=c(SI=fb),
              param=c(W0=Wi[1],SIdiet=fd,ri=ri[1],ro=ro[1], beta=2/3,deltai=par[1],deltao=par[2])))
       +sbLOSS(c(0,TEF[2]),runisodyn(time=c(0,t[2]),state0=c(SI=fb),
              param=c(W0=Wi[2],SIdiet=fd,ri=ri[2],ro=ro[2], beta=2/3,deltai=par[1],deltao=par[2])))
       +sbLOSS(c(0,TEF[3]),runisodyn(time=c(0,t[3]),state0=c(SI=fb),
              param=c(W0=Wi[3],SIdiet=fd,ri=ri[3],ro=ro[4], beta=2/3,deltai=par[1],deltao=par[2])))
       +sbLOSS(c(0,TEF[4]),runisodyn(time=c(0,t[4]),state0=c(SI=fb),
              param=c(W0=Wi[4],SIdiet=fd,ri=ri[4],ro=ro[4], beta=2/3,deltai=par[1],deltao=par[2])))
        )}


#----------------------------------------------------------------------
#DSE simulation and parameter calibrations of Carp case study using DEB
#----------------------------------------------------------------------

#Data retrived from Gaye-Siesseger et al.(2004) using Plot reader

CarpObs<-list(t=c(56,56,56,56),
               Wi=c(30.5,31.4,28.5,28.3),
               Wf=c(31.8,65.1,109.4,152.1),
               kg=c(0.00074535,0.01301994,0.024019764,0.030029936),
               TEF=c(1.7,1.4,1.2,1.0),
               fb=0,fd=0)

#----------------
#DEB parameters
#--------------

# Retrieved from Add my pet collection using ESM 2

Lm<-6.3743
sm<-1.2236
k_M<-0.0698
T_A<-8000
T_ref<-293.15
w<-99.37
g<-0.0337
beta<-2/3

T_exp<-27
Tcor<-Temp_cor(T_exp,T_A,T_ref)
ts<-seq(0,56,1)



#Optimization of the functional scaled response (f) on growth data

f<-NULL
for (i in 1:4){
resf<-optim(par = c(0.16), fn = min.func_f_DEB, lower=0, upper=1.5, method='Brent',
            data = list(W=CarpObs$Wf[i],W0=CarpObs$Wi[i],g=g,k_M=k_M,Tcor=Tcor,Lm=Lm,sm=sm,w=w,t=CarpObs$t[i]),control=list(maxit=2000))
f<-cbind(f,resf$par)}


#Plot of body mass dynamics
marge<-c(4, 5, 3, 1)#c(bottom, left, top, right)
s<-1.5
par(mar = marge, col='black',col.axis = 'black',col.lab='black')
par(cex.axis=s,cex.lab=s,las = 1,tck=-0.02)
Color<-c("slategray1","slategray2","slategray3","slategray4")
leg<-c("L1","L2","L3","L4")
Xname1<-c(expression(paste("Time (d)")))
Yname1<-c(expression(paste("Body mass (g)")))

plot(NULL,xlim=c(0,60),ylim=c(0,200),type="l",xlab=Xname1,ylab=Yname1)
for (i in 1:4){lines(ts,IGM_DEB(W0=CarpObs$Wi[i],f[i],g,k_M,Tcor,Lm,sm,w,ts),lwd=5,cex=s,col=Color[i])
               points(c(0,56),c(CarpObs$Wi[i],CarpObs$Wf[i]),pch=19,cex=2,bg=1,col=Color[i])}

axis(side=2, lwd = 2, col.ticks = 'black')
axis(side = 1, lwd = 2,col.ticks = 'black')
box(lwd = 2)
legend("topleft", inset = .02, legend = leg, col=Color, lty = 1, lwd=2, bg=NA, cex = 1, bty="n")
mtext("A.", side=3, line=1, adj=0.0, cex=1.4, font=2)

dev.copy(png,'Fig5CarpW.png')
dev.off()


# calculation of IsoDyn parameters (ri and ro) from DEB parameters using the optimzed f values 
ro<-c(r_o(f=f[1],g,k_M,Tcor),r_o(f=f[2],g,k_M,Tcor),r_o(f=f[3],g,k_M,Tcor),r_o(f=f[4],g,k_M,Tcor))
Wm<-c(W_m(f=f[1],Lm,sm,w),W_m(f=f[2],Lm,sm,w),W_m(f=f[3],Lm,sm,w),W_m(f=f[4],Lm,sm,w))
ri<-ro*Wm^(1/3)

resEo2<-optim(par = c(0.8,-1), fn = min.func_SI4, data=CarpObs,control=list(maxit=2000))

Ei<-resEo2$par[1]#1
Eo<-resEo2$par[2]#-1.5

marge<-c(4, 5, 3, 1)#c(bottom, left, top, right)
s<-1.5
par(mar = marge, col='black',col.axis = 'black',col.lab='black')
par(cex.axis=s,cex.lab=s,las = 1,tck=-0.02)
Color<-c("slategray1","slategray2","slategray3","slategray4")
#l1<-c(expression(paste("ri=ro=0.2")),expression(paste("ri=ro=0.1")),expression(paste("ri=0.16,ro=0.24")))
Xname3<-c(expression(paste("Time (d)")))
Yname3<-c(expression(paste(delta^15,"N"[b],"-",delta^15,"N"[d])))

plot(NULL,xlim=c(0,60),ylim=c(0,2),type="l",xlab=Xname3,ylab=Yname3)


for (i in 1:4){
lines(ts,runisodyn(time=ts,state0=c(SI=CarpObs$fb),
        param=c(W0=CarpObs$Wi[i],SIdiet=CarpObs$fd,ri=ri[i],ro=ro[i],beta=2/3,deltai=Ei,deltao=Eo)),lwd=5,cex=s,col=Color[i])
points(c(0,56),c(0,CarpObs$TEF[i]),pch=19,cex=2,bg=1,col=Color[i])}


axis(side=2, lwd = 2, col.ticks = 'black')
axis(side=1, lwd = 2,col.ticks = 'black')
box(lwd = 2)
mtext("B.", side=3, line=1, adj=0.0, cex=1.4, font=2)

dev.copy(png,'Fig5CarpISO.png')
dev.off()


X1<-runisodyn(time=ts,state0=c(SI=CarpObs$fb),
          param=c(W0=CarpObs$Wi[1],SIdiet=CarpObs$fd,ri=ri[1],ro=ro[1],beta=2/3,deltai=Ei,deltao=Eo))
X2<-runisodyn(time=ts,state0=c(SI=CarpObs$fb),
              param=c(W0=CarpObs$Wi[2],SIdiet=CarpObs$fd,ri=ri[2],ro=ro[2],beta=2/3,deltai=Ei,deltao=Eo))
X3<-runisodyn(time=ts,state0=c(SI=CarpObs$fb),
              param=c(W0=CarpObs$Wi[3],SIdiet=CarpObs$fd,ri=ri[3],ro=ro[3],beta=2/3,deltai=Ei,deltao=Eo))
X4<-runisodyn(time=ts,state0=c(SI=CarpObs$fb),
              param=c(W0=CarpObs$Wi[4],SIdiet=CarpObs$fd,ri=ri[4],ro=ro[4],beta=2/3,deltai=Ei,deltao=Eo))

s<-1.5
par(mar = marge, col='black',col.axis = 'black',col.lab='black')
par(cex.axis=s,cex.lab=s,las = 1,tck=-0.02)#c(bottom, left, top, right)
Xname3<-c(expression(paste("Specific growth rate, k"[g]," (1/d)")))
Yname3<-c(expression(paste(delta^15,"N"[b],"-",delta^15,"N"[d])))
plot(CarpObs$kg,c(X1[56],X2[56],X3[56],X4[56]),type="l",xlab=Xname3,ylab=Yname3,lwd=5,cex=s,col="slategray2",ylim=c(0.5,2.5))
points(CarpObs$kg,CarpObs$TEF,lwd=5,col="slategray3",pch=19,cex=s,bg=1)

axis(side=2, lwd = 2, col.ticks = 'black')
axis(side = 1, lwd = 2,col.ticks = 'black')
box(lwd = 2)
mtext("C.", side=3, line=1, adj=0.0, cex=1.4, font=2)

dev.copy(png,'Fig5CarpTDF.png')
dev.off()



