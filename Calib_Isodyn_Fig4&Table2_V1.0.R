#-----------------------------------
#R code to calibrate the Isodyn model
#-----------------------------------

rm(list=ls())# clear the current environment
graphics.off()# clear the current plots

setwd("C:/Users/Sébastien/Desktop/To do/CodeIsodynNEW/GitHub")#your working directory


#Packages to be installed
#install.packages('deSolve')
#install.packages('matrixStats')
#install.packages('lme4')
#install.packages('nls2')

library('lme4') #for Nelder_Mead function
library("deSolve") #for numerical solution of IsoDyn
library("matrixStats") #for ColQuantiles function
library("nls2") #for calibration of the TIM and exponential model

#Compute all needed functions which are put together in a unique file
source("Functions_IsodynV1.0.R")
source("Data_Fig4.R")



#------------------
# Best optimization
#------------------
name<-c('NucheF','Guelinckx','MacAvoy')# names of the treatments

num=2#choose the treatment to use 1=NucheF, 2=Guelinkx, 3=MaCAvoy
data=eval(parse(text=name[num]))
data[["beta"]]=2/3# choice of the allometric coefficient
lower<-c(0.001,0.001,0)#lower boundary for the parameter ri,ro,Ei
upper<-c(0.5,0.5,5)#upper boundary for the parameter ri,ro,Ei
n_best=20#maximum number of loop to reach the global minimum (the best estimation)
crit=0.05#the best evaluation is crit better than the before best evaluation
#Call of the function to perform the calibration


#---------------------------
#Calibration of the Isodyn model
#---------------------------

R<-list(NULL)

model=1 #choice of the model. Model differs on the assumption for Ei and Eo 
#(choice 1 Ei and Eo are equal in absolute values, choice 2 Eo=0 and Ei positive, choice 3 Eo is negative and Ei=0) 

data[["model"]]=model
R[[model]]<-BestNM2(dat=data,lower=lower,upper=upper,n_best=n_best,crit=crit)# calibration for IsoDyn

#R[[1]] to print on screen preliminary results

#-----------------------------
#Calibration of the Time model
#-----------------------------

ResTIM <- nls2(SI~TIM(SI0=mean(SI[tSI==0]),SIDiet=SIdiet,L,Delta,t=tSI), 
              start = expand.grid(L = seq(0.01, 0.05, len = 4),Delta = seq(1, 5, len = 4)), data=data)
y=data$W
x=data$t
ResExpG<-nls2(y~Expgrowth(W0=mean(y[x==0]),k,t=x),start=expand.grid(k = seq(0.05, 0.2, len = 4)))


CVeW<-CV(data=data)[1]# calculation of coefficient of variation for body mass
CVeSI<-CV(data=data)[2]# calculation of coefficient of variation for stable isotopes

#perform bootstrap for parameter se
n_boot=2# number of bootstraps ideal is 500-1000, start with 10 (minimum is 2) but it could last very long with hundreds
ALL<-NULL
ALL<-BOOTNM(dat=data,lower=lower,upper=upper,
          n_best=n_best,crit=crit,n_boot=n_boot,CVeW,CVeSI)# number of bootstraps needed to get a stable sd of parameters

ALL<-ALL[ALL[,1]>0,]#filter unsuccesfull calibration

par_isodyn<-R[[model]]$par#ri,ro,Ei,Eo,Cost,k,Lambda_iso,kg_iso,TEF_iso
par_tim<-c(coef(summary(ResExpG))[1,1],coef(summary(ResTIM))[1,1],coef(summary(ResTIM))[2,1])
se_par_tim<-c(coef(summary(ResExpG))[1,2],coef(summary(ResTIM))[1,2],coef(summary(ResTIM))[2,2])


Res_Tab_par<-par_tab(par_isodyn,par_tim,se_par_tim,matBOOT=ALL)# store the results
Res_GOF_model<-GOF2(par=par_isodyn,par2=par_tim,data=data)# calculate the goodness of fit

Plot_Isodyn(data=data,par=list(ri=R[[model]]$par[1],ro=R[[model]]$par[2],Ei=R[[model]]$par[3],Eo=R[[model]]$par[4]),
            par2=par_tim,BootALL=ALL,name[num])#plot the results, automatically save the plot

#Save results
save(Res_Tab_par, file=paste(paste(name[num],'_Res_tab_par','Model_',data$model,sep=""),"RData",sep="."))
save(Res_GOF_model, file=paste(paste(name[num],'_Res_GOF','Model_',data$model,sep=""),"RData",sep="."))
save(R, file=paste(paste(name[num],'_BestISO_','Model_',data$model,sep=""),"RData",sep="."))
save(ALL, file=paste(paste(name[num],'_BOOT_','Model_',data$model,sep=""),"RData",sep="."))





