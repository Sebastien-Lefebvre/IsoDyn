#########################################################
#Collection of funtions to run IsoDyn and its calibration
#########################################################


#---------------------------------------------------------
#Function to calculte mean or sd by column
#--------------------------------------------------------
Meanby<-function(obs){Meanby<-data.matrix(cbind(as.numeric(names(tapply(obs[,2],obs[,1],mean))),as.numeric(tapply(obs[,2],obs[,1],mean))))}
Sdby<-function(obs){Meanby<-data.matrix(cbind(as.numeric(names(tapply(obs[,2],obs[,1],sd))),as.numeric(tapply(obs[,2],obs[,1],sd))))}


#---------------------------------------------------------
#Function of the exponential growth model
#--------------------------------------------------------
Expgrowth<-function(W0,k,t){
  W<-W0*exp(k*t)
  return(W)}

#---------------------------------------------------------
#Function of the time incorporation model
#--------------------------------------------------------
TIM<-function(SI0,SIDiet,L,Delta,t){
  SI_inf<-SIDiet+Delta
  SI<-SI_inf-(SI_inf-SI0)*exp(-L*t)
  return(SI)}


#---------------------------------------------------------
#Function of the Individual growth model in body mass
#--------------------------------------------------------
IGMW<-function(W0,ri,ro,beta,t){#
  W<-(ri/ro-(ri/ro-W0^(1-beta))*exp(-ro*(1-beta)*t))^(1/(1-beta))
  return(W)}  

#---------------------------------------------------------
#Function of the IsoDyn model
#--------------------------------------------------------
isodyn <- function(t,state,parameters){# the model equations of Isodyn:ODEs 
  with(as.list(c(state,parameters)),{  # unpack the state variables, parameters
    #state variables ODEs
    W <- IGMW(W0,ri,ro,beta,t) # body mass at time t using IGMW function
    dSI<- ri*W^(beta-1)*(SIdiet-SI+deltai)-deltao*ro #rate of change of stable isotopes of the consumer
    # the output, packed as a list
    list(c(dSI))# the rate of change
  })}  # end of model

#---------------------------------------------------------
#Function to run IsoDyn model
#--------------------------------------------------------

runisodyn <-function(time,state0,param){#  time, initial values of state variables and parameters 
    out<-as.data.frame(lsoda(state0,time,isodyn,param)) # ode or lsoda are integrator
  return(out$SI)}

#---------------------------------------------------------
#Function to run the loss function
#--------------------------------------------------------
sbLOSS<-function(obs,pred){sum(1/length(obs)*(pred-obs)^2/(mean(pred)^2+mean(obs)^2))}
# dimensionless loss function accounting for the difference between predictions and observations

#---------------------------------------------------------
#Function to calculate the relative error
#--------------------------------------------------------
RelE<-function(obs,pred){sum(1/length(obs)*abs(pred-obs)/mean(obs))}

#---------------------------------------------------------
#Function to be minimized by the Nelder_Mead function
#--------------------------------------------------------
min.funcNM <- function(par) {
  # when Ei=-EO
  deltai=par[3]
  deltao=-par[3]
  if (data$model==2){deltao=0}
  if (data$model==3){deltai=0}
  
  #loss func for body mass
  Cost_W<-sbLOSS(d$W,IGMW(W0=mean(d$W[d$t==0]),ri=par[1],ro=par[2],beta=d$beta,t=d$t))
  #loss func for SI with Ei=-Eo then only 3 parameters to calibrate (ri,ro and Ei)
  Cost_SI<-sbLOSS(d$SI,runisodyn(d$tSI,state0=c(SI=mean(d$SI[d$tSI==0])),param=c(W0=mean(d$W[d$t==0])
                    ,SIdiet=d$SIdiet,ri=par[1],ro=par[2],beta=d$beta,deltai=deltai,deltao=deltao)))
  return(Cost_W+Cost_SI)}


#---------------------------------------------------------
#Function to calibrate IsoDyn model
#--------------------------------------------------------


BestNM2<-function(dat,lower,upper,n_best,crit){#function to calibrate ri, ro and Ei=-Eo
  
  ALL<-NULL#Matrix to store the results
  k=0
  ri<-0;ro<-0;Ei<-0;Eo<-0;Cost<-0
  
  ri_start<-seq(lower[1],upper[1],0.05)
  ro_start<-seq(lower[2],upper[2],0.05)
  Ei_start<-seq(lower[3],upper[3],0.1)

  init<-expand.grid(ri=ri_start,ro=ro_start,Ei=Ei_start)
  
  d<<-dat#make d a global variable accessible outside the function
  repeat{#loop to get n estimates with randomized initial values
    k=k+1
    rd1<-floor(runif(1, min=1, max=length(init[,1])))#random selection in init

    result<-Nelder_Mead(min.funcNM, par=as.numeric(init[rd1,]), lower = lower, upper = upper)
    result<-Nelder_Mead(min.funcNM, par=c(result$par[1],result$par[2],result$par[3]), lower = lower, upper = upper)

    if (result$convergence==0) {ALL<-rbind(ALL,c(result$par,result$fval))}#implement the result in the matrix ALL
    print(paste("n_best=",k))
    
    #Filtering the results by deleting parameter values equal to the boundaries
    for (j in 1:3){ALL<-subset(ALL,ALL[,j]>lower[j]);ALL<-subset(ALL,ALL[,j]<upper[j])}
    
    if (length(ALL[,4])>1){ALL<-ALL[order(ALL[,4]),]} #sorting results  
    
    if (length(ALL[,4])>1) {if ((ALL[2,4]-ALL[1,4])/ALL[1,4] < crit){break}} #stop calibration if optimisation is satisfactory
    
    if(k==n_best){ break }
  } # end of repeat
  
  
  if (length(ALL)>0){#assign the best result if fitting occurred
    ri<-ALL[1,1]
    ro<-ALL[1,2]
    Ei<-ALL[1,3]
    Eo<--ALL[1,3]
    if (d$model==2) {Eo=0}
    if (d$model==3) {Ei=0}
    Cost<-ALL[1,4]
  }
  
  # Calculate estra parameters
  Lambda_iso<-sum(ri*IGMW(W0=mean(d$W[data$t==0]),
                               ri=ri,ro=ro,beta=d$beta,t=d$t)^(d$beta-1))/length(d$t)
  kg_iso<-Lambda_iso-ro
  TEF_iso<-Ei-Eo
  
  return(list(par=c(ri,ro,Ei,Eo,Cost,k,Lambda_iso,kg_iso,TEF_iso),ALL=ALL))}



#---------------------------------------------------------
#Function to perform bootstrap
#--------------------------------------------------------


BOOTNM<-function(data,lower,upper,n_best,crit,n_boot,CVeW,CVeSI){ 
  
  ALLBOOT<-NULL
  
  
  for (j in 1:n_boot){
    print(paste("n_boot=",j))
    res<-NULL
    databoot<-data#creat a new list of observations
    #randomize residual and generate new W and Si obs with Monte Carlo simulation based on centered log-normally distributed residuals
    databoot$W<-databoot$W*(1+(rlnorm(length(databoot$W),meanlog=0,sdlog=CVeW)-1))
    #randomize residual and generate new W obs
    databoot$SI<-databoot$SI*(1+(rlnorm(length(databoot$SI),meanlog=0,sdlog=CVeSI)-1))#randomize residual and generate new N obs
    
    res<-BestNM2(dat=databoot,lower=lower,upper=upper,n_best=n_best,crit=crit)

    ri<-res$par[1]
    ro<-res$par[2]
    Ei<-res$par[3]
    Eo<-res$par[4]
    Cost<-res$par[5]
    n_loop<-res$par[6]
    Lambda_iso<-res$par[7]
    kg_iso<-res$par[8]
    TDF_iso<-res$par[9]
    
    #Store the results
    ALLBOOT<-rbind(ALLBOOT,c(ri,ro,Ei,Eo,Cost,n_loop,Lambda_iso,kg_iso,TDF_iso))
    
  }#end of bootstrap loop
  
  return(ALLBOOT)}


#---------------------------------------------------------
#Function to calculate coefficient of variation
#--------------------------------------------------------


CV<-function(data){
  CVeW<-data$CVeW
  CVeSI<-data$CVeSI
  
  if (data$CVeW==0){CVeW<-0.15}# expert knowledge max 0.2
  if (data$CVeSI==0){CVeSI<-0.05}# expert knowledge max 0.1
  
  if (data$CVeW==1){
  CVeW<-mean(na.omit(Sdby(cbind(data$t,data$W))[,2]/Meanby(cbind(data$t,data$W))[,2]))}
  
  if (data$CVeSI==1){
  CVeSI<-mean(na.omit(Sdby(cbind(data$tSI,data$SI))[,2]/Meanby(cbind(data$tSI,data$SI))[,2]))}
  
  return(c(CVeW,CVeSI))}

#---------------------------------------------------------
#Function to assemble all parameter values
#--------------------------------------------------------

par_tab<-function(par,par2,se_par2,matBOOT){

se_kc<-sqrt((se_par2[1]/par2[1])^2+(se_par2[2]/par2[2])^2)*(par2[2]-par2[1])

CIlow<-colQuantiles(as.matrix(matBOOT), probs = 0.025)
CIhigh<-colQuantiles(as.matrix(matBOOT), probs = 0.975)

CI_low_kgtim<-par2[1]+se_par2[1]*qt(p=0.025,df=length(data$W)-1,lower.tail=T)#df is n obs minus n parameter of the model
CI_low_Ltim<-par2[2]+se_par2[2]*qt(p=0.025,df=length(data$SI)-2,lower.tail=T)
CI_low_TDFtim<-par2[3]+se_par2[3]*qt(p=0.025,df=length(data$SI)-2,lower.tail=T)
CI_low_kctim<-(par2[2]-par2[1])+se_kc*qt(p=0.025,df=length(data$SI)+length(data$W)-3,lower.tail=T)

CI_high_kgtim<-par2[1]+se_par2[1]*qt(p=0.975,df=length(data$W)-1,lower.tail=T)
CI_high_Ltim<-par2[2]+se_par2[2]*qt(p=0.975,df=length(data$SI)-2,lower.tail=T)
CI_high_TDFtim<-par2[3]+se_par2[3]*qt(p=0.975,df=length(data$SI)-2,lower.tail=T)
CI_high_kctim<-(par2[2]-par2[1])+se_kc*qt(p=0.975,df=length(data$SI)+length(data$W)-3,lower.tail=T)

res<-data.frame(Parameter=c('ri','ro','Ei','Eo','Lambda_ISO','kg_ISO','TDF_ISO','kg_TIM','Lambda_TIM','TDF_TIM','kc_TIM'),
                Estimates=c(par[1],par[2],par[3],par[4],par[7],par[8],par[9],par2[1],par2[2],par2[3],par2[2]-par2[1]),
                CIlow=c(CIlow[1],CIlow[2],CIlow[3],CIlow[4],CIlow[7],CIlow[8],CIlow[9],CI_low_kgtim,CI_low_Ltim,CI_low_TDFtim,CI_low_kctim),
                CIhigh=c(CIhigh[1],CIhigh[2],CIhigh[3],CIhigh[4],CIhigh[7],CIhigh[8],CIhigh[9],CI_high_kgtim,CI_high_Ltim,CI_high_TDFtim,CI_high_kctim))
return(res)}


#---------------------------------------------------------
#Function to calculte the Goodness of fit
#--------------------------------------------------------


GOF2<-function(par,par2,data){
  
  Pred_W<-IGMW(W0=mean(data$W[data$t==0]),ri=par[1],ro=par[2],beta=data$beta,t=data$t)
  Pred_SI<-runisodyn(data$tSI,state0=c(SI=mean(data$SI[data$tSI==0])),param=c(W0=mean(data$W[data$t==0])
                                                                              ,SIdiet=data$SIdiet,ri=par[1],ro=par[2],beta=data$beta,deltai=par[3],deltao=par[4]))
  Pred_W_class<-Expgrowth(W0=mean(data$W[data$t==0]),k=par2[1],t=data$t)
  Pred_SI_class<-TIM(SI0=mean(data$SI[data$tSI==0]),SIDiet=data$SIdiet,L=par2[2],Delta=par2[3],t=data$tSI)
  
  SSE_W<-sbLOSS(data$W,Pred_W)
  RE_W<-RelE(data$W,Pred_W)
  SSE_SI<-sbLOSS(data$SI,Pred_SI)
  RE_SI<-RelE(data$SI,Pred_SI)
  
  SSE_W_class<-sbLOSS(data$W,Pred_W_class)
  RE_W_class<-RelE(data$W,Pred_W_class)
  SSE_SI_class<-sbLOSS(data$SI,Pred_SI_class)
  RE_SI_class<-RelE(data$SI,Pred_SI_class)
  
  res<-data.frame(Model=c('Body_mass','Stable isotopes','Mean'),RE_IsoDyn=c(RE_W,RE_SI,mean(c(RE_W,RE_SI)))
                  ,SSE_IsoDyn=c(SSE_W,SSE_SI,mean(c(SSE_W,SSE_SI))),RE_class=c(RE_W_class,RE_SI_class,mean(c(RE_W_class,RE_SI_class)))
                  ,SSE_class=c(SSE_W_class,SSE_SI_class,mean(c(SSE_W_class,SSE_SI_class))))
  return(res)
}

#---------------------------------------------------------
#Function to plot the body mass and stable isotope results
#--------------------------------------------------------
Plot_Isodyn<-function(data,par,par2,BootALL,name){
  Predboot_SI<-NULL
  Predboot_W<-NULL
  CIpred_SI<-NULL
  CIpred_W <-NULL
  
  tsimuW<-seq(0,data$t[length(data$t)],1)
  tsimuSI<-seq(0,data$tSI[length(data$tSI)],1)
  
  if (length(BootALL[,1])>0){
  for (i in 1:length(BootALL[,1])){
    Predboot_SI<-rbind(Predboot_SI,runisodyn(time=tsimuSI,state0=c(SI=mean(data$SI[data$tSI==0])),
                                             param=c(W0=mean(data$W[data$t==0]),SIdiet=data$SIdiet,ri=BootALL[i,1],ro=BootALL[i,2],beta=data$beta,deltai=BootALL[i,3],deltao=-BootALL[i,3])))
    Predboot_W<-rbind(Predboot_W,IGMW(W0=mean(data$W[data$t==0]),ri=BootALL[i,1],ro=BootALL[i,2],beta=data$beta,t=tsimuW))
  }
  
  CIpred_SI<-rbind(colQuantiles(as.matrix(Predboot_SI), probs = 0.95),colQuantiles(as.matrix(Predboot_SI), probs = 0.05))
  CIpred_W<-rbind(colQuantiles(as.matrix(Predboot_W), probs = 0.95),colQuantiles(as.matrix(Predboot_W), probs = 0.05)) 
  }
  
  Color<-c(rgb(0.78,0.89,1,0.2),rgb(0.73,0.83,0.93,0.2),rgb(0.62,0.71,0.8,0.2),rgb(0.42,0.48,0.55,0.2))#
  Color2<-c(rgb(0.78,0.89,1,1),rgb(0.73,0.83,0.93,1),rgb(0.62,0.71,0.8,1),rgb(0.42,0.48,0.55,1))#
  marge<-c(4, 5.5, 3, 1)#c(bottom, left, top, right)
  s<-1.5
  Yname_W<-c(expression(paste("Body mass (g)")))
  Yname_SI<-c(expression(paste(delta^15,"N"[m]," (","\211",")")))#"N"[m]
  par(mar = marge, col='black',col.axis = 'black',col.lab='black')
  par(cex.axis=s,cex.lab=s,las = 1,tck=-0.02)
  XLIM_W = c(0,data$t[length(data$t)]*1.1)
  YLIM_W = c(mean(data$W[data$t==0])*0.9,data$W[length(data$W)]*1.1)
  XLIM_SI = c(0,data$tSI[length(data$tSI)]*1.1)
  MinSI = min(data$SIdiet-1,mean(data$SI[data$tSI==0]),mean(data$SI[data$tSI==data$tSI[length(data$tSI)]]))
  MaxSI = max(mean(data$SI[data$tSI==0]),mean(data$SI[data$tSI==data$tSI[length(data$tSI)]]))
  YLIM_SI = c(MinSI*0.9,MaxSI*1.1)
  
  
  #Body mass
  plot(NULL,xlim=XLIM_W,ylim=YLIM_W,xlab='Time (d)',ylab=Yname_W,cex=s)
  if(length(CIpred_W)>0){polygon(c(tsimuW,rev(tsimuW)),c(CIpred_W[1,],rev(CIpred_W[2,])),col=Color[3],border=NA)}
  points(Meanby(cbind(data$t,data$W))[,1],Meanby(cbind(data$t,data$W))[,2],pch=1,col=Color2[4],cex=s)
  segments(Meanby(cbind(data$t,data$W))[,1],Meanby(cbind(data$t,data$W))[,2]+Sdby(cbind(data$t,data$W))[,2],
           Meanby(cbind(data$t,data$W))[,1],Meanby(cbind(data$t,data$W))[,2]-Sdby(cbind(data$t,data$W))[,2],
           col=Color2[4],cex=s)
  lines(tsimuW,IGMW(W0=mean(data$W[data$t==0]),ri=par$ri,ro=par$ro,beta=data$beta,t=tsimuW),col=Color2[4],lwd=3,cex=s)
  lines(data$t,Expgrowth(W0=mean(data$W[data$t==0]),k=par2[1],t=data$t),col=Color2[4],lwd=3,lty=2,cex=s)
  axis(side=2, lwd = 2, col.ticks = 'black')
  axis(side=1, lwd = 2, col.ticks = 'black')
  box(lwd = 2)

  dev.copy(png,paste(paste(name,'_Body_mass','_Model_',data$model,sep=""),"png",sep="."))
  dev.off()
  
  #SI
  plot(NULL,xlim=XLIM_SI,ylim=YLIM_SI,xlab='Time (d)',ylab=Yname_SI,cex=s)
  abline(a=data$SIdiet,b=0, lty="longdash",lwd=2,col="grey")#color[4]
  if(length(CIpred_SI)>0){polygon(c(tsimuSI,rev(tsimuSI)),c(CIpred_SI[1,],rev(CIpred_SI[2,])),col=Color[4],border=NA)}
  points(Meanby(cbind(data$tSI,data$SI))[,1],Meanby(cbind(data$tSI,data$SI))[,2],pch=1,col=Color2[4],cex=s)
  segments(Meanby(cbind(data$tSI,data$SI))[,1],Meanby(cbind(data$tSI,data$SI))[,2]+Sdby(cbind(data$tSI,data$SI))[,2],
           Meanby(cbind(data$tSI,data$SI))[,1],Meanby(cbind(data$tSI,data$SI))[,2]-Sdby(cbind(data$tSI,data$SI))[,2],
           col=Color2[4],cex=s)
  lines(tsimuSI,runisodyn(time=tsimuSI,state0=c(SI=mean(data$SI[data$tSI==0])),
                        param=c(W0=mean(data$W[data$t==0]),SIdiet=data$SIdiet,ri=par$ri,ro=par$ro,beta=data$beta,deltai=par$Ei,deltao=par$Eo))
        ,col=Color2[4],lwd=3,cex=s)
  lines(tsimuSI,TIM(SI0=mean(data$SI[data$tSI==0]),SIDiet=data$SIdiet,L=par2[2],Delta=par2[3],t=tsimuSI)
        ,col=Color2[4],lwd=3,lty=2, cex=s)
  axis(side=2, lwd = 2, col.ticks = 'black')
  axis(side = 1, lwd = 2,col.ticks = 'black')
  box(lwd = 2)

  dev.copy(png,paste(paste(name,'_SI','_Model_',data$model,sep=""),"png",sep="."))
  dev.off()
  
  }

