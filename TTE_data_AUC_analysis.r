library('Hmisc')
library('FME')

#Data_CFU<-read.csv("Documents/My_Files/Projects/TTE_is_PK-PD_robust_variable/Data/HRZE_mono_drug_conc.csv")
Data_CFU<-read.csv("~/Documents/Projects/TB/TTE_PK_modeling/HRZE_mono_drug_conc.csv")

errbar(Data_CFU$Hrs,Data_CFU$Mean_Rif,Data_CFU$Mean_Rif+Data_CFU$SD_rif,Data_CFU$Mean_Rif-Data_CFU$SD_rif,bg='grey40',pch=21, ylab='Concentration',xlab='Time',lwd=2,cex=1.35,cex.lab=1.35,font=2,font.lab=2,ylim=c(0,50),xlim=c(0,25))
errbar(Data_CFU$Hrs,Data_CFU$Mean_Inh,Data_CFU$Mean_Inh+Data_CFU$SD_inh,Data_CFU$Mean_Inh-Data_CFU$SD_inh,add=T,bg='orange',pch=21,lwd=2,cex=1.35,cex.lab=1.35,font=2,font.lab=2)
errbar(Data_CFU$Hrs,Data_CFU$Mean_Emb,Data_CFU$Mean_Emb+Data_CFU$SD_emb,Data_CFU$Mean_Emb-Data_CFU$SD_emb,add=T,bg='firebrick',pch=21,lwd=2,cex=1.35,cex.lab=1.35,font=2,font.lab=2)
errbar(Data_CFU$Hrs,Data_CFU$Mean_PZA,Data_CFU$Mean_PZA+Data_CFU$SD_pza,Data_CFU$Mean_PZA-Data_CFU$SD_pza,add=T,bg='cyan',pch=21,lwd=2,cex=1.35,cex.lab=1.35,font=2,font.lab=2)

lines(Data_CFU$Mean_Rif~Data_CFU$Hrs,subset=complete.cases(Data_CFU$Hrs,Data_CFU$Mean_Rif),pch=21,lwd=2,cex=1.35,cex.lab=1.35,font=2,font.lab=2)
lines(Data_CFU$Mean_Inh~Data_CFU$Hrs,subset=complete.cases(Data_CFU$Hrs,Data_CFU$Mean_Inh),pch=21,lwd=2,cex=1.35,cex.lab=1.35,font=2,font.lab=2)
lines(Data_CFU$Mean_Emb~Data_CFU$Hrs,subset=complete.cases(Data_CFU$Hrs,Data_CFU$Mean_Emb),pch=21,lwd=2,cex=1.35,cex.lab=1.35,font=2,font.lab=2)
lines(Data_CFU$Mean_PZA~Data_CFU$Hrs,subset=complete.cases(Data_CFU$Hrs,Data_CFU$Mean_PZA),pch=21,lwd=2,cex=1.35,cex.lab=1.35,font=2,font.lab=2)

#PK model for the 4 drugs in monotherapy
Model_Standard_Drugs_PK<- function(pars_HFS,times,state) {
  derivs <- function(times,state,pars_HFS) { # returns rate of change
    with(as.list(c(state,pars_HFS)), {
      if((24*times)%%24<=1)
      {dose_rif=input_rif
       dose_inh=input_inh
       dose_emb=input_emb
       dose_pza=input_pza}
      else {dose_rif=dose_inh=dose_emb=dose_pza=0.0}
      dDrif=dose_rif-in_rif*Drif
      dDinh=dose_inh-in_inh*Dinh
      dDemb=dose_emb-in_emb*Demb
      dDpza=dose_pza-in_pza*Dpza
      
      dRif=in_rif*Drif-out_rif*Rif
      dInh=in_inh*Dinh-out_inh*Inh
      dEmb=in_emb*Demb-out_emb*Emb
      dPza=in_pza*Dpza-out_pza*PZA
      return(list(c(dDrif,dDinh,dDemb,dDpza,dRif,dInh,dEmb,dPza)))
    })
  }
  out<-lsodes(y=state,times=times,func=derivs,parms=pars_HFS)
  return(out) 
}

Model<-Model_Standard_Drugs_PK
Time<-seq(0,35,by=1/24)
pars_HFS<-list(in_rif=0.25,in_inh=0.25,in_emb=0.25,in_pza=0.25,out_rif=0.25,out_inh=0.25,out_emb=0.25,out_pza=0.25)
#p_control<-coef(Fit_C)
#pars_MAC[c('r','B0')] <- Est_control$C0
Initial_HFS<-c(Drif=600,Dinh=800,Demb=1500,Dpza=300,Rif=0,Inh=0,dEmb=0,Pza=0)
out <- as.data.frame(Model(pars_HFS,Time,Initial_HFS))
#plot(Time,log10(out$B+out$M),type='l',lwd=2,col='red')
plot(Time,log10(out$B+1),type='l',lwd=2,col='red')
lines(Time,log10(out$M+1),type='l',lwd=2,col='red')



