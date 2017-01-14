library('FME')
library('asbio')
library('Hmisc')

Data_All<- as.data.frame(read.csv("~/Desktop/Asthma/Data/Whole_Blood.csv"))
Act_score<- as.data.frame(read.csv("~/Desktop/Asthma/Data/Formatted_ACT2.csv"))


Partition_groups<-function(Data_All,Var_colmn_N)
{

Data_Marker<-Data_All[,c(1,2,Var_colmn_N)]
colnames(Data_Marker)[3]<-c('Marker')

id=0
Data<-Data_Marker[Data_Marker$ID==id,]
B_line<-as.numeric(as.vector(Data$Marker))
D_vals<-B_line[!(is.na(B_line))]

par_est<-ci.median(D_vals, conf = 0.95)
B_median<-as.numeric(par_est$ci[1])
CI_half<-(as.numeric(par_est$ci[3])-as.numeric(par_est$ci[2]))*0.5

B_mean<-mean(D_vals)
Sd_val<-sd(D_vals)

Patient_IDs<-c(5,6,7,9,seq(10,23,by=1),seq(25,50,by=1))
IDs<-Patient_IDs[!duplicated(Patient_IDs)]

Grp1<-numeric(length(IDs))
Grp2<-numeric(length(IDs))
Grp3<-numeric(length(IDs))
Grp4<-numeric(length(IDs))
Grp5<-numeric(length(IDs))
Grp6<-numeric(length(IDs))
Grp7<-numeric(length(IDs))
Grp8<-numeric(length(IDs))

j=1
for(i in IDs)
{
#i=6

Data_vary1<-Data_Marker[Data_Marker$ID==i,]$Marker
Data_time<-Data_Marker[Data_Marker$ID==i,]$Days
#Data_vary<-approx(Data_time,Data_vary1,xout=Data_time,method='linear',n=length(Data_time),rule=1,f=0)$y
Data_vary<-approxExtrap(Data_time,Data_vary1,xout=Data_time,method="linear",n=length(Data_time),rule=2,f=0,na.rm=FALSE)$y

End<-Data_vary[length(Data_vary)]
Start<-Data_vary[1]

#Dff_Start<-abs(Start-B_median)
#Dff_End<-abs(End-B_median)
#Dff_mean<-mean(Data_vary-B_median)
#Sd_val<-CI_half

Dff_Start<-abs(Start-B_mean)
Dff_End<-abs(End-B_mean)
Dff_mean<-mean(Data_vary-B_mean)
Sd_val<-Sd_val

if((Dff_Start>Sd_val)&&(Start>B_mean))
{
  if(Dff_End>Sd_val)
     {
      if(Dff_mean>Sd_val)
          {Grp1[j]<-i} #Non responder
      else{Grp2[j]<-i} #Weak non responder        
     }
  if(Dff_End<Sd_val)
     {
      if(Dff_mean>Sd_val)
          {Grp3[j]<-i} #Weak responder
      else{Grp4[j]<-i} #Responder        
     }
}

if(Dff_Start<=Sd_val)
{
  if(Dff_End>Sd_val)
     {
      if(Dff_mean>Sd_val)
          {Grp5[j]<-i} #Non responder
      else{Grp6[j]<-i} #Indifferent        
      }
  if(Dff_End<Sd_val)
     {
      if(Dff_mean>Sd_val)
          {Grp7[j]<-i} # Weak responder
      else{Grp8[j]<-i} #Indifferent        
     }
}
j=j+1
}

Grp1[Grp1==0] <- NA
Group1<-Grp1[!(is.na(Grp1))]

Grp2[Grp2==0] <- NA
Group2<-Grp2[!(is.na(Grp2))]

Grp3[Grp3==0] <- NA
Group3<-Grp3[!(is.na(Grp3))]

Grp4[Grp4==0] <- NA
Group4<-Grp4[!(is.na(Grp4))]

Grp5[Grp5==0] <- NA
Group5<-Grp5[!(is.na(Grp5))]

Grp6[Grp6==0] <- NA
Group6<-Grp6[!(is.na(Grp6))]

Grp7[Grp7==0] <- NA
Group7<-Grp7[!(is.na(Grp7))]

Grp8[Grp8==0] <- NA
Group8<-Grp8[!(is.na(Grp8))]

Group_Check<-setdiff(IDs,c(Group1,Group2,Group3,Group4,Group5,Group6,Group7,Group8))
if(length(Group_Check)==0)
{
print('All grouped')
}
else
    {print(Group_Check)
     print('No group assigned')}

Selected_Groups<-rbind.fill(data.frame(Group1),data.frame(Group2),data.frame(Group3),data.frame(Group4),data.frame(Group5),data.frame(Group6),data.frame(Group7),data.frame(Group8))

return(Selected_Groups)
}


DF_Cells<-Partition_groups(Data_All,Var_colmn_N=8)
Cell_Data<-Data_All[,c(1,2,Var_colmn_N=8)]

IDA<-c(DF_Cells$Group1[!(is.na(DF_Cells$Group1))])
Plot_Cell_group_patterns(IDA,name='Neutrophils_High-A',Cell_Data,Data_All)

ID2<-c(DF_Cells$Group2[!(is.na(DF_Cells$Group2))])
ID5<-c(DF_Cells$Group5[!(is.na(DF_Cells$Group5))])
ID6<-c(DF_Cells$Group6[!(is.na(DF_Cells$Group6))])
IDB<-c(ID2,ID5,ID6)
Plot_Cell_group_patterns(IDB,name='Neutrophils_High-B',Cell_Data,Data_All)

ID3<-c(DF_Cells$Group3[!(is.na(DF_Cells$Group3))])
ID4<-c(DF_Cells$Group4[!(is.na(DF_Cells$Group4))])
IDC1<-c(ID3,ID4)
#Plot_Cell_group_patterns(ID3,name='Neutrophils_Controlled',Cell_Data,Data_All)

ID7<-c(DF_Cells$Group7[!(is.na(DF_Cells$Group7))])
ID8<-c(DF_Cells$Group8[!(is.na(DF_Cells$Group8))])
IDother<-c(18)

IDC<-c(ID7,ID8,IDother,IDC1)
Plot_Cell_group_patterns(IDC,name='Neutrophils_Normal-C',Cell_Data,Data_All)



Plot_Cell_group_patterns<-function(Group_ID,name,Cell_Data,Data_All)
{

Data_Marker<-Data_All[,c(1,2,Var_colmn_N=8)]
colnames(Data_Marker)[3]<-c('Marker')

id=0
Data<-Data_Marker[Data_Marker$ID==id,]
B_line<-as.numeric(as.vector(Data$Marker))
D_vals<-B_line[!(is.na(B_line))]

par_est<-ci.median(D_vals, conf = 0.95)
B_median<-as.numeric(par_est$ci[1])
CI_half<-(as.numeric(par_est$ci[3])-as.numeric(par_est$ci[2]))*0.5

B_mean<-mean(D_vals)
Sd_val<-sd(D_vals)

#B_mean<-B_median
#Sd_val<-CI_half


path<-paste("~/Desktop/Asthma/Data/Cell_var_",name,".svg")
svg(path,width=9,height=7,pointsize=18)
No_id<-length(Group_ID)
plot(Cell_Data[Cell_Data$ID==5,]$Days,Cell_Data[Cell_Data$ID==5,c(3)],type='n',xlim=c(0,200),ylim=c(0,max(Cell_Data[,3],na.rm=T)/B_mean),xlab='Time in days',ylab=name,cex.lab=1.25,cex.axis=1.25)
mtext(side = 3, text =c('IDs:',Group_ID), at = seq(-15,15+No_id*6,length=length(Group_ID)+1),cex=0.55, line = 0.25,las=1,font=2)

C_center<-seq(-10,208.0,by=0.1)
yCC<-rep((B_mean+Sd_val)/B_mean,length(C_center))
yCW<-rep((B_mean-Sd_val)/B_mean,length(C_center))

polygon(c(C_center, rev(C_center)), c(yCW, rev(yCC)),
        col = "yellow", lty='dashed',lwd=0.1)


k=1
grp_col=sample(rainbow(No_id))
for (i in Group_ID)
{
Input_data<-Cell_Data[Cell_Data$ID==i,]
points(Input_data$Days,Input_data[,3]/B_mean,type='o',pch=16,col='black',cex=1.1,lwd=2.5)
k=k+1
}
abline(h=B_mean/B_mean,lwd=4,lty='dashed',col='firebrick')
abline(h=(B_mean+Sd_val)/B_mean,lwd=4,lty='dashed',col='darkorange')
abline(h=(B_mean-Sd_val)/B_mean,lwd=3,lty='dashed',col='darkorange')
dev.off()
}

