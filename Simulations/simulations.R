# Simulations Section 4

#Loading packages
library(robustbetareg)
library(ggplot2)
library(gridExtra)
library(xtable)
#caminho="/figuras"
source("contamination-functions.R", encoding = "UTF-8")
#setwd(caminho)

X2 <- read.table("Matrix_X.txt", header = FALSE, sep = "", dec = ".")[1:40,]
z=as.matrix(rep(1,length(X2)))
x=as.matrix(cbind(z,X2))

# Scenario A
cr=0.1
B=c(-1,-2)
G=c(5)
N=c(40,80,160,320)
M=1000
est.LSMLE=est.SMLE=est.MLE=est.LMDPDE=est.MDPDE=NULL
pb = txtProgressBar(min = 0, max = M*length(N), initial = 0, style = 3)
stepi=1
tempo.inicio=Sys.time()
for(j in 1:length(N)){
  tempo.amostra=Sys.time()
  for(i in 1:M){
    Seed=Seed+1
    amostra=contaminacaoA(x%*%B,z%*%G,ep=(dim(x)[1]*cr),seed=Seed)
    
    
    y=amostra$sample
    
    est1=tryCatch(robustbetareg(y~x[,2]|1),error=function(e) NULL)#,control = robustbetareg.control(L=0.03)
    if(!is.null(est1)){
      est.LSMLE=rbind(est.LSMLE,c(as.numeric(do.call("c",est1$coefficients)),as.numeric(do.call("c",est1$std.error)),length(x[,2]),est1$Tuning,1,est1$converged))
    }
    
    est2=tryCatch(robustbetareg(y~x[,2]|1,type = c("LMDPDE")),error=function(e) NULL)
    if(!is.null(est2)){
      est.LMDPDE=rbind(est.LMDPDE,c(as.numeric(do.call("c",est2$coefficients)),as.numeric(do.call("c",est2$std.error)),length(x[,2]),est2$Tuning,1,est2$converged))
    }
    
    est3=tryCatch(robustbetareg(y~x[,2]|1,type = c("SMLE")),error=function(e) NULL)
    if(!is.null(est3)){
      est.SMLE=rbind(est.SMLE,c(as.numeric(do.call("c",est3$coefficients)),as.numeric(do.call("c",est3$std.error)),length(x[,2]),est3$Tuning,1,est3$converged))
    }
    
    est4=tryCatch(robustbetareg(y~x[,2]|1,type = c("MDPDE")),error=function(e) NULL)
    if(!is.null(est4)){
      est.MDPDE=rbind(est.MDPDE,c(as.numeric(do.call("c",est4$coefficients)),as.numeric(do.call("c",est4$std.error)),length(x[,2]),est4$Tuning,1,est4$converged))
    }
    
    est5=tryCatch(betareg(y~x[,2]|1),error=function(e) NULL)
    if(!is.null(est5)){
      est.MLE=rbind(est.MLE,c(as.numeric(do.call("c",est5$coefficients)),as.numeric(sqrt(diag(est5$vcov))),length(x[,2]),0,1,1))
    }
    
    y2=amostra$y
    
    est6=tryCatch(robustbetareg(y2~x[,2]|1),error=function(e) NULL)
    if(!is.null(est6)){
      est.LSMLE=rbind(est.LSMLE,c(as.numeric(do.call("c",est6$coefficients)),as.numeric(do.call("c",est6$std.error)),length(x[,2]),est6$Tuning,0,est6$converged))
    }
    
    est7=tryCatch(robustbetareg(y2~x[,2]|1,type = c("LMDPDE")),error=function(e) NULL)
    if(!is.null(est7)){
      est.LMDPDE=rbind(est.LMDPDE,c(as.numeric(do.call("c",est7$coefficients)),as.numeric(do.call("c",est7$std.error)),length(x[,2]),est7$Tuning,0,est7$converged))
    }
    
    est8=tryCatch(robustbetareg(y2~x[,2]|1,type = c("SMLE")),error=function(e) NULL)
    if(!is.null(est8)){
      est.SMLE=rbind(est.SMLE,c(as.numeric(do.call("c",est8$coefficients)),as.numeric(do.call("c",est8$std.error)),length(x[,2]),est8$Tuning,0,est8$converged))
    }
    
    est9=tryCatch(robustbetareg(y2~x[,2]|1,type = c("MDPDE")),error=function(e) NULL)
    if(!is.null(est9)){
      est.MDPDE=rbind(est.MDPDE,c(as.numeric(do.call("c",est9$coefficients)),as.numeric(do.call("c",est9$std.error)),length(x[,2]),est9$Tuning,0,est9$converged))
    }
    
    est10=tryCatch(betareg(y2~x[,2]|1),error=function(e) NULL)
    if(!is.null(est10)){
      est.MLE=rbind(est.MLE,c(as.numeric(do.call("c",est10$coefficients)),as.numeric(sqrt(diag(est10$vcov))),length(x[,2]),0,0,1))
    }
    
    setTxtProgressBar(pb,stepi)
    stepi=stepi+1
    save(est.LSMLE,est.LMDPDE,est.SMLE,est.MDPDE,est.MLE,Seed,file="Est-Robst.RData")
  }
  x=rbind(x,x)
  z=rbind(z,z)
  tempo.amostra=Sys.time()-tempo.amostra
  print(paste0("Tempo amostra N=",N[j]," ",tempo.amostra))
  save(est.LSMLE,est.LMDPDE,est.SMLE,est.MDPDE,est.MLE,Seed,file="Est-Robst.RData")
}
Sys.time()-tempo.inicio

nomesCol=c("B1","B2","G1","ep.B1","ep.B2","ep.G1","N","alpha","Cont","convergence")
est.LSMLE=data.frame(est.LSMLE)
names(est.LSMLE)=nomesCol
est.LSMLE$N <- as.factor(est.LSMLE$N)
est.LSMLE=subset(est.LSMLE,convergence==1)
est.LSMLE$Cont = factor(est.LSMLE$Cont,levels=c(0,1),labels = c("Non-Contaminated","Contaminated"))

est.LMDPDE=data.frame(est.LMDPDE)
names(est.LMDPDE)=nomesCol
est.LMDPDE$N <- as.factor(est.LMDPDE$N)
est.LMDPDE$Cont = factor(est.LMDPDE$Cont,levels=c(0,1),labels = c("Non-Contaminated","Contaminated"))

est.SMLE=data.frame(est.SMLE)
names(est.SMLE)=nomesCol
est.SMLE$N <- as.factor(est.SMLE$N)
est.SMLE$Cont = factor(est.SMLE$Cont,levels=c(0,1),labels = c("Non-Contaminated","Contaminated"))

est.MDPDE=data.frame(est.MDPDE)
names(est.MDPDE)=nomesCol
est.MDPDE$N <- as.factor(est.MDPDE$N)
est.MDPDE$Cont = factor(est.MDPDE$Cont,levels=c(0,1),labels = c("Non-Contaminated","Contaminated"))

est.MLE=data.frame(est.MLE)
names(est.MLE)=nomesCol
est.MLE$N <- as.factor(est.MLE$N)
est.MLE$Cont = factor(est.MLE$Cont,levels=c(0,1),labels = c("Non-Contaminated","Contaminated"))

b1.min=min(est.MLE[,1])
b1.max=max(est.MLE[,1])
b2.min=min(est.MLE[,2])
b2.max=max(est.MLE[,2])
g1.min=min(est.MLE[,3])
g1.max=max(est.MLE[,3])

#LSMLE
pdf(file = paste0(caminho,"/A10-LSMLE-B1",10*G[1],".pdf"),width=6,height=6)
ggplot(est.LSMLE, aes(x=N, y=B1,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=B[1], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(beta[1]),x="Sample Size",title = expression('LSMLE for '~ beta[1]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(b1.min, b1.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/A10-LSMLE-B2",10*G[1],".pdf"),width=6,height=6)
ggplot(est.LSMLE, aes(x=N, y=B2,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=B[2], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(beta[2]),x="Sample Size",title = expression('LSMLE for '~ beta[2]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(b2.min, b2.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/A10-LSMLE-G1",10*G[1],".pdf"),width=6,height=6)
ggplot(est.LSMLE, aes(x=N, y=G1,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=G[1], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(gamma[1]),x="Sample Size",title = expression('LSMLE for '~ gamma[1]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(g1.min, g1.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/A10-LSMLE-Tuning",10*G[1],".pdf"),width=6,height=6)
ggplot(est.LSMLE, aes(x=N, y=alpha,fill=Cont))+geom_boxplot(fill="white")+ labs(fill = "Tuning",y=expression(alpha),x="Sample Size",title = expression('Optimal'~ alpha~'for LSMLE'))+facet_wrap(~Cont)+theme(plot.title=element_text(size=25,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=25),axis.title.y = element_blank(),axis.title.x=element_text(size=25,margin = margin(t=20)))
dev.off()

#SMLE
pdf(file = paste0(caminho,"/A10-SMLE-B1",10*G[1],".pdf"),width=6,height=6)
ggplot(est.SMLE, aes(x=N, y=B1,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=B[1], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(beta[1]),x="Sample Size",title = expression('SMLE for '~ beta[1]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(b1.min, b1.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/A10-SMLE-B2",10*G[1],".pdf"),width=6,height=6)
ggplot(est.SMLE, aes(x=N, y=B2,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=B[2], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(beta[2]),x="Sample Size",title = expression('SMLE for '~ beta[2]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(b2.min, b2.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/A10-SMLE-G1",10*G[1],".pdf"),width=6,height=6)
ggplot(est.SMLE, aes(x=N, y=G1,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=G[1], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(gamma[1]),x="Sample Size",title = expression('SMLE for '~ gamma[1]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(g1.min, g1.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/A10-SMLE-Tuning",10*G[1],".pdf"),width=6,height=6)
ggplot(est.SMLE, aes(x=N, y=alpha,fill=Cont))+geom_boxplot(fill="white")+ labs(fill = "Tuning",y=expression(alpha),x="Sample Size",title = expression('Optimal'~ alpha~'for SMLE'))+facet_wrap(~Cont)+theme(plot.title=element_text(size=25,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=25),axis.title.y = element_blank(),axis.title.x=element_text(size=25,margin = margin(t=20)))
dev.off()

#MLE
pdf(file = paste0(caminho,"/A10-MLE-B1",10*G[1],".pdf"),width=6,height=6)
ggplot(est.MLE, aes(x=N, y=B1,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=B[1], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(beta[1]),x="Sample Size",title = expression('MLE for '~ beta[1]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(b1.min, b1.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/A10-MLE-B2",10*G[1],".pdf"),width=6,height=6)
ggplot(est.MLE, aes(x=N, y=B2,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=B[2], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(beta[2]),x="Sample Size",title = expression('MLE for '~ beta[2]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(b2.min, b2.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/A10-MLE-G1",10*G[1],".pdf"),width=6,height=6)
ggplot(est.MLE, aes(x=N, y=G1,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=G[1], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(gamma[1]),x="Sample Size",title = expression('MLE for '~ gamma[1]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(g1.min, g1.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()

#LMDPDE
pdf(file = paste0(caminho,"/A10-LMDPDE-B1",10*G[1],".pdf"),width=6,height=6)
ggplot(est.LMDPDE, aes(x=N, y=B1,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=B[1], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(beta[1]),x="Sample Size",title = expression('LMDPDE for '~ beta[1]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(b1.min, b1.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/A10-LMDPDE-B2",10*G[1],".pdf"),width=6,height=6)
ggplot(est.LMDPDE, aes(x=N, y=B2,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=B[2], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(beta[2]),x="Sample Size",title = expression('LMDPDE for '~ beta[2]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(b2.min, b2.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/A10-LMDPDE-G1",10*G[1],".pdf"),width=6,height=6)
ggplot(est.LMDPDE, aes(x=N, y=G1,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=G[1], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(gamma[1]),x="Sample Size",title = expression('LMDPDE for '~ gamma[1]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(g1.min, g1.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/A10-LMDPDE-Tuning",10*G[1],".pdf"),width=6,height=6)
ggplot(est.LMDPDE, aes(x=N, y=alpha,fill=Cont))+geom_boxplot(fill="white")+ labs(fill = "Tuning",y=expression(alpha),x="Sample Size",title = expression('Optimal'~ alpha~'for LMDPDE'))+facet_wrap(~Cont)+theme(plot.title=element_text(size=25,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=25),axis.title.y = element_blank(),axis.title.x=element_text(size=25,margin = margin(t=20)))
dev.off()

#MDPDE
pdf(file = paste0(caminho,"/A10-MDPDE-B1",10*G[1],".pdf"),width=6,height=6)
ggplot(est.MDPDE, aes(x=N, y=B1,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=B[1], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(beta[1]),x="Sample Size",title = expression('MDPDE for '~ beta[1]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(b1.min, b1.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/A10-MDPDE-B2",10*G[1],".pdf"),width=6,height=6)
ggplot(est.MDPDE, aes(x=N, y=B2,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=B[2], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(beta[2]),x="Sample Size",title = expression('MDPDE for '~ beta[2]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(b2.min, b2.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/A10-MDPDE-G1",10*G[1],".pdf"),width=6,height=6)
ggplot(est.MDPDE, aes(x=N, y=G1,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=G[1], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(gamma[1]),x="Sample Size",title = expression('MDPDE for '~ gamma[1]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(g1.min, g1.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/A10-MDPDE-Tuning",10*G[1],".pdf"),width=6,height=6)
ggplot(est.MDPDE, aes(x=N, y=alpha,fill=Cont))+geom_boxplot(fill="white")+ labs(fill = "Tuning",y=expression(alpha),x="Sample Size",title = expression('Optimal'~ alpha~'for MDPDE'))+facet_wrap(~Cont)+theme(plot.title=element_text(size=25,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=25),axis.title.y = element_blank(),axis.title.x=element_text(size=25,margin = margin(t=20)))
dev.off()


#Saving
z=as.matrix(rep(1,length(X2)))
x=cbind(z,X2)
parametros=list(B=B,G=G,Seed=2022,x=x,z=z)
save(est.LMDPDE,est.LSMLE,est.MDPDE,est.SMLE,est.MLE,parametros,file = "SimCont-CenarioA10.RData")
getwd()

# Scenario B

X2 <- read.table("Matrix_X.txt", header = FALSE, sep = "", dec = ".")[1:40,]
z=as.matrix(rep(1,length(X2)))
x=as.matrix(cbind(z,X2))

cr=0.05
B=c(-1.0,-5.5)
G=c(5.0)
N=c(40,80,160,320)
M=1000
est.LSMLE=est.SMLE=est.MLE=est.LMDPDE=est.MDPDE=NULL
pb = txtProgressBar(min = 0, max = M*length(N), initial = 0, style = 3)
stepi=1
temp=NULL
tempo.inicio=Sys.time()
for(j in 1:length(N)){
  tempo.amostra=Sys.time()
  for(i in 1:M){
    Seed=Seed+1
    amostra=contaminacaoB(x%*%B,z%*%G,ep=(dim(x)[1]*cr),seed=Seed,x=x)
    
    y=amostra$sample 
    
    est1=tryCatch(robustbetareg(y~x[,2]|1),error=function(e) NULL)#,control = robustbetareg.control(L=0.03)
    if(!is.null(est1)){
      est.LSMLE=rbind(est.LSMLE,c(as.numeric(do.call("c",est1$coefficients)),as.numeric(do.call("c",est1$std.error)),length(x[,2]),est1$Tuning,1,est1$converged))
      temp=c(temp,est1$Tuning)
    }
    
    est2=tryCatch(robustbetareg(y~x[,2]|1,type = c("LMDPDE")),error=function(e) NULL)
    if(!is.null(est2)){
      est.LMDPDE=rbind(est.LMDPDE,c(as.numeric(do.call("c",est2$coefficients)),as.numeric(do.call("c",est2$std.error)),length(x[,2]),est2$Tuning,1,est2$converged))
    }
    
    est3=tryCatch(robustbetareg(y~x[,2]|1,type = c("SMLE")),error=function(e) NULL)
    if(!is.null(est3)){
      est.SMLE=rbind(est.SMLE,c(as.numeric(do.call("c",est3$coefficients)),as.numeric(do.call("c",est3$std.error)),length(x[,2]),est3$Tuning,1,est3$converged))
    }
    
    est4=tryCatch(robustbetareg(y~x[,2]|1,type = c("MDPDE")),error=function(e) NULL)
    if(!is.null(est4)){
      est.MDPDE=rbind(est.MDPDE,c(as.numeric(do.call("c",est4$coefficients)),as.numeric(do.call("c",est4$std.error)),length(x[,2]),est4$Tuning,1,est4$converged))
    }
    
    est5=tryCatch(betareg(y~x[,2]|1),error=function(e) NULL)
    if(!is.null(est5)){
      est.MLE=rbind(est.MLE,c(as.numeric(do.call("c",est5$coefficients)),as.numeric(sqrt(diag(est5$vcov))),length(x[,2]),0,1,1))
    }
    
    ########
    
    y2=amostra$y 
    
    est6=tryCatch(robustbetareg(y2~x[,2]|1),error=function(e) NULL)
    if(!is.null(est6)){
      est.LSMLE=rbind(est.LSMLE,c(as.numeric(do.call("c",est6$coefficients)),as.numeric(do.call("c",est6$std.error)),length(x[,2]),est6$Tuning,0,est6$converged))
    }
    
    est7=tryCatch(robustbetareg(y2~x[,2]|1,type = c("LMDPDE")),error=function(e) NULL)
    if(!is.null(est7)){
      est.LMDPDE=rbind(est.LMDPDE,c(as.numeric(do.call("c",est7$coefficients)),as.numeric(do.call("c",est7$std.error)),length(x[,2]),est7$Tuning,0,est7$converged))
    }
    
    est8=tryCatch(robustbetareg(y2~x[,2]|1,type = c("SMLE")),error=function(e) NULL)
    if(!is.null(est8)){
      est.SMLE=rbind(est.SMLE,c(as.numeric(do.call("c",est8$coefficients)),as.numeric(do.call("c",est8$std.error)),length(x[,2]),est8$Tuning,0,est8$converged))
    }
    
    est9=tryCatch(robustbetareg(y2~x[,2]|1,type = c("MDPDE")),error=function(e) NULL)
    if(!is.null(est9)){
      est.MDPDE=rbind(est.MDPDE,c(as.numeric(do.call("c",est9$coefficients)),as.numeric(do.call("c",est9$std.error)),length(x[,2]),est9$Tuning,0,est9$converged))
    }
    
    est10=tryCatch(betareg(y2~x[,2]|1),error=function(e) NULL)
    if(!is.null(est10)){
      est.MLE=rbind(est.MLE,c(as.numeric(do.call("c",est10$coefficients)),as.numeric(sqrt(diag(est10$vcov))),length(x[,2]),0,0,1))
    }
    
    setTxtProgressBar(pb,stepi)
    stepi=stepi+1
    save(est.LSMLE,est.LMDPDE,est.SMLE,est.MDPDE,est.MLE,Seed,file="Est-Robst.RData")
  }
  x=rbind(x,x)
  z=rbind(z,z)
  tempo.amostra=Sys.time()-tempo.amostra
  print(paste0("Tempo amostra N=",N[j]," ",tempo.amostra))
  save(est.LSMLE,est.LMDPDE,est.SMLE,est.MDPDE,est.MLE,Seed,file="Est-Robst.RData")
}

nomesCol=c("B1","B2","G1","ep.B1","ep.B2","ep.G1","N","alpha","Cont","convergence")
est.LSMLE=data.frame(est.LSMLE)
names(est.LSMLE)=nomesCol
est.LSMLE$N <- as.factor(est.LSMLE$N)
est.LSMLE=subset(est.LSMLE,convergence==1)
est.LSMLE$Cont = factor(est.LSMLE$Cont,levels=c(0,1),labels = c("Non-Contaminated","Contaminated"))

est.LMDPDE=data.frame(est.LMDPDE)
names(est.LMDPDE)=nomesCol
est.LMDPDE$N <- as.factor(est.LMDPDE$N)
est.LMDPDE$Cont = factor(est.LMDPDE$Cont,levels=c(0,1),labels = c("Non-Contaminated","Contaminated"))

est.SMLE=data.frame(est.SMLE)
names(est.SMLE)=nomesCol
est.SMLE$N <- as.factor(est.SMLE$N)
est.SMLE$Cont = factor(est.SMLE$Cont,levels=c(0,1),labels = c("Non-Contaminated","Contaminated"))

est.MDPDE=data.frame(est.MDPDE)
names(est.MDPDE)=nomesCol
est.MDPDE$N <- as.factor(est.MDPDE$N)
est.MDPDE$Cont = factor(est.MDPDE$Cont,levels=c(0,1),labels = c("Non-Contaminated","Contaminated"))

est.MLE=data.frame(est.MLE)
names(est.MLE)=nomesCol
est.MLE$N <- as.factor(est.MLE$N)
est.MLE$Cont = factor(est.MLE$Cont,levels=c(0,1),labels = c("Non-Contaminated","Contaminated"))


b1.min=min(est.MLE[,1])
b1.max=max(est.MLE[,1])
b2.min=min(est.MLE[,2])
b2.max=max(est.MLE[,2])
g1.min=min(est.MLE[,3])
g1.max=max(est.MLE[,3])

#LSMLE
pdf(file = paste0(caminho,"/B10-LSMLE-B1",10*G[1],".pdf"),width=6,height=6)
ggplot(est.LSMLE, aes(x=N, y=B1,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=B[1], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(beta[1]),x="Sample Size",title = expression('LSMLE for '~ beta[1]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(b1.min, b1.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/B10-LSMLE-B2",10*G[1],".pdf"),width=6,height=6)
ggplot(est.LSMLE, aes(x=N, y=B2,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=B[2], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(beta[2]),x="Sample Size",title = expression('LSMLE for '~ beta[2]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(b2.min, b2.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/B10-LSMLE-G1",10*G[1],".pdf"),width=6,height=6)
ggplot(est.LSMLE, aes(x=N, y=G1,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=G[1], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(gamma[1]),x="Sample Size",title = expression('LSMLE for '~ gamma[1]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(g1.min, g1.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/B10-LSMLE-Tuning",10*G[1],".pdf"),width=6,height=6)
ggplot(est.LSMLE, aes(x=N, y=alpha,fill=Cont))+geom_boxplot(fill="white")+ labs(fill = "Tuning",y=expression(alpha),x="Sample Size",title = expression('Optimal'~ alpha~'for LSMLE'))+facet_wrap(~Cont)+theme(plot.title=element_text(size=25,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=25),axis.title.y = element_blank(),axis.title.x=element_text(size=25,margin = margin(t=20)))
dev.off()

#SMLE
pdf(file = paste0(caminho,"/B10-SMLE-B1",10*G[1],".pdf"),width=6,height=6)
ggplot(est.SMLE, aes(x=N, y=B1,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=B[1], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(beta[1]),x="Sample Size",title = expression('SMLE for '~ beta[1]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(b1.min, b1.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/B10-SMLE-B2",10*G[1],".pdf"),width=6,height=6)
ggplot(est.SMLE, aes(x=N, y=B2,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=B[2], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(beta[2]),x="Sample Size",title = expression('SMLE for '~ beta[2]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(b2.min, b2.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/B10-SMLE-G1",10*G[1],".pdf"),width=6,height=6)
ggplot(est.SMLE, aes(x=N, y=G1,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=G[1], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(gamma[1]),x="Sample Size",title = expression('SMLE for '~ gamma[1]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(g1.min, g1.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/B10-SMLE-Tuning",10*G[1],".pdf"),width=6,height=6)
ggplot(est.SMLE, aes(x=N, y=alpha,fill=Cont))+geom_boxplot(fill="white")+ labs(fill = "Tuning",y=expression(alpha),x="Sample Size",title = expression('Optimal'~ alpha~'for SMLE'))+facet_wrap(~Cont)+theme(plot.title=element_text(size=25,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=25),axis.title.y = element_blank(),axis.title.x=element_text(size=25,margin = margin(t=20)))
dev.off()

#MLE
pdf(file = paste0(caminho,"/B10-MLE-B1",10*G[1],".pdf"),width=6,height=6)
ggplot(est.MLE, aes(x=N, y=B1,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=B[1], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(beta[1]),x="Sample Size",title = expression('MLE for '~ beta[1]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(b1.min, b1.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/B10-MLE-B2",10*G[1],".pdf"),width=6,height=6)
ggplot(est.MLE, aes(x=N, y=B2,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=B[2], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(beta[2]),x="Sample Size",title = expression('MLE for '~ beta[2]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(b2.min, b2.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/B10-MLE-G1",10*G[1],".pdf"),width=6,height=6)
ggplot(est.MLE, aes(x=N, y=G1,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=G[1], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(gamma[1]),x="Sample Size",title = expression('MLE for '~ gamma[1]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(g1.min, g1.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()

#LMDPDE
pdf(file = paste0(caminho,"/B10-LMDPDE-B1",10*G[1],".pdf"),width=6,height=6)
ggplot(est.LMDPDE, aes(x=N, y=B1,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=B[1], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(beta[1]),x="Sample Size",title = expression('LMDPDE for '~ beta[1]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(b1.min, b1.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/B10-LMDPDE-B2",10*G[1],".pdf"),width=6,height=6)
ggplot(est.LMDPDE, aes(x=N, y=B2,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=B[2], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(beta[2]),x="Sample Size",title = expression('LMDPDE for '~ beta[2]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(b2.min, b2.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/B10-LMDPDE-G1",10*G[1],".pdf"),width=6,height=6)
ggplot(est.LMDPDE, aes(x=N, y=G1,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=G[1], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(gamma[1]),x="Sample Size",title = expression('LMDPDE for '~ gamma[1]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(g1.min, g1.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/B10-LMDPDE-Tuning",10*G[1],".pdf"),width=6,height=6)
ggplot(est.LMDPDE, aes(x=N, y=alpha,fill=Cont))+geom_boxplot(fill="white")+ labs(fill = "Tuning",y=expression(alpha),x="Sample Size",title = expression('Optimal'~ alpha~'for LMDPDE'))+facet_wrap(~Cont)+theme(plot.title=element_text(size=25,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=25),axis.title.y = element_blank(),axis.title.x=element_text(size=25,margin = margin(t=20)))
dev.off()

#MDPDE
pdf(file = paste0(caminho,"/B10-MDPDE-B1",10*G[1],".pdf"),width=6,height=6)
ggplot(est.MDPDE, aes(x=N, y=B1,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=B[1], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(beta[1]),x="Sample Size",title = expression('MDPDE for '~ beta[1]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(b1.min, b1.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/B10-MDPDE-B2",10*G[1],".pdf"),width=6,height=6)
ggplot(est.MDPDE, aes(x=N, y=B2,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=B[2], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(beta[2]),x="Sample Size",title = expression('MDPDE for '~ beta[2]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(b2.min, b2.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/B10-MDPDE-G1",10*G[1],".pdf"),width=6,height=6)
ggplot(est.MDPDE, aes(x=N, y=G1,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=G[1], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(gamma[1]),x="Sample Size",title = expression('MDPDE for '~ gamma[1]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(g1.min, g1.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/B10-MDPDE-Tuning",10*G[1],".pdf"),width=6,height=6)
ggplot(est.MDPDE, aes(x=N, y=alpha,fill=Cont))+geom_boxplot(fill="white")+ labs(fill = "Tuning",y=expression(alpha),x="Sample Size",title = expression('Optimal'~ alpha~'for MDPDE'))+facet_wrap(~Cont)+theme(plot.title=element_text(size=25,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=25),axis.title.y = element_blank(),axis.title.x=element_text(size=25,margin = margin(t=20)))
dev.off()

#Saving
z=as.matrix(rep(1,length(X2)))
x=cbind(z,X2)
parametros=list(B=B,G=G,Seed=2022,x=x,z=z)
save(est.LMDPDE,est.LSMLE,est.MDPDE,est.SMLE,est.MLE,parametros,file = "SimCont-CenarioB10.RData")
save(est.LMDPDE,est.LSMLE,est.MDPDE,est.SMLE,est.MLE,parametros,file = "SimCont-B.RData")


# Scenario C

X2 <- read.table("Matrix_X.txt", header = FALSE, sep = "", dec = ".")[1:40,]
x=z=cbind(rep(1,length(X2)),X2)

cr=0.05
B=c(-3,7.5)
G=c(1,2)
N=c(40,80,160,320)
M=1000
est.LSMLE=est.SMLE=est.MLE=est.LMDPDE=est.MDPDE=NULL
pb = txtProgressBar(min = 0, max = M*length(N), initial = 0, style = 3)
stepi=1
temp=NULL
tempo.inicio=Sys.time()
for(j in 1:length(N)){
  tempo.amostra=Sys.time()
  for(i in 1:M){
    Seed=Seed+1
    amostra=contaminacaoC(x%*%B,z%*%G,ep=(dim(x)[1]*cr),seed=Seed)
    
    y=amostra$sample 
    
    est1=tryCatch(robustbetareg(y~x[,2]|x[,2]),error=function(e) NULL)#,control = robustbetareg.control(L=0.03)
    if(!is.null(est1)){
      est.LSMLE=rbind(est.LSMLE,c(as.numeric(do.call("c",est1$coefficients)),as.numeric(do.call("c",est1$std.error)),length(x[,2]),est1$Tuning,1,est1$converged))
      temp=c(temp,est1$Tuning)
    }
    
    est2=tryCatch(robustbetareg(y~x[,2]|x[,2],type = c("LMDPDE")),error=function(e) NULL)
    if(!is.null(est2)){
      est.LMDPDE=rbind(est.LMDPDE,c(as.numeric(do.call("c",est2$coefficients)),as.numeric(do.call("c",est2$std.error)),length(x[,2]),est2$Tuning,1,est2$converged))
    }
    
    est3=tryCatch(robustbetareg(y~x[,2]|x[,2],type = c("SMLE")),error=function(e) NULL)
    if(!is.null(est3)){
      est.SMLE=rbind(est.SMLE,c(as.numeric(do.call("c",est3$coefficients)),as.numeric(do.call("c",est3$std.error)),length(x[,2]),est3$Tuning,1,est3$converged))
    }
    
    est4=tryCatch(robustbetareg(y~x[,2]|x[,2],type = c("MDPDE")),error=function(e) NULL)
    if(!is.null(est4)){
      est.MDPDE=rbind(est.MDPDE,c(as.numeric(do.call("c",est4$coefficients)),as.numeric(do.call("c",est4$std.error)),length(x[,2]),est4$Tuning,1,est4$converged))
    }
    
    est5=tryCatch(betareg(y~x[,2]|x[,2]),error=function(e) NULL)
    if(!is.null(est5)){
      est.MLE=rbind(est.MLE,c(as.numeric(do.call("c",est5$coefficients)),as.numeric(sqrt(diag(est5$vcov))),length(x[,2]),0,1,1))
    }
    
    y2=amostra$y 
    
    est6=tryCatch(robustbetareg(y2~x[,2]|x[,2]),error=function(e) NULL)
    if(!is.null(est6)){
      est.LSMLE=rbind(est.LSMLE,c(as.numeric(do.call("c",est6$coefficients)),as.numeric(do.call("c",est6$std.error)),length(x[,2]),est6$Tuning,0,est6$converged))
    }
    
    est7=tryCatch(robustbetareg(y2~x[,2]|x[,2],type = c("LMDPDE")),error=function(e) NULL)
    if(!is.null(est7)){
      est.LMDPDE=rbind(est.LMDPDE,c(as.numeric(do.call("c",est7$coefficients)),as.numeric(do.call("c",est7$std.error)),length(x[,2]),est7$Tuning,0,est7$converged))
    }
    
    est8=tryCatch(robustbetareg(y2~x[,2]|x[,2],type = c("SMLE")),error=function(e) NULL)
    if(!is.null(est8)){
      est.SMLE=rbind(est.SMLE,c(as.numeric(do.call("c",est8$coefficients)),as.numeric(do.call("c",est8$std.error)),length(x[,2]),est8$Tuning,0,est8$converged))
    }
    
    est9=tryCatch(robustbetareg(y2~x[,2]|x[,2],type = c("MDPDE")),error=function(e) NULL)
    if(!is.null(est9)){
      est.MDPDE=rbind(est.MDPDE,c(as.numeric(do.call("c",est9$coefficients)),as.numeric(do.call("c",est9$std.error)),length(x[,2]),est9$Tuning,0,est9$converged))
    }
    
    est10=tryCatch(betareg(y2~x[,2]|x[,2]),error=function(e) NULL)
    if(!is.null(est10)){
      est.MLE=rbind(est.MLE,c(as.numeric(do.call("c",est10$coefficients)),as.numeric(sqrt(diag(est10$vcov))),length(x[,2]),0,0,1))
    }
    
    setTxtProgressBar(pb,stepi)
    stepi=stepi+1
    save(est.LSMLE,est.LMDPDE,est.SMLE,est.MDPDE,est.MLE,Seed,file="Est-Robst.RData")
  }
  x=rbind(x,x)
  z=rbind(z,z)
  tempo.amostra=Sys.time()-tempo.amostra
  print(paste0("Tempo amostra N=",N[j]," ",tempo.amostra))
  save(est.LSMLE,est.LMDPDE,est.SMLE,est.MDPDE,est.MLE,Seed,file="Est-Robst.RData")
}
Sys.time()-tempo.inicio

nomesCol=c("B1","B2","G1","G2","ep.B1","ep.B2","ep.G1","ep.G2","N","alpha","Cont","convergence")
est.LSMLE=data.frame(est.LSMLE)
names(est.LSMLE)=nomesCol
est.LSMLE$N <- as.factor(est.LSMLE$N)
est.LSMLE=subset(est.LSMLE,convergence==1)
est.LSMLE$Cont = factor(est.LSMLE$Cont,levels=c(0,1),labels = c("Non-Contaminated","Contaminated"))

est.LMDPDE=data.frame(est.LMDPDE)
names(est.LMDPDE)=nomesCol
est.LMDPDE$N <- as.factor(est.LMDPDE$N)
est.LMDPDE$Cont = factor(est.LMDPDE$Cont,levels=c(0,1),labels = c("Non-Contaminated","Contaminated"))

est.SMLE=data.frame(est.SMLE)
names(est.SMLE)=nomesCol
est.SMLE$N <- as.factor(est.SMLE$N)
est.SMLE$Cont = factor(est.SMLE$Cont,levels=c(0,1),labels = c("Non-Contaminated","Contaminated"))

est.MDPDE=data.frame(est.MDPDE)
names(est.MDPDE)=nomesCol
est.MDPDE$N <- as.factor(est.MDPDE$N)
est.MDPDE$Cont = factor(est.MDPDE$Cont,levels=c(0,1),labels = c("Non-Contaminated","Contaminated"))

est.MLE=data.frame(est.MLE)
names(est.MLE)=nomesCol
est.MLE$N <- as.factor(est.MLE$N)
est.MLE$Cont = factor(est.MLE$Cont,levels=c(0,1),labels = c("Non-Contaminated","Contaminated"))

delta=2
b1.min=min(est.MLE[,1])-delta
b1.max=max(est.MLE[,1])+delta
b2.min=min(est.MLE[,2])-delta
b2.max=max(est.MLE[,2])+delta
g1.min=min(est.MLE[,3])-delta
g1.max=max(est.MLE[,3])+delta/2
g2.min=min(est.MLE[,4])-delta
g2.max=max(est.MLE[,4])+delta

#LSMLE
pdf(file = paste0(caminho,"/C-LSMLE-B1.pdf"),width=6,height=6)
ggplot(est.LSMLE, aes(x=N, y=B1,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=B[1], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(beta[1]),x="Sample Size",title = expression('LSMLE for '~ beta[1]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(b1.min, b1.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/C-LSMLE-B2.pdf"),width=6,height=6)
ggplot(est.LSMLE, aes(x=N, y=B2,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=B[2], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(beta[2]),x="Sample Size",title = expression('LSMLE for '~ beta[2]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(b2.min, b2.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/C-LSMLE-G1.pdf"),width=6,height=6)
ggplot(est.LSMLE, aes(x=N, y=G1,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=G[1], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(gamma[1]),x="Sample Size",title = expression('LSMLE for '~ gamma[1]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(g1.min, g1.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/C-LSMLE-G2.pdf"),width=6,height=6)
ggplot(est.LSMLE, aes(x=N, y=G2,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=G[2], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(gamma[2]),x="Sample Size",title = expression('LSMLE for '~ gamma[2]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(g2.min, g2.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/C-LSMLE-Tuning.pdf"),width=6,height=6)
ggplot(est.LSMLE, aes(x=N, y=alpha,fill=Cont))+geom_boxplot(fill="white")+ labs(fill = "Tuning",y=expression(alpha),x="Sample Size",title = expression('Optimal'~ alpha~'for LSMLE'))+facet_wrap(~Cont)+theme(plot.title=element_text(size=25,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=25),axis.title.y = element_blank(),axis.title.x=element_text(size=25,margin = margin(t=20)))
dev.off()

#SMLE
pdf(file = paste0(caminho,"/C-SMLE-B1.pdf"),width=6,height=6)
ggplot(est.SMLE, aes(x=N, y=B1,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=B[1], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(beta[1]),x="Sample Size",title = expression('SMLE for '~ beta[1]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(b1.min, b1.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/C-SMLE-B2.pdf"),width=6,height=6)
ggplot(est.SMLE, aes(x=N, y=B2,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=B[2], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(beta[2]),x="Sample Size",title = expression('SMLE for '~ beta[2]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(b2.min, b2.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/C-SMLE-G1.pdf"),width=6,height=6)
ggplot(est.SMLE, aes(x=N, y=G1,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=G[1], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(gamma[1]),x="Sample Size",title = expression('SMLE for '~ gamma[1]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(g1.min, g1.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/C-SMLE-G2.pdf"),width=6,height=6)
ggplot(est.SMLE, aes(x=N, y=G2,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=G[2], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(gamma[2]),x="Sample Size",title = expression('SMLE for '~ gamma[2]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(g2.min, g2.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/C-SMLE-Tuning.pdf"),width=6,height=6)
ggplot(est.SMLE, aes(x=N, y=alpha,fill=Cont))+geom_boxplot(fill="white")+ labs(fill = "Tuning",y=expression(alpha),x="Sample Size",title = expression('Optimal'~ alpha~'for SMLE'))+facet_wrap(~Cont)+theme(plot.title=element_text(size=25,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=25),axis.title.y = element_blank(),axis.title.x=element_text(size=25,margin = margin(t=20)))
dev.off()

#MLE
pdf(file = paste0(caminho,"/C-MLE-B1.pdf"),width=6,height=6)
ggplot(est.MLE, aes(x=N, y=B1,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=B[1], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(beta[1]),x="Sample Size",title = expression('MLE for '~ beta[1]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(b1.min, b1.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/C-MLE-B2.pdf"),width=6,height=6)
ggplot(est.SMLE, aes(x=N, y=B2,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=B[2], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(beta[2]),x="Sample Size",title = expression('MLE for '~ beta[2]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(b2.min, b2.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/C-MLE-G1.pdf"),width=6,height=6)
ggplot(est.MLE, aes(x=N, y=G1,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=G[1], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(gamma[1]),x="Sample Size",title = expression('MLE for '~ gamma[1]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(g1.min, g1.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/C-MLE-G2.pdf"),width=6,height=6)
ggplot(est.MLE, aes(x=N, y=G2,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=G[2], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(gamma[2]),x="Sample Size",title = expression('MLE for '~ gamma[2]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(g2.min, g2.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()

#LMDPDE
pdf(file = paste0(caminho,"/C-LMDPDE-B1.pdf"),width=6,height=6)
ggplot(est.LMDPDE, aes(x=N, y=B1,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=B[1], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(beta[1]),x="Sample Size",title = expression('LMDPDE for '~ beta[1]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(b1.min, b1.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/C-LMDPDE-B2.pdf"),width=6,height=6)
ggplot(est.LMDPDE, aes(x=N, y=B2,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=B[2], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(beta[2]),x="Sample Size",title = expression('LMDPDE for '~ beta[2]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(b2.min, b2.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/C-LMDPDE-G1.pdf"),width=6,height=6)
ggplot(est.LMDPDE, aes(x=N, y=G1,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=G[1], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(gamma[1]),x="Sample Size",title = expression('LMDPDE for '~ gamma[1]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(g1.min, g1.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/C-LMDPDE-G2.pdf"),width=6,height=6)
ggplot(est.LMDPDE, aes(x=N, y=G2,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=G[2], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(gamma[2]),x="Sample Size",title = expression('LMDPDE for '~ gamma[2]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(g2.min, g2.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/C-LMDPDE-Tuning.pdf"),width=6,height=6)
ggplot(est.LMDPDE, aes(x=N, y=alpha,fill=Cont))+geom_boxplot(fill="white")+ labs(fill = "Tuning",y=expression(alpha),x="Sample Size",title = expression('Optimal'~ alpha~'for LMDPDE'))+facet_wrap(~Cont)+theme(plot.title=element_text(size=25,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=25),axis.title.y = element_blank(),axis.title.x=element_text(size=25,margin = margin(t=20)))
dev.off()

#MDPDE
pdf(file = paste0(caminho,"/C-MDPDE-B1.pdf"),width=6,height=6)
ggplot(est.MDPDE, aes(x=N, y=B1,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=B[1], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(beta[1]),x="Sample Size",title = expression('MDPDE for '~ beta[1]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(b1.min, b1.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/C-MDPDE-B2.pdf"),width=6,height=6)
ggplot(est.MDPDE, aes(x=N, y=B2,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=B[2], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(beta[2]),x="Sample Size",title = expression('MDPDE for '~ beta[2]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(b2.min, b2.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/C-MDPDE-G1.pdf"),width=6,height=6)
ggplot(est.MDPDE, aes(x=N, y=G1,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=G[1], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(gamma[1]),x="Sample Size",title = expression('MDPDE for '~ gamma[1]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(g1.min, g1.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/C-MDPDE-G2.pdf"),width=6,height=6)
ggplot(est.MDPDE, aes(x=N, y=G2,fill=Cont))+geom_boxplot(fill="white")+geom_hline(yintercept=G[2], linetype="dashed", color = "red")+ labs(fill = "Tuning",y=expression(gamma[2]),x="Sample Size",title = expression('MDPDE for '~ gamma[2]))+facet_wrap(~Cont)+coord_cartesian(ylim = c(g2.min, g2.max))+theme(plot.title=element_text(size=20,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=18),axis.title.y = element_blank(),axis.title.x=element_text(size=20,margin = margin(t=20)))
dev.off()
pdf(file = paste0(caminho,"/C-MDPDE-Tuning.pdf"),width=6,height=6)
ggplot(est.MDPDE, aes(x=N, y=alpha,fill=Cont))+geom_boxplot(fill="white")+ labs(fill = "Tuning",y=expression(alpha),x="Sample Size",title = expression('Optimal'~ alpha~'for MDPDE'))+facet_wrap(~Cont)+theme(plot.title=element_text(size=25,hjust=0.5), strip.text.x = element_text(size = 19),axis.text=element_text(size=25),axis.title.y = element_blank(),axis.title.x=element_text(size=25,margin = margin(t=20)))
dev.off()

x=z=cbind(rep(1,length(X2)),X2)
parametros=list(B=B,G=G,Seed=2022,x=x,z=z)
save(est.LMDPDE,est.LSMLE,est.MDPDE,est.SMLE,est.MLE,parametros,file = "SimCont-CenarioC.RData")
save(est.LMDPDE,est.LSMLE,est.MDPDE,est.SMLE,est.MLE,parametros,file = "SimCont-C.RData")


## Simulation for Figure 2
# load the functions in RobustEstimators.R script

# Scenario B

sim.convergencia = function(M, n){
  
  set.seed(2021 + 04*01)
  
  X = matrix( c( rep(1, 40), runif(40)), ncol = 2, byrow = FALSE )
  
  if(n == 80){
    X <- rbind(X, X)
  }
  par.mu = c(-1, -5.5)
  pred.mu = X%*%par.mu
  mu = exp(pred.mu)/(1+exp(pred.mu)) 
  
  Z = matrix( c( rep(1,n)), ncol = 1, byrow = FALSE )
  par.phi = 5
  pred.phi = Z%*%par.phi
  phi = exp(pred.phi) 
  
  
  alpha <- seq(0, 1, 0.05)
  
  conv.SMLE   <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  conv.LSMLE  <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  conv.MDPDE  <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  conv.LMDPDE <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  
  conv.SMLE.cont   <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  conv.LSMLE.cont  <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  conv.MDPDE.cont  <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  conv.LMDPDE.cont <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  
  
  NAsd.SMLE   <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  NAsd.LSMLE  <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  NAsd.MDPDE  <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  NAsd.LMDPDE <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  
  NAsd.SMLE.cont   <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  NAsd.LSMLE.cont  <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  NAsd.MDPDE.cont  <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  NAsd.LMDPDE.cont <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  
  Total.SMLE   <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  Total.LSMLE  <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  Total.MDPDE  <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  Total.LMDPDE <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  
  Total.SMLE.cont   <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  Total.LSMLE.cont  <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  Total.MDPDE.cont  <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  Total.LMDPDE.cont <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  
  i = 1
  k = 0
  
  while(i <= M){
    tryCatch({
      k = k + 1
      set.seed(3*k+2021)
      y  = rbeta(n, mu*phi, (1-mu)*phi)
      
      #contaminando os dados
      y.cont <- y
      m <- round(n*0.05)
      amos1 <- order(X[,2])[1:m]
      y.cont[amos1] <- rbeta(m, 0.002*phi[amos1], (1 - 0.002)*phi[amos1])
      
      
      for(l in 1:length(alpha)){
        # SMLE
        fit.SMLE <- SMLE_BETA(y, X, Z, qoptimal = FALSE, q0 = 1 - alpha[l])
        fit.SMLE.cont <- SMLE_BETA(y.cont, X, Z, qoptimal = FALSE, q0 = 1 - alpha[l])
        # LSMLE
        fit.LSMLE <- LSMLE_BETA(y, X, Z, qoptimal = FALSE, q0 = 1 - alpha[l])
        fit.LSMLE.cont <- LSMLE_BETA(y.cont, X, Z, qoptimal = FALSE, q0 = 1 - alpha[l])
        # MDPDE
        fit.MDPDE <- MDPDE_BETA(y, X, Z, qoptimal = FALSE, q0 = 1 - alpha[l])
        fit.MDPDE.cont <- MDPDE_BETA(y.cont, X, Z, qoptimal = FALSE, q0 = 1 - alpha[l])
        # LMDPDE
        fit.LMDPDE <- LMDPDE_BETA(y, X, Z, qoptimal = FALSE, q0 = 1 - alpha[l])
        fit.LMDPDE.cont <- LMDPDE_BETA(y.cont, X, Z, qoptimal = FALSE, q0 = 1 - alpha[l])
        
        
        conv.SMLE[l, i]        <- ifelse(is.null(fit.SMLE$Warning), 1, 0) # 1 conv; 0 n?o convergiu
        conv.SMLE.cont[l, i]   <- ifelse(is.null(fit.SMLE.cont$Warning), 1, 0) 
        conv.LSMLE[l, i]       <- ifelse(is.null(fit.LSMLE$Warning), 1, 0) 
        conv.LSMLE.cont[l, i]  <- ifelse(is.null(fit.LSMLE.cont$Warning), 1, 0) 
        conv.MDPDE[l, i]       <- ifelse(is.null(fit.MDPDE$Warning), 1, 0) 
        conv.MDPDE.cont[l, i]  <- ifelse(is.null(fit.MDPDE.cont$Warning), 1, 0) 
        conv.LMDPDE[l, i]      <- ifelse(is.null(fit.LMDPDE$Warning), 1, 0) 
        conv.LMDPDE.cont[l, i] <- ifelse(is.null(fit.LMDPDE.cont$Warning), 1, 0) 
        
        NAsd.SMLE[l, i]        <- ifelse(anyNA(fit.SMLE[[3]], fit.SMLE[[4]]), 0, 1) # 1 funcinou; 0 erro
        NAsd.SMLE.cont[l, i]   <- ifelse(anyNA(fit.SMLE.cont[[3]], fit.SMLE.cont[[4]]), 0, 1)
        NAsd.LSMLE[l, i]       <- ifelse(anyNA(fit.LSMLE[[3]], fit.LSMLE[[4]]), 0, 1)
        NAsd.LSMLE.cont[l, i]  <- ifelse(anyNA(fit.LSMLE.cont[[3]], fit.LSMLE.cont[[4]]), 0, 1)
        NAsd.MDPDE[l, i]       <- ifelse(anyNA(fit.MDPDE[[3]], fit.MDPDE[[4]]), 0, 1)
        NAsd.MDPDE.cont[l, i]  <- ifelse(anyNA(fit.MDPDE.cont[[3]], fit.MDPDE.cont[[4]]), 0, 1) 
        NAsd.LMDPDE[l, i]      <- ifelse(anyNA(fit.LMDPDE[[3]], fit.LMDPDE[[4]]), 0, 1)
        NAsd.LMDPDE.cont[l, i] <- ifelse(anyNA(fit.LMDPDE.cont[[3]], fit.LMDPDE.cont[[4]]), 0, 1)
        
        
        Total.SMLE[l, i]        <- ifelse(sum(conv.SMLE[l, i] + NAsd.SMLE[l, i]) == 2, 1, 0) # 1 funcinou; 0 erro
        Total.SMLE.cont[l, i]   <- ifelse(sum(conv.SMLE.cont[l, i] + NAsd.SMLE.cont[l, i]) == 2, 1, 0)
        Total.LSMLE[l, i]       <- ifelse(sum(conv.LSMLE[l, i] + NAsd.LSMLE[l, i]) == 2, 1, 0)
        Total.LSMLE.cont[l, i]  <- ifelse(sum(conv.LSMLE.cont[l, i] + NAsd.LSMLE.cont[l, i]) == 2, 1, 0)
        Total.MDPDE[l, i]       <- ifelse(sum(conv.MDPDE[l, i] + NAsd.MDPDE[l, i]) == 2, 1, 0)
        Total.MDPDE.cont[l, i]  <- ifelse(sum(conv.MDPDE.cont[l, i] + NAsd.MDPDE.cont[l, i]) == 2, 1, 0)
        Total.LMDPDE[l, i]      <- ifelse(sum(conv.LMDPDE[l, i] + NAsd.LMDPDE[l, i]) == 2, 1, 0)
        Total.LMDPDE.cont[l, i] <- ifelse(sum(conv.LMDPDE.cont[l, i] + NAsd.LMDPDE.cont[l, i]) == 2, 1, 0)
        
      }
      
      i = i + 1
      
      print(k)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  
  # Salvando resultados
  
  convergencia <- rbind(
    apply(conv.SMLE, 1, mean),
    apply(conv.LSMLE, 1, mean), 
    apply(conv.MDPDE, 1, mean) ,
    apply(conv.LMDPDE, 1, mean),
    apply(conv.SMLE.cont, 1, mean),  
    apply(conv.LSMLE.cont, 1, mean), 
    apply(conv.MDPDE.cont, 1, mean),
    apply(conv.LMDPDE.cont, 1, mean))
  colnames(convergencia) <- alpha
  rownames(convergencia) <- c("SMLE","LSMLE", "MDPDE", "LMDPDE", "SMLE.cont", "LSMLE.cont", "MDPDE.cont", "LMDPDE.cont")
  
  NAsd <- rbind(
    apply(NAsd.SMLE, 1, mean),
    apply(NAsd.LSMLE, 1, mean), 
    apply(NAsd.MDPDE, 1, mean) ,
    apply(NAsd.LMDPDE, 1, mean),
    apply(NAsd.SMLE.cont, 1, mean),  
    apply(NAsd.LSMLE.cont, 1, mean), 
    apply(NAsd.MDPDE.cont, 1, mean),
    apply(NAsd.LMDPDE.cont, 1, mean))
  colnames(NAsd) <- alpha
  rownames(NAsd) <- c("SMLE","LSMLE", "MDPDE", "LMDPDE", "SMLE.cont", "LSMLE.cont", "MDPDE.cont", "LMDPDE.cont")
  
  Total <- rbind(
    apply(Total.SMLE, 1, mean),
    apply(Total.LSMLE, 1, mean), 
    apply(Total.MDPDE, 1, mean) ,
    apply(Total.LMDPDE, 1, mean),
    apply(Total.SMLE.cont, 1, mean),  
    apply(Total.LSMLE.cont, 1, mean), 
    apply(Total.MDPDE.cont, 1, mean),
    apply(Total.LMDPDE.cont, 1, mean))
  colnames(Total) <- alpha
  rownames(Total) <- c("SMLE","LSMLE", "MDPDE", "LMDPDE", "SMLE.cont", "LSMLE.cont", "MDPDE.cont", "LMDPDE.cont")
  
  
  result <- list(
    convergencia = convergencia,
    NAsd = NAsd,
    Total = Total,
    conv.SMLE   = conv.SMLE,   
    conv.LSMLE  = conv.LSMLE,
    conv.MDPDE  = conv.MDPDE,  
    conv.LMDPDE = conv.LMDPDE, 
    conv.SMLE.cont   = conv.SMLE.cont,   
    conv.LSMLE.cont  = conv.LSMLE.cont,  
    conv.MDPDE.cont  = conv.MDPDE.cont,  
    conv.LMDPDE.cont = conv.LMDPDE.cont, 
    NAsd.SMLE   = NAsd.SMLE,  
    NAsd.LSMLE  = NAsd.LSMLE, 
    NAsd.MDPDE  = NAsd.MDPDE, 
    NAsd.LMDPDE = NAsd.LMDPDE,
    NAsd.SMLE.cont   = NAsd.SMLE.cont,  
    NAsd.LSMLE.cont  = NAsd.LSMLE.cont, 
    NAsd.MDPDE.cont  = NAsd.MDPDE.cont, 
    NAsd.LMDPDE.cont = NAsd.LMDPDE.cont,
    Total.SMLE   = Total.SMLE,   
    Total.LSMLE  = Total.LSMLE,  
    Total.MDPDE  = Total.MDPDE,  
    Total.LMDPDE = Total.LMDPDE, 
    Total.SMLE.cont   = Total.SMLE.cont,  
    Total.LSMLE.cont  = Total.LSMLE.cont, 
    Total.MDPDE.cont  = Total.MDPDE.cont, 
    Total.LMDPDE.cont = Total.LMDPDE.cont
  )
  return(result)
}

cenB <- sim.convergencia(1000, 40)

save(cenB, file = "cenB_falhas.RData")

alpha <- seq(0, 1, 0.05)
plot(alpha, 1-cenB$Total[1,], type="l", ylim=c(0,1), ylab="failure rate", las=1, lty=3, xlab=expression(alpha)) # SMLE
lines(alpha, 1-cenB$Total[2,], lty=1) #LSMLE
lines(alpha, 1-cenB$Total[3,], lty=4) #MDPDE
lines(alpha, 1-cenB$Total[4,], lty=2) #LMDPDE
legend("topleft", c("MDPDE", "SMLE", "LMDPDE", "LSMLE"), lty = c(4, 3, 2, 1), bty = "n", cex = 0.9)

plot(alpha, 1-cenB$Total[5,], type="l", ylim=c(0,1), ylab="failure rate", las=1, lty=3, xlab=expression(alpha)) # SMLE
lines(alpha, 1-cenB$Total[6,], lty=1) #LSMLE
lines(alpha, 1-cenB$Total[7,], lty=4) #MDPDE
lines(alpha, 1-cenB$Total[8,], lty=2) #LMDPDE
legend("topleft", c("MDPDE", "SMLE", "LMDPDE", "LSMLE"), lty = c(4, 3, 2, 1), bty = "n", cex = 0.9)

# Scenario C

sim.convergencia = function(M, n){
  
  set.seed(2021 + 04*01)
  
  X = matrix( c( rep(1, 40), runif(40)), ncol = 2, byrow = FALSE )
  
  if(n == 80){
    X <- rbind(X, X)
  }
  par.mu = c(-3, 7.5)
  pred.mu = X%*%par.mu
  mu = exp(pred.mu)/(1+exp(pred.mu)) 
  
  Z = X
  par.phi = c(1, 2)
  pred.phi = Z%*%par.phi
  phi = exp(pred.phi) 
  
  
  alpha <- seq(0, 1, 0.05)
  
  conv.SMLE   <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  conv.LSMLE  <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  conv.MDPDE  <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  conv.LMDPDE <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  
  conv.SMLE.cont   <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  conv.LSMLE.cont  <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  conv.MDPDE.cont  <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  conv.LMDPDE.cont <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  
  
  NAsd.SMLE   <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  NAsd.LSMLE  <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  NAsd.MDPDE  <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  NAsd.LMDPDE <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  
  NAsd.SMLE.cont   <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  NAsd.LSMLE.cont  <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  NAsd.MDPDE.cont  <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  NAsd.LMDPDE.cont <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  
  Total.SMLE   <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  Total.LSMLE  <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  Total.MDPDE  <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  Total.LMDPDE <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  
  Total.SMLE.cont   <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  Total.LSMLE.cont  <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  Total.MDPDE.cont  <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  Total.LMDPDE.cont <- matrix(c(rep(0, length(alpha)*M)), nrow = length(alpha))
  
  i = 1
  k = 0
  
  while(i <= M){
    tryCatch({
      k = k + 1
      set.seed(3*k+2021)
      y  = rbeta(n, mu*phi, (1-mu)*phi)
      
      #contaminando os dados
      y.cont <- y
      m <- round(n*0.05)
      amos1 <- order(mu)[(n-m+1):n]
      mu.linha <- exp(par.mu[1])/(1+exp(par.mu[1]))
      phi.linha <- exp(par.phi[1] + par.phi[2])
      y.cont[amos1] <- rbeta(m, mu.linha*phi.linha, (1 - mu.linha)*phi.linha)
      
      
      for(l in 1:length(alpha)){
        # SMLE
        fit.SMLE <- SMLE_BETA(y, X, Z, qoptimal = FALSE, q0 = 1 - alpha[l])
        fit.SMLE.cont <- SMLE_BETA(y.cont, X, Z, qoptimal = FALSE, q0 = 1 - alpha[l])
        # LSMLE
        fit.LSMLE <- LSMLE_BETA(y, X, Z, qoptimal = FALSE, q0 = 1 - alpha[l])
        fit.LSMLE.cont <- LSMLE_BETA(y.cont, X, Z, qoptimal = FALSE, q0 = 1 - alpha[l])
        # MDPDE
        fit.MDPDE <- MDPDE_BETA(y, X, Z, qoptimal = FALSE, q0 = 1 - alpha[l])
        fit.MDPDE.cont <- MDPDE_BETA(y.cont, X, Z, qoptimal = FALSE, q0 = 1 - alpha[l])
        # LMDPDE
        fit.LMDPDE <- LMDPDE_BETA(y, X, Z, qoptimal = FALSE, q0 = 1 - alpha[l])
        fit.LMDPDE.cont <- LMDPDE_BETA(y.cont, X, Z, qoptimal = FALSE, q0 = 1 - alpha[l])
        
        
        conv.SMLE[l, i]        <- ifelse(is.null(fit.SMLE$Warning), 1, 0) # 1 conv; 0 n?o convergiu
        conv.SMLE.cont[l, i]   <- ifelse(is.null(fit.SMLE.cont$Warning), 1, 0) 
        conv.LSMLE[l, i]       <- ifelse(is.null(fit.LSMLE$Warning), 1, 0) 
        conv.LSMLE.cont[l, i]  <- ifelse(is.null(fit.LSMLE.cont$Warning), 1, 0) 
        conv.MDPDE[l, i]       <- ifelse(is.null(fit.MDPDE$Warning), 1, 0) 
        conv.MDPDE.cont[l, i]  <- ifelse(is.null(fit.MDPDE.cont$Warning), 1, 0) 
        conv.LMDPDE[l, i]      <- ifelse(is.null(fit.LMDPDE$Warning), 1, 0) 
        conv.LMDPDE.cont[l, i] <- ifelse(is.null(fit.LMDPDE.cont$Warning), 1, 0) 
        
        NAsd.SMLE[l, i]        <- ifelse(anyNA(fit.SMLE[[3]], fit.SMLE[[4]]), 0, 1) # 1 funcinou; 0 erro
        NAsd.SMLE.cont[l, i]   <- ifelse(anyNA(fit.SMLE.cont[[3]], fit.SMLE.cont[[4]]), 0, 1)
        NAsd.LSMLE[l, i]       <- ifelse(anyNA(fit.LSMLE[[3]], fit.LSMLE[[4]]), 0, 1)
        NAsd.LSMLE.cont[l, i]  <- ifelse(anyNA(fit.LSMLE.cont[[3]], fit.LSMLE.cont[[4]]), 0, 1)
        NAsd.MDPDE[l, i]       <- ifelse(anyNA(fit.MDPDE[[3]], fit.MDPDE[[4]]), 0, 1)
        NAsd.MDPDE.cont[l, i]  <- ifelse(anyNA(fit.MDPDE.cont[[3]], fit.MDPDE.cont[[4]]), 0, 1) 
        NAsd.LMDPDE[l, i]      <- ifelse(anyNA(fit.LMDPDE[[3]], fit.LMDPDE[[4]]), 0, 1)
        NAsd.LMDPDE.cont[l, i] <- ifelse(anyNA(fit.LMDPDE.cont[[3]], fit.LMDPDE.cont[[4]]), 0, 1)
        
        
        Total.SMLE[l, i]        <- ifelse(sum(conv.SMLE[l, i] + NAsd.SMLE[l, i]) == 2, 1, 0) # 1 funcinou; 0 erro
        Total.SMLE.cont[l, i]   <- ifelse(sum(conv.SMLE.cont[l, i] + NAsd.SMLE.cont[l, i]) == 2, 1, 0)
        Total.LSMLE[l, i]       <- ifelse(sum(conv.LSMLE[l, i] + NAsd.LSMLE[l, i]) == 2, 1, 0)
        Total.LSMLE.cont[l, i]  <- ifelse(sum(conv.LSMLE.cont[l, i] + NAsd.LSMLE.cont[l, i]) == 2, 1, 0)
        Total.MDPDE[l, i]       <- ifelse(sum(conv.MDPDE[l, i] + NAsd.MDPDE[l, i]) == 2, 1, 0)
        Total.MDPDE.cont[l, i]  <- ifelse(sum(conv.MDPDE.cont[l, i] + NAsd.MDPDE.cont[l, i]) == 2, 1, 0)
        Total.LMDPDE[l, i]      <- ifelse(sum(conv.LMDPDE[l, i] + NAsd.LMDPDE[l, i]) == 2, 1, 0)
        Total.LMDPDE.cont[l, i] <- ifelse(sum(conv.LMDPDE.cont[l, i] + NAsd.LMDPDE.cont[l, i]) == 2, 1, 0)
        
      }
      
      i = i + 1
      
      print(k)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  
  convergencia <- rbind(
    apply(conv.SMLE, 1, mean),
    apply(conv.LSMLE, 1, mean), 
    apply(conv.MDPDE, 1, mean) ,
    apply(conv.LMDPDE, 1, mean),
    apply(conv.SMLE.cont, 1, mean),  
    apply(conv.LSMLE.cont, 1, mean), 
    apply(conv.MDPDE.cont, 1, mean),
    apply(conv.LMDPDE.cont, 1, mean))
  colnames(convergencia) <- alpha
  rownames(convergencia) <- c("SMLE","LSMLE", "MDPDE", "LMDPDE", "SMLE.cont", "LSMLE.cont", "MDPDE.cont", "LMDPDE.cont")
  
  NAsd <- rbind(
    apply(NAsd.SMLE, 1, mean),
    apply(NAsd.LSMLE, 1, mean), 
    apply(NAsd.MDPDE, 1, mean) ,
    apply(NAsd.LMDPDE, 1, mean),
    apply(NAsd.SMLE.cont, 1, mean),  
    apply(NAsd.LSMLE.cont, 1, mean), 
    apply(NAsd.MDPDE.cont, 1, mean),
    apply(NAsd.LMDPDE.cont, 1, mean))
  colnames(NAsd) <- alpha
  rownames(NAsd) <- c("SMLE","LSMLE", "MDPDE", "LMDPDE", "SMLE.cont", "LSMLE.cont", "MDPDE.cont", "LMDPDE.cont")
  
  Total <- rbind(
    apply(Total.SMLE, 1, mean),
    apply(Total.LSMLE, 1, mean), 
    apply(Total.MDPDE, 1, mean) ,
    apply(Total.LMDPDE, 1, mean),
    apply(Total.SMLE.cont, 1, mean),  
    apply(Total.LSMLE.cont, 1, mean), 
    apply(Total.MDPDE.cont, 1, mean),
    apply(Total.LMDPDE.cont, 1, mean))
  colnames(Total) <- alpha
  rownames(Total) <- c("SMLE","LSMLE", "MDPDE", "LMDPDE", "SMLE.cont", "LSMLE.cont", "MDPDE.cont", "LMDPDE.cont")
  
  
  result <- list(
    convergencia = convergencia,
    NAsd = NAsd,
    Total = Total,
    conv.SMLE   = conv.SMLE,   
    conv.LSMLE  = conv.LSMLE,
    conv.MDPDE  = conv.MDPDE,  
    conv.LMDPDE = conv.LMDPDE, 
    conv.SMLE.cont   = conv.SMLE.cont,   
    conv.LSMLE.cont  = conv.LSMLE.cont,  
    conv.MDPDE.cont  = conv.MDPDE.cont,  
    conv.LMDPDE.cont = conv.LMDPDE.cont, 
    NAsd.SMLE   = NAsd.SMLE,  
    NAsd.LSMLE  = NAsd.LSMLE, 
    NAsd.MDPDE  = NAsd.MDPDE, 
    NAsd.LMDPDE = NAsd.LMDPDE,
    NAsd.SMLE.cont   = NAsd.SMLE.cont,  
    NAsd.LSMLE.cont  = NAsd.LSMLE.cont, 
    NAsd.MDPDE.cont  = NAsd.MDPDE.cont, 
    NAsd.LMDPDE.cont = NAsd.LMDPDE.cont,
    Total.SMLE   = Total.SMLE,   
    Total.LSMLE  = Total.LSMLE,  
    Total.MDPDE  = Total.MDPDE,  
    Total.LMDPDE = Total.LMDPDE, 
    Total.SMLE.cont   = Total.SMLE.cont,  
    Total.LSMLE.cont  = Total.LSMLE.cont, 
    Total.MDPDE.cont  = Total.MDPDE.cont, 
    Total.LMDPDE.cont = Total.LMDPDE.cont
  )
  return(result)
}

cenC <- sim.convergencia(1000, 40)

load("cenC_falhas.RData")

pdf(file = "failureCenC_nonCont.pdf",   width = 5,  height = 5) 

alpha <- seq(0, 0.5, 0.05)
length(alpha)
plot(alpha, 1-cenC$Total[1,1:length(alpha)], type="l", ylim=c(0,1), ylab="failure rate", las=1, lty=1, xlab=expression(alpha)) # SMLE
abline(h=0, lty=2, col="gray")
abline(h=1, lty=2, col="gray")
points(alpha, 1-cenC$Total[1,1:length(alpha)], pch = 17)
lines(alpha,  1-cenC$Total[2,1:length(alpha)], lty=1) #LSMLE
points(alpha, 1-cenC$Total[2,1:length(alpha)], pch = 16)
lines(alpha,  1-cenC$Total[3,1:length(alpha)], lty=1) #MDPDE
points(alpha, 1-cenC$Total[3,1:length(alpha)], pch=15)
lines(alpha,  1-cenC$Total[4,1:length(alpha)], lty=1) #LMDPDE
points(alpha, 1-cenC$Total[4,1:length(alpha)], pch = 8)
#legend("topleft", c("MDPDE", "SMLE", "LMDPDE", "LSMLE"), lty = c(4, 3, 2, 1), bty = "n", cex = 0.9)
dev.off()


pdf(file = "failureCenC_Cont.pdf",   width = 5,  height = 5) 

alpha <- seq(0, 0.5, 0.05)
length(alpha)
plot(alpha, 1-cenC$Total[5,1:length(alpha)], type="l", ylim=c(0,1), ylab="failure rate", las=1, lty=1, xlab=expression(alpha)) # SMLE
abline(h=0, lty=2, col="gray")
abline(h=1, lty=2, col="gray")
points(alpha, 1-cenC$Total[5,1:length(alpha)], pch = 17)
lines(alpha,  1-cenC$Total[6,1:length(alpha)], lty=1) #LSMLE
points(alpha, 1-cenC$Total[6,1:length(alpha)], pch = 16)
lines(alpha,  1-cenC$Total[7,1:length(alpha)], lty=1) #MDPDE
points(alpha, 1-cenC$Total[7,1:length(alpha)], pch=15)
lines(alpha,  1-cenC$Total[8,1:length(alpha)], lty=1) #LMDPDE
points(alpha, 1-cenC$Total[8,1:length(alpha)], pch = 8)
#legend("topleft", c("MDPDE", "SMLE", "LMDPDE", "LSMLE"), lty = c(4, 3, 2, 1), bty = "n", cex = 0.9)
dev.off()


plot(alpha, 1-cenB$Total[5,], type="l", ylim=c(0,1), ylab="failure rate", las=1, lty=3, xlab=expression(alpha)) # SMLE
abline(h=0, lty=2, col="gray")
abline(h=1, lty=2, col="gray")
lines(alpha, 1-cenB$Total[6,], lty=1) #LSMLE
lines(alpha, 1-cenB$Total[7,], lty=4) #MDPDE
lines(alpha, 1-cenB$Total[8,], lty=2) #LMDPDE
legend("topleft", c("MDPDE", "SMLE", "LMDPDE", "LSMLE"), lty = c(4, 3, 2, 1), bty = "n", cex = 0.9)






