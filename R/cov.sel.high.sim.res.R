cov.sel.high.sim.res<-function(object){

N<-object[[1]]$N      
Setting<-object[[1]]$Setting  
rep<-object[[1]]$rep 
Models<-object[[1]]$Models
type<-object[[1]]$type 
varnames<-object[[1]]$varnames
if(is.null(object[[1]]$est_psm)){betahat<-FALSE}else{betahat<-TRUE}
if(betahat==TRUE){resmat<-matrix(NA,ncol=92,nrow=rep)}else{resmat<-matrix(NA,ncol=20,nrow=rep)}

covarcol<-1:(length(varnames)-2)  
ycol<-length(varnames)-1
Tcol<-length(varnames)

##True ATE   
if(Models=="Binary"){
  if(Setting==1){beta<-0.1235}
  if(Setting==2){beta<-0.0826}
}else{
  beta<-2
}
if(Setting==1){ 
  uc1<-varnames[c(1,2,7)]
  uc2<-varnames[c(1,2,8)]
  XTS<-varnames[c(1,2,3,4,7)]
  Q1S<-Q0S<-QS<-varnames[c(1,2,7)]
  X1S<-X0S<-XYS<-varnames[c(1,2,5,6,8)]
  Z1S<-Z0S<-ZS<-varnames[c(1,2,8)]
  XTYS<-varnames[1:8]
}

if(Setting==2){
  uc1<-varnames[c(1,2,4,7)]
  uc2<-varnames[c(1,2,4,8)]
  XTS<-varnames[c(1,2,3,4,7)]
  Q1S<-Q0S<-QS<-varnames[c(1,2,4,7)]
  X1S<-X0S<-XYS<-varnames[c(1,2,5,6,8)]
  Z1S<-Z0S<-ZS<-varnames[c(1,2,8)]
  XTYS<-varnames[1:8]
  collider<-"x9"
}
  


for(i in 1:rep){
     XT<-object[[i]]$X.T
     Q<-object[[i]]$Q
     XY<-object[[i]]$X.Y
     Z<-object[[i]]$Z
     XTY<-object[[i]]$X.TY
     cards<-numeric(5)
     for(j in 1:5){
                  cards[j]<-object[[i]]$cardinalities[[j]]
     }
     
 
    #Common cause criterion (Waernbaum, de Luna and Richardson)
    ##Algorithm 1
    ##Subset X.T
    if(length(XT)==0){
          XTuc<-SinXT<-XTeqS<-0
          }else{
                XTuc1<-length(which(match(uc1,XT)!="NA"))
                XTuc2<-length(which(match(uc2,XT)!="NA"))
                if(Setting==1){
                        XTuc<-ifelse(XTuc1==length(uc1) || XTuc2==length(uc2),1,0)
                        }else{ 
                              XTuc0<-ifelse(XTuc1==length(uc1) || XTuc2==length(uc2),1,0)
                              XTuc<-ifelse((XTuc0==1) && (length(which(match(collider,XT)!="NA")))==0,1,0)
                              }
                SinXT<-ifelse(length(which(match(XTS,XT)!="NA"))==length(XTS),1,0)
                XTeqS<-ifelse(SinXT==1 && length(XT)==length(XTS),1,0)
                }
    ###Subset Q
    if(length(Q)==0){
          Quc<-SinQ<-QeqS<-0
          }else{
                Quc1<-length(which(match(uc1,Q)!="NA"))
                Quc2<-length(which(match(uc2,Q)!="NA"))
                if(Setting==1){
                        Quc<-ifelse(Quc1==length(uc1) || Quc2==length(uc2),1,0)
                        }else{ 
                              Quc0<-ifelse(Quc1==length(uc1) || Quc2==length(uc2),1,0)
                              Quc<-ifelse((Quc0==1) && (length(which(match(collider,Q)!="NA")))==0,1,0)
                              }
              SinQ<-ifelse(length(which(match(QS,Q)!="NA"))==length(QS),1,0)
              QeqS<-ifelse(SinQ==1 && length(Q)==length(QS),1,0)
              }


    ##Algortithm 2
    ##Subset X.Y
    if(length(XY)==0){
          XYuc<-SinXY<-XYeqS<-0
          }else{
                XYuc1<-length(which(match(uc1,XY)!="NA"))
                XYuc2<-length(which(match(uc2,XY)!="NA"))
                if(Setting==1){
                        XYuc<-ifelse(XYuc1==length(uc1) || XYuc2==length(uc2),1,0)
                        }else{
                              XYuc0<-ifelse(XYuc1==length(uc1) || XYuc2==length(uc2),1,0)
                              XYuc<-ifelse((XYuc0==1) && (length(which(match(collider,XY)!="NA")))==0,1,0)
                              }
                SinXY<-ifelse(length(which(match(XYS,XY)!="NA"))==length(XYS),1,0)
                XYeqS<-ifelse(SinXY==1 && length(XY)==length(XYS),1,0)
                }

    ##Subset Z
    if(length(Z)==0){
          Zuc<-SinZ<-ZeqS<-0
          }else{
                Zuc1<-length(which(match(uc1,Z)!="NA"))
                Zuc2<-length(which(match(uc2,Z)!="NA"))
                if(Setting==1){
                        Zuc<-ifelse(Zuc1==length(uc1) || Zuc2==length(uc2),1,0)
                        }else{ 
                              Zuc0<-ifelse(Zuc1==length(uc1) || Zuc2==length(uc2),1,0)
                              Zuc<-ifelse((Zuc0==1) && (length(which(match(collider,Z)!="NA")))==0,1,0)
                              }
                SinZ<-ifelse(length(which(match(ZS,Z)!="NA"))==length(ZS),1,0)
                ZeqS<-ifelse(SinZ==1 && length(Z)==length(ZS),1,0)
              }
                  
    #Disjunctive cause criterion (VanderWeele and Shpitser)
    ##Subset X.D
    if(length(XTY)==0){
          XTYuc<-SinXTY<-XTYeqS<-0
          }else{             
                XTYuc1<-length(which(match(uc1,XTY)!="NA"))
                XTYuc2<-length(which(match(uc2,XTY)!="NA"))
                if(Setting==1){
                        XTYuc<-ifelse(XTYuc1==length(uc1) || XTYuc2==length(uc2),1,0)
                        }else{ 
                              XTYuc0<-ifelse(XTYuc1==length(uc1) || XTYuc2==length(uc2),1,0)
                              XTYuc<-ifelse((XTYuc0==1) && (length(which(match(collider,XTY)!="NA")))==0,1,0)
                              }
                        SinXTY<-ifelse(length(which(match(XTYS,XTY)!="NA"))==length(XTYS),1,0)
                        XTYeqS<-ifelse(SinXTY==1 && length(XTY)==length(XTYS),1,0)
                }
 
                  
     if(betahat==TRUE){
     betahatest_psm<-betahatest_tmle<-betahatse_psm<-betahatse_tmle<-numeric(6)
     for(j in 1:6){
       betahatest_psm[j]<-object[[i]]$est_psm[[j]]
       betahatse_psm[j]<-object[[i]]$se_psm[[j]]
       betahatest_tmle[j]<-object[[i]]$est_tmle[[j]]
       betahatse_tmle[j]<-object[[i]]$se_tmle[[j]]
     }
     
     
    ##Confidence intervals for ATE estimates                
    ciL_psm<-betahatest_psm-1.96*betahatse_psm
    ciU_psm<-betahatest_psm+1.96*betahatse_psm
    ciwidth_psm<-ciU_psm-ciL_psm
    ciL_tmle<-betahatest_tmle-1.96*betahatse_tmle
    ciU_tmle<-betahatest_tmle+1.96*betahatse_tmle
    ciwidth_tmle<-ciU_tmle-ciL_tmle
    ## Coverage indicator
    cimat1<-matrix(c(ciL_psm,ciU_psm),ncol=2)
    cimat2<-matrix(c(ciL_tmle,ciU_tmle),ncol=2)
    cifunc<-function(cimat){
      cicov<-ifelse(cimat[1]<beta && cimat[2]>beta,1,0)
      }
    betahat_cicov_psm<-apply(cimat1,1,cifunc) 
    betahat_cicov_tmle<-apply(cimat2,1,cifunc) 
    
    #Matrix with results for every replication
    resmat[i,]<-c(XTuc,Quc,XYuc,Zuc,XTYuc,
                  SinXT,SinQ,SinXY,SinZ,SinXTY,
                  XTeqS,QeqS,XYeqS,ZeqS,XTYeqS,
                  cards, betahatest_psm, betahatse_psm, 
                  betahat_cicov_psm, ciL_psm, ciU_psm,
                  betahatest_tmle, betahatse_tmle,
                  betahat_cicov_tmle, ciL_tmle, ciU_tmle,ciwidth_psm,ciwidth_tmle)
    
   
     }else{
       
       resmat[i,]<-c(XTuc,Quc,XYuc,Zuc,XTYuc,
                     SinXT,SinQ,SinXY,SinZ,SinXTY,
                     XTeqS,QeqS,XYeqS,ZeqS,XTYeqS,
                     cards)
       
     }
      
    }
#Summarizes the simulation results            

##Proportion of replications where the three different subset condtions are fullfilled
ss_uc<-matrix(colMeans(resmat)[1:15],nrow=5)
##Median cardinalities
med_cards<-apply(resmat[,16:20],2,median) 

if(betahat==TRUE){
##ATE bias
betahat_bias_psm<-colMeans(resmat[,21:26],na.rm=TRUE)-beta
##ATE SD
betahat_sd_psm<-apply(resmat[,21:26],2,sd,na.rm=TRUE) 
##ATE MSE
betahat_mse_psm<-betahat_bias_psm^2+apply(resmat[,21:26],2,var,na.rm=TRUE)
##Mean CI coverage
betahat_meancoverage_psm<-colMeans(resmat,na.rm=TRUE)[33:38]
##Mean lower CI
mean_ciL_psm<-colMeans(resmat[,39:44],na.rm=TRUE)
##Mean upper CI
mean_ciU_psm<-colMeans(resmat[,45:50],na.rm=TRUE)
##Mean CI width
mean_ciwidth_psm<-colMeans(resmat[,81:86],na.rm=TRUE)


##ATE bias
betahat_bias_tmle<-colMeans(resmat[,51:56],na.rm=TRUE)-beta
##ATE SD
betahat_sd_tmle<-apply(resmat[,51:56],2,sd,na.rm=TRUE) 
##ATE MSE
betahat_mse_tmle<-betahat_bias_tmle^2+apply(resmat[,51:56],2,var,na.rm=TRUE)
##Mean CI coverage
betahat_meancoverage_tmle<-colMeans(resmat,na.rm=TRUE)[63:68]
##Mean lower CI
mean_ciL_tmle<-colMeans(resmat[,69:74],na.rm=TRUE)
##Mean upper CI
mean_ciU_tmle<-colMeans(resmat[,75:80],na.rm=TRUE)
##Mean CI width
mean_ciwidth_tmle<-colMeans(resmat[,87:92],na.rm=TRUE)



##List containing simulation results                              
summary_resmat<-list(Subset_selection=ss_uc,Median_cardinality=med_cards,
                     Betahat_bias_psm=betahat_bias_psm,Betahat_sd_psm=betahat_sd_psm,
                     Betahat_mse_psm=betahat_mse_psm,Betahat_CI_coverage_psm=betahat_meancoverage_psm,Betahat_CI_width_psm=mean_ciwidth_psm,
                     Betahat_mean_lower_CI_psm=mean_ciL_psm,Betahat_mean_upper_CI_psm=mean_ciU_psm,
                     Betahat_bias_tmle=betahat_bias_tmle,Betahat_sd_tmle=betahat_sd_tmle,
                     Betahat_mse_tmle=betahat_mse_tmle,Betahat_CI_coverage_tmle=betahat_meancoverage_tmle,Betahat_CI_width_tmle=mean_ciwidth_tmle,
                     Betahat_mean_lower_CI_tmle=mean_ciL_tmle,Betahat_mean_upper_CI_tmle=mean_ciU_tmle)
##Table for producing LaTeX table with simulation results                              
xtab11<-data.frame(matrix(c( N, rep(" ",5), " ", type, rep(" ", 4)), ncol=2))
xtab21<-data.frame(matrix(c("$X$", "$\\hat{X}_{\\rightarrow T}$", "$\\hat{Q}_{\\rightarrow T}$" ,"$ \\hat{X}_{\\rightarrow Y}$" , "$\\hat{Z}_{\\rightarrow Y}$",
                            "$\\hat{X}_{\\rightarrow T, Y}$" ), ncol=1))
if(Setting==1){
  xtab31<-data.frame(rbind(c(1,1,1),summary_resmat$Subset_selection)*100)
}else{
  xtab31<-data.frame(rbind(c(0,1,0),summary_resmat$Subset_selection)*100)
  
}
xtab41<-data.frame(matrix(c(100,summary_resmat$Median_cardinality,summary_resmat$Betahat_bias_psm,
                            summary_resmat$Betahat_sd_psm,summary_resmat$Betahat_mse_psm,
                            summary_resmat$Betahat_bias_tmle,summary_resmat$Betahat_sd_tmle,summary_resmat$Betahat_mse_tmle),ncol=7))
xtab1<-cbind(xtab11,xtab21, xtab31,xtab41)
digits1<-c(rep(0,4),rep(1,3),0,rep(3,6))
align1<-c(rep("l",4),rep("r",10))
colnames(xtab1)<-c("n","Method","$\\hat{S}$" ," $Y_t \\perp\\!\\!\\!\\perp T \\mid \\hat{S}$",
                   "$S \\subseteq \\hat{S}$","$S=\\hat{S}$",  "\\#", "Bias", "SD" ,"MSE", "Bias", "SD" ,"MSE")
xtab_res1<-xtable(x=xtab1,digits=digits1,align=align1)

##Table for producing LaTeX table with simulation results                              
xtab12<-data.frame(matrix(c(N, rep(" ",5), " ", type, rep(" ", 4)), ncol=2))
xtab22<-data.frame(matrix(c("$X$", "$\\hat{X}_{\\rightarrow T}$", "$\\hat{Q}_{\\rightarrow T}$" ,"$ \\hat{X}_{\\rightarrow Y}$" , "$\\hat{Z}_{\\rightarrow Y}$",
                            "$\\hat{X}_{\\rightarrow T, Y}$" ), ncol=1))
if(Setting==1){
  xtab32<-data.frame(rbind(c(1,1,1),summary_resmat$Subset_selection)*100)[,1]
}else{
  xtab32<-data.frame(rbind(c(0,1,0),summary_resmat$Subset_selection)*100)[,1]
  
}
xtab42<-data.frame(matrix(c(summary_resmat$Betahat_CI_coverage_psm*100,summary_resmat$Betahat_CI_width_psm, summary_resmat$Betahat_mean_lower_CI_psm,
                            summary_resmat$Betahat_mean_upper_CI_psm,
                            summary_resmat$Betahat_CI_coverage_tmle*100,summary_resmat$Betahat_CI_width_tmle, summary_resmat$Betahat_mean_lower_CI_tmle,
                            summary_resmat$Betahat_mean_upper_CI_tmle),ncol=8))
xtab2<-cbind(xtab12,xtab22, xtab32,xtab42)
digits2<-c(rep(0,4),rep(1,3),rep(3,2),1,1,rep(3,2))
align2<-c(rep("l",4),rep("r",9))
colnames(xtab2)<-c("n","Method","$\\hat{S}$" ," $Y_t \\perp\\!\\!\\!\\perp T \\mid \\hat{S}$", "CP","CIW" ,"CIL", "CIU", "CP" ,"CIW","CIL", "CIU")
xtab_res2<-xtable(x=xtab2,digits=digits2,align=align2)

}else{
  summary_resmat<-list(Subset_selection=ss_uc,Median_cardinality=med_cards)
  ##Table for producing LaTeX table with simulation results                              
  xtab11<-data.frame(matrix(c( N, rep(" ",5), " ", type, rep(" ", 4)), ncol=2))
  xtab21<-data.frame(matrix(c("$X$", "$\\hat{X}_{\\rightarrow T}$", "$\\hat{Q}_{\\rightarrow T}$" ,"$ \\hat{X}_{\\rightarrow Y}$" , "$\\hat{Z}_{\\rightarrow Y}$",
                              "$\\hat{X}_{\\rightarrow T, Y}$" ), ncol=1))
  if(Setting==1){
    xtab31<-data.frame(rbind(c(1,1,1),summary_resmat$Subset_selection)*100)
  }else{
    xtab31<-data.frame(rbind(c(0,1,0),summary_resmat$Subset_selection)*100)
    
  }
  xtab41<-data.frame(matrix(c(100,summary_resmat$Median_cardinality),ncol=1))
  xtab1<-cbind(xtab11,xtab21, xtab31,xtab41)
  digits1<-c(rep(0,4),rep(1,3),0)
  align1<-c(rep("l",4),rep("r",4))
  colnames(xtab1)<-c("n","Method","$\\hat{S}$" ," $Y_t \\perp\\!\\!\\!\\perp T \\mid \\hat{S}$",
                     "$S \\subseteq \\hat{S}$","$S=\\hat{S}$",  "\\#")
  xtab_res1<-xtable(x=xtab1,digits=digits1,align=align1)
  
  
  xtab_res2<-NULL
  
  
  
  
  
}
                               
return(list(resmat=resmat,summary_resmat=summary_resmat,xtable1=xtab_res1,xtable2=xtab_res2))
}
