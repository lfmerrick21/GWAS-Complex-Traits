#A function that will filter a genotype matrix based on maf and missingness
#Calculates the proportion of missing data for every marker in a genotype matrix (mising data is NA)
calc_missrate <- function(gt_mat)
{
  col_func <- function(gt_col)
  {
    missrate <- sum(is.na(gt_col)) / length(gt_col)
    return(missrate)
  }
  
  missrate_vect <- apply(gt_mat, 2, col_func)
  
  return(missrate_vect)
}

# Calculates the minor allele frequency for every marker in a genotype matrix (coded as c(-1,0,1))
calc_maf_apply <- function(gt_mat, encoding = c(-1, 0, 1))
{
  col_func1 <- function(gt_col)
  {
    allele1_ct <- (sum(gt_col == -1, na.rm = T) * 2) + sum(gt_col == 0, na.rm = T)
    allele2_ct <- (sum(gt_col == 1, na.rm = T) * 2) + sum(gt_col == 0, na.rm = T)
    
    maf <- min(c(allele1_ct, allele2_ct)) / (sum(!is.na(gt_col))*2)
  }
  
  col_func2 <- function(gt_col)
  {
    allele1_ct <- (sum(gt_col == 0, na.rm = T) * 2) + sum(gt_col == 1, na.rm = T)
    allele2_ct <- (sum(gt_col == 2, na.rm = T) * 2) + sum(gt_col == 1, na.rm = T)
    
    maf <- min(c(allele1_ct, allele2_ct)) / (sum(!is.na(gt_col))*2)
  }
  
  if (all(encoding == c(-1, 0, 1)))
  {
    maf_vect <- apply(gt_mat, 2, col_func1)
  } else if (all(encoding == c(0, 1, 2)))
  {
    maf_vect <- apply(gt_mat, 2, col_func2)
  } else{
    print('Encoding not recognized, returning NULL')
    maf_vect <- NULL
  }
  
  return(maf_vect)
}

# This is a function that will split data into a list of k-folds
make_CV_sets <- function(list_length, k = 5)
{
  rand_values <- rnorm(list_length)
  k_quantiles <- quantile(rand_values, 0:k/k)
  k_assign <- cut(rand_values, k_quantiles, include.lowest = T, labels = F)
  
  cv_list <- list()
  for (i in 1:k)
  {
    fold_assignment <- k_assign != i
    cv_list[[i]] <- fold_assignment
  }
  return(cv_list)
}
#Sommer_MM_GWAS
sommer_MTMM=function(Y=Gen_Table_JMP_DP18[,c(1,2,4)],SNP_INFO=GBS_qam_adj18$map,model=ansx18,X=myGM18,A=A,MAFX=NULL){
  colnames(SNP_INFO)<-c('SNP','Chr','Pos')
  ## preparing genotype and Kinship data 
  Y_<-Y
  X<-X[which(rownames(X)%in%Y_[,1]),]
  Y<-Y_[which(Y_[,1]%in%rownames(X)),]
  ### enter here filtering step for 2 different kinships for different subsets .
  
  Y1<-(Y[,2])
  Y2<-(Y[,3])
  names(Y1)<-Y[,1]
  names(Y2)<-Y[,1]
  Y1_<-na.omit(Y1)
  Y2_<-na.omit(Y2)
  
  ecot_id1<-as.integer(names(Y1_))
  ecot_id2<-as.integer(names(Y2_))
  
  K<-as.matrix(A[which(rownames(A)%in%names(Y1)),which(colnames(A)%in%names(Y1))])
  K1<-as.matrix(A[which(rownames(A)%in%ecot_id1),which(colnames(A)%in%ecot_id1)])
  K2<-as.matrix(A[which(rownames(A)%in%ecot_id2),which(colnames(A)%in%ecot_id2)])
  
  n<-nrow(Y)
  n1<-length(ecot_id1)
  n2<-length(ecot_id2)
  
  
  # combining the traits
  
  Y_ok<-c(Y1,Y2)
  
  #environment
  
  Env<-c(rep(0,n),rep(1,n))
  
  #standardize kinship 
  
  K_stand<-(n-1)/sum((diag(n)-matrix(1,n,n)/n)*K)*K
  K_inv<<-solve(K_stand)
  #K_stand1<-(n1-1)/sum((diag(n1)-matrix(1,n1,n1)/n1)*K1)*K1
  #K_stand2<-(n2-1)/sum((diag(n2)-matrix(1,n2,n2)/n2)*K2)*K2
  
  count2<-function(x) {length(which(x==2))}
  count1<-function(x) {length(which(x==1))}
  
  AC_2 <- data.frame(colnames(X),apply(X,2,count2))
  colnames(AC_2)<-c('SNP','AC_2')
  AC_1 <- data.frame(colnames(X),apply(X,2,count1))
  colnames(AC_1)<-c('SNP','AC_1')
  
  MAF_1<-data.frame(AC_2,AC_1$AC_1,AC_0=nrow(X)-AC_1$AC_1-AC_2$AC_2)
  
  MAF_2<-data.frame(MAF_1,MAC=apply(MAF_1[,c(2,4)],1,min))
  
  MAF_3<-data.frame(MAF_2,MAF=(MAF_2$MAC/nrow(X)))
  
  MAF_ok<<-merge(SNP_INFO,MAF_3,by='SNP')
  
  rm(AC_1,MAF_1,MAF_2,MAF_3)
  
  #Filter for MAF
  
  
  
  
  if(!is.null(MAFX)){
    MAF<-subset(MAF_ok,MAF==0)[,1]
    X_ok<-X[,!colnames(X) %in% MAF]
  }else{
    X_ok<-X
  }
  
  Xo<-rep(1,2*n)
  Xo1<-rep(1,n1)
  ex1<-as.matrix(Xo1)
  Xo2<-rep(1,n2)
  ex2<-as.matrix(Xo2)
  
  
  varcov<-model$sigma[,,1]
  ve<-model$sigma[,,2]
  
  rho_g<-varcov[1,2]/sqrt(varcov[1,1]*varcov[2,2])
  rho_e<-ve[1,2]/sqrt(ve[1,1]*ve[2,2])
  rho_p<-(varcov[1,2]+ve[1,2])/sqrt((varcov[1,1]+ve[1,1])*(varcov[2,2]+ve[2,2]))
  pearson<-cor(Y1,Y2)
  heriti1<-varcov[1,1]/(varcov[1,1]+ve[1,1])
  heriti2<-varcov[2,2]/(varcov[2,2]+ve[2,2])
  correlation<-list(pearson=pearson,gen_cor=rho_g,env_cor=rho_e,phen_cor=rho_p,h1_joint=heriti1,h2_joint=heriti2,converge=model$convergence)
  
  K_comb<-kronecker(varcov,K_stand)
  rownames(K_comb)<-c(ecot_id1,ecot_id2)
  colnames(K_comb)<-c(ecot_id1,ecot_id2)
  
  I_comb<-kronecker(ve,diag(n))
  rownames(I_comb)<-rownames(K_comb)
  colnames(I_comb)<-colnames(K_comb)
  
  bigK<-K_comb+I_comb
  
  
  M<-solve(chol(bigK))
  
  
  # scaling of the SNPs to interpret GxE results
  X_ok1<-X_ok*sqrt(varcov[1,1])
  X_ok2<-X_ok*sqrt(varcov[2,2])
  
  Y_t<-crossprod(M,Y_ok)
  cof_t<-crossprod(M,cbind(Xo,Env))
  
  RSS_env<-sum(lsfit(cof_t,Y_t,intercept = FALSE)$residuals^2)
  
  nbchunks<-5
  m<-ncol(X_ok)
  
  RSS_full<-list()
  RSS_glob<-list()
  
  for (j in 1:(nbchunks-1)){
    X_<-rbind(X_ok1[,((j-1)*round(m/nbchunks)+1):(j*round(m/nbchunks))],X_ok2[,((j-1)*round(m/nbchunks)+1):(j*round(m/nbchunks))])
    X_t<-array(dim=c(nrow(X_),ncol(X_),2))
    X_t[,,1]<-crossprod(M,X_)
    X_t[,,2]<-crossprod(M,(Env*X_))
    RSS_full[[j]]<-apply(X_t,2,function(x){sum(lsfit(cbind(cof_t,x),Y_t,intercept=FALSE)$residuals^2)})
    RSS_glob[[j]]<-apply(X_t[,,1],2,function(x){sum(lsfit(cbind(cof_t,x),Y_t,intercept=FALSE)$residuals^2)})
    rm(X_,X_t)}
  X_<-rbind(X_ok1[,((j)*round(m/nbchunks)+1):m],X_ok2[,((j)*round(m/nbchunks)+1):m])
  X_t<-array(dim=c(nrow(X_),ncol(X_),2))
  X_t[,,1]<-crossprod(M,X_)
  X_t[,,2]<-crossprod(M,(Env*X_))
  RSS_full[[nbchunks]]<-apply(X_t,2,function(x){sum(lsfit(cbind(cof_t,x),Y_t,intercept=FALSE)$residuals^2)})
  RSS_glob[[nbchunks]]<-apply(X_t[,,1],2,function(x){sum(lsfit(cbind(cof_t,x),Y_t,intercept=FALSE)$residuals^2)})
  rm(X_,X_t,j)
  
  RSS_full<-unlist(RSS_full)
  RSS_glob<-unlist(RSS_glob)
  
  #nb parameters in each models
  #env
  par_env<-ncol(cof_t)
  par_glob<-par_env+1
  par_full<-par_glob+1
  
  ##FTESTS
  #FULL vs NULL
  F_full<-(rep(RSS_env,m)/RSS_full-1)*(2*n-par_full)/(par_full-par_env)
  F_ge<-(RSS_glob/RSS_full-1)*(2*n-par_full)/(par_full-par_glob)
  F_glob<-(rep(RSS_env,m)/RSS_glob-1)*(2*n-par_glob)/(par_glob-par_env)
  
  pval_full<-pf(F_full,par_full-par_env,2*n-par_full,lower.tail=FALSE)
  pval_ge<-pf(F_ge,par_full-par_glob,2*n-par_full,lower.tail=FALSE)
  pval_glob<-pf(F_glob,par_glob-par_env,2*n-par_glob,lower.tail=FALSE)
  
  #outputs
  outfull<<-data.frame('SNP'=colnames(X_ok),pval=pval_full)
  outge<<-data.frame('SNP'=colnames(X_ok),pval=pval_ge)
  outglob<<-data.frame('SNP'=colnames(X_ok),pval=pval_glob)
  
  
  
  sommer_scores=data.frame(t(model$scores))
  py=sommer_scores[,3]
  if(max(py,na.rm=TRUE)>1){
    pvaly<-10^-py
  }else{
    pvaly<-py
  }
  sommer_scores$Y.score=pvaly
  pc=sommer_scores[,4]
  if(max(pc,na.rm=TRUE)>1){
    pvalc<-10^-pc
  }else{
    pvalc<-pc
  }
  sommer_scores$DP_BLUP.score=pvalc
  sommer_scores=cbind(SNP=rownames(sommer_scores),sommer_scores)
  
  data.out<-merge(MAF_ok,cbind(outfull,outge[,2],outglob[,2]),by='SNP')
  data.outt<-left_join(data.out,sommer_scores,by='SNP')
  colnames(data.outt)[9:11]<-c('pval_full','pval_trait_specific','pval_trait_common')
  data.out_<-data.outt[order(data.outt[,3]),]
  out<-data.out_[order(data.out_[,2]),]
  results<-list(phenotype=Y,pvals=out,statistics=correlation,kinship=K)
  return(results)
}

#just need phenotypic vector and GWAS results

Marker_Effects<-function(Pheno=NULL,GWAS=NULL,alpha=0.05,correction="Bonferonni",messages=TRUE,model="BLINK"){
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(dplyr, compiler)
  #The downloaded link at: http://cran.r-project.org/package=scatterplot3d
  #source("http://zzlab.net/GAPIT/emma.txt")
  #source("http://zzlab.net/GAPIT/gapit_functions.txt")
  library(GAPIT3)
  #Input your phenotypic vector
  myY_train <- Pheno
  if(messages==TRUE){print(paste0("Marker effects are being calculated and using the ",correction," correction."))}
  #Correction Method
  #Bonferonni Correction
  #FDR Correction
  if(model=="BLINK"){
    if(correction=="FDR"){
      GWAS1=GWAS$GWAS
      GWAS1$Rs=GWAS1$effect
      GWASFDR=GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure(GWAS1)
      sig_mar <- which(GWASFDR$PWIP$`FDR_Adjusted_P-values` <= alpha)
      sig_markers <-as.numeric(c(rownames(GWASFDR$PWIP)[sig_mar]))
    }else{
      sig_markers <- which(GWAS$GWAS$P.value <= alpha/length(GWAS$GWAS$P.value))
    }}
  
  if(model=="FarmCPU"){
    if(correction=="FDR"){
      GWAS1=GWAS$GWAS
      GWAS1$Rs=GWAS1$effect
      GWASFDR=GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure(GWAS1)
      sig_mar <- which(GWASFDR$PWIP$`FDR_Adjusted_P-values` <= alpha)
      sig_markers <-as.numeric(c(rownames(GWASFDR$PWIP)[sig_mar]))
    }else{
      sig_markers <- which(GWAS$GWAS$P.value <= alpha/length(GWAS$GWAS$P.value))
    }}
  
  if(model=="SUPER"){
    if(correction=="FDR"){
      GWAS1=GWAS$GWAS
      GWAS1=GWAS1[,!names(GWAS$GWAS) %in% c("effect")]
      GWASFDR=GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure(GWAS1)
      sig_mar <- which(GWASFDR$PWIP$`FDR_Adjusted_P-values` <= alpha)
      sig_markers <-as.numeric(c(rownames(GWASFDR$PWIP)[sig_mar]))
    }else{
      sig_markers <- which(GWAS$GWAS$P.value <= alpha/length(GWAS$GWAS$P.value))
    }}
  
  if(model=="MLM"){
    if(correction=="FDR"){
      GWAS1=GWAS$GWAS
      GWAS1=GWAS1[,!names(GWAS$GWAS) %in% c("effect")]
      GWASFDR=GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure(GWAS1)
      sig_mar <- which(GWASFDR$PWIP$`FDR_Adjusted_P-values` <= alpha)
      sig_markers <-as.numeric(c(rownames(GWASFDR$PWIP)[sig_mar]))
    }else{
      sig_markers <- which(GWAS$GWAS$P.value <= alpha/length(GWAS$GWAS$P.value))
    }}
  
  
  
  #Regular P-value
  if(correction=="P-value"){sig_markers <- which(GWAS$GWAS$P.value <= alpha)}
  PC_mat <-GWAS$PCA[,-1]
  MAS_mat<-GWAS$GD[,-1]
  MAS_mat<-data.frame(MAS_mat[,sig_markers])
  if(length(sig_markers)==1){names(MAS_mat)<-GWAS$GWAS[sig_markers,]$SNP[[1]]}
  if(length(sig_markers)==0){print("Your GWAS is terrible and Halle thinks you suck.")}
  if(length(sig_markers)==0){stop("Your GWAS and correction method found zero significant markers.")}
  if(messages==TRUE){print(paste0("There are ",length(sig_markers)," significant marker(s)."))}
  GLM_data <- data.frame(myY_train,PC_mat,MAS_mat)
  names(GLM_data)[1] <- "Y"
  MAS_model <- lm(Y ~ ., data = GLM_data)
  Marker_Results=data.frame()
  cov_length=ncol(data.frame(myY_train))+ncol(PC_mat)
  for(i in cov_length+1:ncol(as.data.frame(MAS_mat)) ){
    #Null Model
    GLM_data_null <- data.frame(myY_train,PC_mat)
    names(GLM_data_null)[1] <- "Y"
    MAS_model_null <- lm(Y ~ ., data = GLM_data_null)
    R2null=summary(MAS_model_null)$r.squared
    
    #Model with 1 marker
    SNP_Name=names(GLM_data)[i]
    GLM_data_w1 <- data.frame(myY_train,PC_mat,GLM_data[,i])
    names(GLM_data_w1)[1] <- "Y"
    names(GLM_data_w1)[cov_length+1] <- SNP_Name
    MAS_model_w1 <- lm(Y ~ ., data = GLM_data_w1)
    R2w1=summary(MAS_model_w1)$r.squared
    R2w1m=R2w1-R2null #r2 for significant marker
    
    #Model without marker
    GLM_data_sm1 <- GLM_data[,-i]
    MAS_model_sm1 <- lm(Y ~ ., data = GLM_data_sm1)
    R2full=summary(MAS_model)$r.squared
    R2wo=summary(MAS_model_sm1)$r.squared
    R2wo1m=summary(MAS_model)$r.squared-summary(MAS_model_sm1)$r.squared #r2 for significant marker
    Full_Effects=MAS_model$coefficients[i]
    SNP_Effects=MAS_model_w1$coefficients[cov_length+1]
    
    results=data.frame(SNP_Name,R2null,R2w1,R2w1m,R2full,R2wo,R2wo1m,Full_Effects,SNP_Effects)
    names(results)=c("SNP","R2.null.model","R2.add.single.marker","R2.of.single.marker","R2.full.model","R2.minus.single.marker","R2.diff","Effects.marker.full.model","Effects.single.marker")
    Marker_Results=rbind(Marker_Results,results)
    
    if(correction=="FDR"){
      GWAS_Table=left_join(GWASFDR$PWIP[sig_mar,!names(GWASFDR$PWIP) %in% c("effect","Rs")],Marker_Results,by="SNP")
      GWAS_Table <- GWAS_Table[order(GWAS_Table$SNP),]
    }else{
      GWAS_Table=left_join(GWAS$GWAS[sig_markers,!names(GWAS$GWAS) %in% c("effect")],Marker_Results,by="SNP")
    }
    
    if(model=="BLINK"){
      if(correction=="FDR"){
        GWAS_Table=left_join(GWASFDR$PWIP[sig_mar,!names(GWASFDR$PWIP) %in% c("effect","Rs")],Marker_Results,by="SNP")
        GWAS_Table <- GWAS_Table[order(GWAS_Table$SNP),]
      }else{
        GWAS_Table=left_join(GWAS$GWAS[sig_markers,!names(GWAS$GWAS) %in% c("effect")],Marker_Results,by="SNP")
      }}
    
    if(model=="FarmCPU"){
      if(correction=="FDR"){
        GWAS_Table=left_join(GWASFDR$PWIP[sig_mar,!names(GWASFDR$PWIP) %in% c("effect","Rs")],Marker_Results,by="SNP")
        GWAS_Table <- GWAS_Table[order(GWAS_Table$SNP),]
      }else{
        GWAS_Table=left_join(GWAS$GWAS[sig_markers,!names(GWAS$GWAS) %in% c("effect")],Marker_Results,by="SNP")
      }}
    
    if(model=="SUPER"){
      if(correction=="FDR"){
        GWAS_Table=left_join(GWASFDR$PWIP[sig_mar,!names(GWASFDR$PWIP) %in% c("effect")],Marker_Results,by="SNP")
        GWAS_Table <- GWAS_Table[order(GWAS_Table$SNP),]
      }else{
        GWAS_Table=left_join(GWAS$GWAS[sig_markers,!names(GWAS$GWAS) %in% c("effect")],Marker_Results,by="SNP")
      }}
    
    if(model=="MLM"){
      if(correction=="FDR"){
        GWAS_Table=left_join(GWASFDR$PWIP[sig_mar,!names(GWASFDR$PWIP) %in% c("effect")],Marker_Results,by="SNP")
        GWAS_Table <- GWAS_Table[order(GWAS_Table$SNP),]
      }else{
        GWAS_Table=left_join(GWAS$GWAS[sig_markers,!names(GWAS$GWAS) %in% c("effect")],Marker_Results,by="SNP")
      }}
    
    
    
    
  }
  if(messages==TRUE){print(paste0("Wow! Your signficant markers and covariates accounted for ",round(R2full*100,2),"% of the variation."))}
  return(GWAS_Table)
}


Marker_Effects_Allele<-function(Pheno=NULL,GWAS=NULL,alpha=0.05,correction="Bonferonni",messages=TRUE,Markers=NULL,model="BLINK"){
  
  #Pheno=GBS_qam_adj22$pheno[,c(22)];GWAS=myGAPIT_MLM_qam_adj_22;alpha=0.05;correction="FDR";messages=TRUE;Markers=EM_GBS_SNPs;model="MLM"
  #Pheno=GBS_qam_adj24$pheno[,c(24)];GWAS=FarmCPU_qam_adj_24;alpha=0.05;correction="Bonferonni";messages=TRUE;Markers=EM_GBS_SNPs
  
  
  
  corr_sig=Marker_Effects(Pheno=Pheno,GWAS=GWAS,alpha=alpha,correction=correction,messages=messages,model=model)
  names(Markers)[1]<-"SNP"
  GWAS_Results=left_join(corr_sig,Markers[,1:3],by="SNP")
  
  if(model=="BLINK"){
    if(correction=="FDR"){
      GWAS_Results=GWAS_Results[,c(1:3,16,17,4:15)]
    }else{
      GWAS_Results=GWAS_Results[,c(1:3,15,16,4:14)]
    }}
  
  if(model=="FarmCPU"){
    if(correction=="FDR"){
      GWAS_Results=GWAS_Results[,c(1:3,16,17,4:15)]
    }else{
      GWAS_Results=GWAS_Results[,c(1:3,15,16,4:14)]
    }}
  
  if(model=="SUPER"){
    if(correction=="FDR"){
      GWAS_Results=GWAS_Results[,c(1:3,17,18,4:16)]
    }else{
      GWAS_Results=GWAS_Results[,c(1:3,17,18,4:16)]
    }}
  
  if(model=="MLM"){
    if(correction=="FDR"){
      GWAS_Results=GWAS_Results[,c(1:3,18,19,4:6,9:17)]
    }else{
      GWAS_Results=GWAS_Results[,c(1:3,17,18,4:6,9:16)]
    }}
  return(GWAS_Results)
}


manhattan_plot <- function(GWAS,cutoff=NULL,QTN_index = NULL, trait = NULL,labels=NULL,model="unknown",Nchr=22){
  #GWAS=myGAPIT_MLM_qam_adj_23$GWAS
  #cutoff=NULL
  #QTN_index=Bon_sig_B1_bl_CV_Col$SNP[2:3]
  #Multi=myGAPIT_MLM_adj_13$GWAS
  #Multi_labels=Bon_sig_M1_bl$SNP
  #Third_Labels=Bon_sig_B1_Col_Len_blup$SNP
  #Trait_one="Diversity Panel"
  #Trait_two="Breeding Lines"
  #model="No Covariates"
  #labels = Bon_sig_M6$SNP
  
  m=nrow(GWAS)
  cutoff.final=(ifelse(
    is.null(cutoff),
    0.05/m,
    cutoff
  ))
  SG_pval <- -log10(cutoff.final)
  
  col.Oceanic=rep(c(  '#EC5f67',    '#FAC863',  '#99C794',    '#6699CC',  '#C594C5'),ceiling(Nchr/5))
  
  don <- GWAS %>% 
    # Compute chromosome size
    group_by(Chromosome) %>% 
    summarise(chr_len=max(Position)) %>% 
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(GWAS, ., by=c("Chromosome"="Chromosome")) %>%
    # Add a cumulative position of each SNP
    arrange(Chromosome, Position) %>%
    mutate( BPcum=Position+tot)
  
  if(!is.null(labels)){
    don=don%>% mutate( is_annotate=ifelse(SNP %in% labels, "yes", "no"))
    don=don%>%mutate( is_highlight=ifelse(P.value<=cutoff.final, "yes", "no"))
  }
  
  #Then we need to prepare the X axis. Indeed we do not want to display the cumulative position of SNP in bp, but just show the chromosome name instead.
  axisdf = don %>% group_by(Chromosome) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  #Cutoff
  
  
  
  manhattan_plot <- ggplot(don, aes(x=BPcum, y=-log10(P.value))) +
    # Show all points
    geom_point( aes(color=as.factor(Chromosome)), alpha=0.8, size=1.3) +
    scale_color_manual(values = col.Oceanic) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$Chromosome, breaks= axisdf$center ) +
    scale_y_continuous(limits = c(0,  12),labels =scales::number_format(accuracy=1)   ) +     # remove space between plot area and x axis
    #scales::number_format(accuracy = 0)   
    #ylim(0,  max(-log10(don$P.value)))+
    # Custom the theme:
    geom_hline(yintercept= SG_pval,color="green")+
    labs(#title = paste("GWAS manhattan plot for trait:", as.character(trait)),
      y = "-log10(p)",
      x = "Chromosome",
      color = "Chromosome",
      tag=model)+
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.tag=element_text(angle=-90,size=10),
      plot.tag.position=c("right")
    )
  if(!is.null(QTN_index)){
    QTN=don%>%filter(SNP %in% QTN_index)
    manhattan_plot=manhattan_plot+geom_vline(xintercept = QTN$BPcum, color = "red") }
  if(!is.null(trait)){manhattan_plot=manhattan_plot+labs(title = paste("GWAS manhattan plot for ", as.character(trait))) }
  if(!is.null(labels)){
    manhattan_plot=manhattan_plot+
      # Add highlighted points
      geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=3) +
      # Add label using ggrepel to avoid overlapping
      geom_text_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=3)
  }
  return(manhattan_plot)
}

manhattan_plot_Multi <- function(GWAS,cutoff=0.05,QTN_index2 = NULL,QTN_index3=NULL,ZZ_label=NULL,ZZ_circle=NULL,Sig_L=NULL,Sig_Multi=NULL, trait = NULL,labels=NULL,model="unknown",Nchr=22,Multi=NULL,Multi_labels=NULL,correction="FDR",model_cor="BLINK",Third_Labels=NULL,Trait_one="Trait 1",Trait_two="Trait 2",Comparison="Population",Multi2=NULL,Multi2_labels=NULL){
  
  #library(compiler)
  #source("http://zzlab.net/GAPIT/GAPIT.library.R")
  #source("http://zzlab.net/GAPIT/gapit_functions.txt")
  library(GAPIT3)
  #find cutoff values using FDR or Bonferonni Cutfoff
  if(!is.null(Sig_L)){
    if(model_cor=="BLINK"){
      if(correction=="FDR"){
        cutoff.final=max(Sig_L$P.value)
      }else{
        cutoff.final=cutoff/length(GWAS$P.value)
      }}
    
    if(model_cor=="MLM"){
      if(correction=="FDR"){
        cutoff.final=max(Sig_L$P.value)
      }else{
        cutoff.final=cutoff/length(GWAS$P.value)
      }}
    
    SG_pval <- -log10(cutoff.final)
  }
  
  if(!is.null(Sig_Multi)){
    if(model_cor=="BLINK"){
      if(correction=="FDR"){
        cutoff.final.Multi=max(Sig_Multi$P.value)
      }else{
        cutoff.final.Multi=cutoff/length(Multi$P.value)
      }}
    
    if(model_cor=="MLM"){
      if(correction=="FDR"){
        cutoff.final.Multi=max(Sig_Multi$P.value)
      }else{
        cutoff.final.Multi=cutoff/length(Multi$P.value)
      }}
    
    SG_pval_Multi <- -log10(cutoff.final.Multi)
  }
  
  if(is.null(Sig_L)){
    cutoff.final=cutoff/length(GWAS$P.value)
    SG_pval <- -log10(cutoff.final)
  }
  if(is.null(Sig_Multi)){
    cutoff.final.Multi=cutoff/length(Multi$P.value)
    SG_pval_Multi <- -log10(cutoff.final.Multi)
    
  }
  
  
  col.Oceanic=rep(c(  '#EC5f67',    '#FAC863',  '#99C794',    '#6699CC',  '#C594C5'),ceiling(Nchr/5))
  
  #GWAS$Chromosome=as.character(GWAS$Chromosome)
  #GWAS$Position=as.numeric(GWAS$Chromosome)
  
  don <- GWAS %>% 
    # Compute chromosome size
    group_by(Chromosome) %>% 
    summarise(chr_len=max(Position)) %>% 
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(GWAS, ., by=c("Chromosome"="Chromosome")) %>%
    # Add a cumulative position of each SNP
    arrange(Chromosome, Position) %>%
    mutate( BPcum=Position+tot)
  
  if(!is.null(labels)){
    don=don%>% mutate( is_annotate=ifelse(SNP %in% labels, "yes", "no"))
    don=don%>%mutate( is_highlight=ifelse(P.value<=cutoff.final, "yes", "no"))
    SNP_don= don %>% filter(SNP %in% labels$SNP)
  }else{
    don$is_annotate=NA
    don$is_highlight=NA
  }
  
  #Then we need to prepare the X axis. Indeed we do not want to display the cumulative position of SNP in bp, but just show the chromosome name instead.
  axisdf = don %>% group_by(Chromosome) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  don$Trait=rep(Trait_one,length(don$SNP))
  
  
  if(!is.null(Multi)){
    #Multi$Chromosome=as.character(Multi$Chromosome)
    #Multi$Position=as.numeric(Multi$Chromosome)
    don_Multi <- Multi %>% 
      # Compute chromosome size
      group_by(Chromosome) %>% 
      summarise(chr_len=max(Position)) %>% 
      # Calculate cumulative position of each chromosome
      mutate(tot=cumsum(chr_len)-chr_len) %>%
      select(-chr_len) %>%
      # Add this info to the initial dataset
      left_join(Multi, ., by=c("Chromosome"="Chromosome")) %>%
      # Add a cumulative position of each SNP
      arrange(Chromosome, Position) %>%
      mutate( BPcum=Position+tot)
    if(!is.null(Multi_labels)){
      don_Multi=don_Multi%>% mutate( is_annotate=ifelse(SNP %in% Multi_labels$SNP, "yes", "no"))
      don_Multi=don_Multi%>%mutate( is_highlight=ifelse(P.value<=cutoff.final.Multi, "yes", "no"))
      SNP_don_Multi= don_Multi %>% filter(SNP %in% Multi_labels$SNP)
      
    }else{
      don_Multi$is_annotate=NA
      don_Multi$is_highlight=NA
    }
  }
  
  if(!is.null(Multi)){
    #Then we need to prepare the X axis. Indeed we do not want to display the cumulative position of SNP in bp, but just show the chromosome name instead.
    axisdf = don_Multi %>% group_by(Chromosome) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
    
    don_Multi$Trait=rep(Trait_two,length(don_Multi$SNP))
    
    
    don_Both=rbind(don,don_Multi)
    don_Both$Chr_Trait=paste(don_Both$Chromosome,don_Both$Trait)
    #don_Both=mutate(don_Both,Chr_Trait,paste(chromosome,Trait))
  }
  
  if(!is.null(Multi_labels)){
    SNP_both=rbind(SNP_don,SNP_don_Multi)
    both_labels=intersect(labels$SNP,Multi_labels$SNP)
    SNP_Multi= don_Both %>% filter(SNP %in% both_labels)
  }else{
    SNP_Multi=data.frame()
  }
  
  if(is.null(Multi)){
    
    manhattan_plot <- ggplot(don, aes(x=BPcum, y=-log10(P.value))) +
      # Show all points
      geom_point( aes(color=as.factor(Chromosome)), alpha=0.8, size=1.3) +
      scale_color_manual(values = col.Oceanic) +
      
      # custom X axis:
      scale_x_continuous( label = axisdf$Chromosome, breaks= axisdf$center ) +
      scale_y_continuous(limits = c(0,  12),labels =scales::number_format(accuracy=1)   ) +
      scale_shape_manual(values=c(21))+# remove space between plot area and x axis
      #scales::number_format(accuracy = 0)   
      #ylim(0,  max(-log10(don$P.value)))+
      # Custom the theme:
      geom_hline(yintercept= SG_pval,color="black")+
      labs(#title = paste("GWAS manhattan plot for trait:", as.character(trait)),
        y = "-log10(p)",
        x = "Chromosome",
        color = "Chromosome",
        tag=model)+
      theme_bw() +
      theme( 
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.tag=element_text(angle=-90,size=10),
        plot.tag.position=c("right")
      )
    if(!is.null(QTN_index3)){
      QTN=don%>%filter(SNP %in% QTN_index3)
      manhattan_plot=manhattan_plot+geom_vline(xintercept = QTN$BPcum, color = "red",linetype="solid") }
    if(!is.null(QTN_index2)){
      QTN=don%>%filter(SNP %in% QTN_index2)
      manhattan_plot=manhattan_plot+geom_vline(xintercept = QTN$BPcum, color = "green",linetype="dashed") }
    if(!is.null(trait)){manhattan_plot=manhattan_plot+labs(title = paste("GWAS manhattan plot for ", as.character(trait))) }
    if(!is.null(labels)){
      manhattan_plot=manhattan_plot+
        # Add highlighted points
        geom_point(data=subset(don, is_highlight=="yes"),size=3,shape=21) +
        # Add label using ggrepel to avoid overlapping
        geom_text( data=SNP_don, aes(label=SNP, size=3), nudge_y = 5)
    }
    
    if(!is.null(Third_Labels)){
      
      don_both_third=don%>%filter(SNP %in% Third_Labels)
      manhattan_plot=manhattan_plot+ geom_label_repel(data=don_both_third, aes(label=SNP), size=3) }
    
    
    if(!is.null(ZZ_label)){
      
      ZZ=don%>%filter(SNP %in% ZZ_label)
      manhattan_plot=manhattan_plot+ geom_text_repel( data=ZZ, aes(label=SNP, size=3,y = 10))  }
    
  }else{
    col.dual=rep(c(  '#194CEC','#B10913', '#349EF7','#F5414C'),ceiling(Nchr*2/4))
    col.dual.shape=rep(c(1,2),ceiling(Nchr*2/2))
    levels(as.factor(don_Both$Chr_Trait))
    manhattan_plot <- ggplot(don_Both, aes(x=BPcum, y=-log10(P.value))) +
      # Show all points
      geom_point( aes(shape= as.factor(Trait), color=as.factor(Chr_Trait)), alpha=0.8, size=1.3) +
      scale_color_manual(values = col.dual) +
      scale_shape_manual(values=c(21, 24))+
      #scale_shape_manual(values=col.dual.shape)+
      
      # custom X axis:
      scale_x_continuous( label = axisdf$Chromosome, breaks= axisdf$center ) +
      scale_y_continuous(limits = c(0,  12),labels =scales::number_format(accuracy=1)   ) +     # remove space between plot area and x axis
      #scales::number_format(accuracy = 0)   
      #ylim(0,  max(-log10(don$P.value)))+
      # Custom the theme:
      geom_hline(yintercept= SG_pval,color="black")+
      geom_hline(yintercept= SG_pval_Multi,color="black",linetype="dashed")+
      
      labs(#title = paste("GWAS manhattan plot for trait:", as.character(trait)),
        y = "-log10(p)",
        x = "Chromosome",
        color = Comparison,
        tag=model)+
      theme_bw() +
      theme( 
        #legend.position=c(1,1), legend.justification=c(1,1),legend.background=element_rect(fill="white", colour="black"),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.tag=element_text(angle=-90,size=10),
        plot.tag.position=c("right"))
    if(!is.null(QTN_index3)){
      QTN=don%>%filter(SNP %in% QTN_index3)
      manhattan_plot=manhattan_plot+geom_vline(xintercept = QTN$BPcum, color = "red",linetype="solid") }
    if(!is.null(QTN_index2)){
      QTN=don%>%filter(SNP %in% QTN_index2)
      manhattan_plot=manhattan_plot+geom_vline(xintercept = QTN$BPcum, color = "green",linetype="dashed") }
    if(!is.null(ZZ_label)){
      ZZ=don%>%filter(SNP %in% ZZ_label)
      manhattan_plot=manhattan_plot+ geom_text_repel( data=ZZ, aes(label=SNP,y = 10)) }
    if(!is.null(ZZ_circle)){
      ZZ=don%>%filter(SNP %in% ZZ_circle)
      manhattan_plot=manhattan_plot+ geom_point(data=ZZ,shape=1, size=5,color="green") }
    if(!is.null(trait)){manhattan_plot=manhattan_plot+labs(title = paste("GWAS manhattan plot for ", as.character(trait))) }
    if(!is.null(labels)){
      if(!is.null(Multi_labels)){
        manhattan_plot=manhattan_plot+
          # Add highlighted points
          geom_point(data=subset(don_Both, is_highlight=="yes"),aes(shape= as.factor(Trait), color=as.factor(Chr_Trait)), size=3) +
          # Add label using ggrepel to avoid overlapping
          geom_text_repel( data=SNP_both, aes(label=SNP, size=3))
      }else{
        manhattan_plot=manhattan_plot+
          # Add highlighted points
          geom_point(data=subset(don_Both, is_highlight=="yes"),aes(shape= as.factor(Trait), color=as.factor(Chr_Trait)), size=3) +
          # Add label using ggrepel to avoid overlapping
          geom_text( data=SNP_don, aes(label=SNP, size=3), nudge_y = 5)
      }
      
      # Add highlighted points
      #geom_point(data=subset(don, is_highlight=="yes"), size=3) +
      # Add label using ggrepel to avoid overlapping
      #geom_text_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=3)+
      
      #geom_point(data=subset(don_Multi, is_highlight=="yes"), size=3) +
      # Add label using ggrepel to avoid overlapping
      #geom_text_repel( data=subset(don_Multi, is_annotate=="yes"), aes(label=SNP), size=3)
    }
    if(!is.null(Third_Labels)){
      
      don_both_third=don%>%filter(SNP %in% Third_Labels)
      manhattan_plot=manhattan_plot+ geom_label_repel(data=don_both_third, aes(label=SNP), size=3) }
    if(!nrow(SNP_Multi)==0){manhattan_plot=manhattan_plot+geom_point(data=SNP_Multi,shape=1, size=5,color="green")
    
    
    }
  }
  
  return(manhattan_plot)
}

manhattan_plot_Single <- function(GWAS,cutoff=NULL,QTN_index = NULL, trait = NULL,labels=NULL,model="unknown",Nchr=22,Multi_labels=NULL,Third_Labels=NULL,Trait_one="Trait 1",Trait_two="Trait 2",Comparison="Population",Multi2_labels=NULL){
  #cutoff=NULL;
  #GWAS=FarmCPU_CV_Col_qam_adj_23$GWAS
  #labels =Bon_sig_F6_CV_Col$SNP
  #model="CL";
  #QTN_index=Bon_sig_B1_bl_CV_Col$SNP[2];
  #Multi_labels=c(Bon_sig_F1_bl_CV$SNP,Bon_sig_F1_bl$SNP);
  #Third_Labels=Bon_sig_B1_Col_Len_blup$SNP;
  
  m=nrow(GWAS)
  cutoff.final=(ifelse(
    is.null(cutoff),
    0.05/m,
    cutoff
  ))
  SG_pval <- -log10(cutoff.final)
  
  col.Oceanic=rep(c(  '#EC5f67',    '#FAC863',  '#99C794',    '#6699CC',  '#C594C5'),ceiling(Nchr/5))
  
  #GWAS$Chromosome=as.character(GWAS$Chromosome)
  #GWAS$Position=as.numeric(GWAS$Chromosome)
  
  don <- GWAS %>% 
    # Compute chromosome size
    group_by(Chromosome) %>% 
    summarise(chr_len=max(Position)) %>% 
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(GWAS, ., by=c("Chromosome"="Chromosome")) %>%
    # Add a cumulative position of each SNP
    arrange(Chromosome, Position) %>%
    mutate( BPcum=Position+tot)
  
  if(!is.null(labels)){
    don=don%>% mutate( is_annotate=ifelse(SNP %in% labels, "yes", "no"))
    don=don%>%mutate( is_highlight=ifelse(P.value<=cutoff.final, "yes", "no"))
  }
  
  #Then we need to prepare the X axis. Indeed we do not want to display the cumulative position of SNP in bp, but just show the chromosome name instead.
  axisdf = don %>% group_by(Chromosome) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  #don$Trait=rep(Trait_one,length(don$SNP))
  
  Multi_labels=unique(Multi_labels)
  both_labels=intersect(labels,Multi_labels)
  
  SNP_Multi= don %>% filter(SNP %in% both_labels)
  
  
  manhattan_plot <- ggplot(don, aes(x=BPcum, y=-log10(P.value))) +
    # Show all points
    geom_point( aes(color=as.factor(Chromosome)), alpha=0.8, size=1.3) +
    scale_color_manual(values = col.Oceanic) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$Chromosome, breaks= axisdf$center ) +
    scale_y_continuous(limits = c(0,  12),labels =scales::number_format(accuracy=1)   ) +     # remove space between plot area and x axis
    #scales::number_format(accuracy = 0)   
    #ylim(0,  max(-log10(don$P.value)))+
    # Custom the theme:
    geom_hline(yintercept= SG_pval,color="green")+
    #geom_hline(yintercept= SG_pval_Multi,color="green")+
    labs(#title = paste("GWAS manhattan plot for trait:", as.character(trait)),
      y = "-log10(p)",
      x = "Chromosome",
      color = Comparison,
      tag=model)+
    theme_bw() +
    theme( 
      #legend.position=c(1,1), legend.justification=c(1,1),legend.background=element_rect(fill="white", colour="black"),
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.tag=element_text(angle=-90,size=10),
      plot.tag.position=c("right")
    )
  if(!is.null(QTN_index)){
    QTN=don%>%filter(SNP %in% QTN_index)
    manhattan_plot=manhattan_plot+geom_vline(xintercept = QTN$BPcum, color = "red") }
  if(!is.null(trait)){manhattan_plot=manhattan_plot+labs(title = paste("GWAS manhattan plot for ", as.character(trait))) }
  if(!is.null(labels)){
    manhattan_plot=manhattan_plot+
      # Add highlighted points
      geom_point(data=subset(don, is_highlight=="yes"),color="orange", size=3) +
      # Add label using ggrepel to avoid overlapping
      geom_text_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=3)}
  
  
  # Add highlighted points
  #geom_point(data=subset(don, is_highlight=="yes"), size=3) +
  # Add label using ggrepel to avoid overlapping
  #geom_text_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=3)+
  
  #geom_point(data=subset(don_Multi, is_highlight=="yes"), size=3) +
  # Add label using ggrepel to avoid overlapping
  #geom_text_repel( data=subset(don_Multi, is_annotate=="yes"), aes(label=SNP), size=3)
  if(!is.null(Multi_labels)){
    
    don_both_second=don%>%filter(SNP %in% Multi_labels)
    manhattan_plot=manhattan_plot+ geom_text_repel(data=don_both_second, aes(label=SNP), size=3) }
  
  if(!is.null(Third_Labels)){
    
    don_both_third=don%>%filter(SNP %in% Third_Labels)
    manhattan_plot=manhattan_plot+ geom_label_repel(data=don_both_third, aes(label=SNP), size=3) }
  if(!nrow(SNP_Multi)==0){manhattan_plot=manhattan_plot+geom_point(data=SNP_Multi, color="green", size=5)}
  
  
  
  
  return(manhattan_plot)
}
