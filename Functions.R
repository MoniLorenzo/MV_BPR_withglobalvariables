###########Functions to include###############################
#Author Lorenzo Moni
##############################################################

################################################################################
######################## MULTI-VIEW BPR#########################################
################################################################################

###FUNCTION TO INCLUDE##########################


##########################
#2. Sample from categorical
#(ountput numeric from 1 to n_cat)
##########################
SampleMult1=function(Probs){
  K=length(Probs)
  return(sample.int(K, 1, prob = Probs))
}


##########################
#3. Sample/Update   Dir
##########################
SampleDir=function(a ){
  l=sum(!is.na(a) )
  r=a
  N=rgamma(l,shape = a)  
  
  r[1:l]= N/sum(N, na.rm = TRUE)  
  
  return(r)
}

##################################
#4. update alpha
##########################
Update_alpha=function(K_nuberCluster, nobs, alphaold, Hyp2alpha){
  a=Hyp2alpha["a"]
  b=Hyp2alpha["b"]
  
  log_epsilon=log(rbeta(1,alphaold+1, nobs))
  
  
  log_pi_eta=log(a+K_nuberCluster-1)-log(nobs*b-nobs*log_epsilon+a+K_nuberCluster-1)
  #draw from ber(explogpi)
  if (log_pi_eta>log(runif(1))){
    #if c==1
    return(rgamma(1,a+K_nuberCluster, b-log_epsilon))
  }else{
    return(rgamma(1,a+K_nuberCluster-1, b-log_epsilon))
  }
  
}

##################################
#5. UPDATE Z_i for irrelevant view Neal's 2nd alg.
##########################
Update_irrelevantview_Zi=function(X_is_ds,#=X[,Gamma_v, drop=FALSE],
                                  n_,#=N, 
                                  alpha_v,#=1,
                                  A_v_,#= NULL,#unname(A_v_Nk[,"A_v"]),
                                  #   Hyper_a=hyp_a_mat,
                                  Prior_pred_prod_philog_i,#=rowSums(LogPrior_pred_phi_i[,Gamma_v, drop=FALSE]),
                                  Internal_state_log,#=lapply(STATE[["view1"]], function(x) list(Phi=log(x$Phi)[Gamma_v, ,drop=FALSE ], n_k=x$n_k)) ,                                #Phi_ks_list=STATE[[paste0("view",view)]][["Phi"]],
                                  Z_old,#=STATE$Clusterall[, view] ,
                                  G0_sample_log){#=function() {log(apply(hyp_a_mat[Gamma_v, ,drop=FALSE ],1, SampleDir)) }
  
  
  
  indexN_var_view_v=seq_along(X_is_ds[1,])#1:ncol(X_is_ds)
  
  # print(nrow(X_is_ds)==n_)
  K_n_=  length(Internal_state_log)
  newlabels=A_v_
  
  for (i in 1:n_) {
    #   print(newlabels)
    #  cat("\n")
    #check if n_-i,k==1 if yes rm the k st k==z_^v[i] component or decrease n_k
    char_Zi=paste0("k",Z_old[i])
    n_k_i=Internal_state_log[[char_Zi]]$n_k
    
    #  cat("s=",s," i-",i, " zi",Z_old[i], " ni=", n_k_i, "\n")
    
    if(n_k_i==1){
      Internal_state_log[[char_Zi]]=NULL
      K_n_=K_n_-1
      #relabeling
      newlabels=unique(Z_old[-i])
      
      Z_old=match(Z_old, newlabels)
      names(Internal_state_log)= paste0("k",
                                        match(as.numeric(gsub("\\D", "", 
                                                              names(Internal_state_log))),newlabels ))
      newlabels=as.numeric(gsub("\\D", "", names(Internal_state_log)))
      #  print("Relabel")
    }else{
      Internal_state_log[[char_Zi]]$n_k=n_k_i-1
    }
    
    #existing component
    Probs_Existing_log=sapply(Internal_state_log, function(x){
      log(x$n_k)+sum(x$Phi[cbind(indexN_var_view_v, X_is_ds[i,])])
    })
    
    #new component
    Prob_new_log= log(alpha_v)+Prior_pred_prod_philog_i[i]
    
    Probs_unorm=c(Probs_Existing_log,Prob_new_log)
    maxi=max(Probs_unorm)
    
    #sample component ie z^v_i
    Z_New_i=SampleMult1(Probs = exp(Probs_unorm-maxi) )
    
    
    if(Z_New_i==K_n_+1){#add
      Z_old[i]=Z_New_i 
      Internal_state_log[[paste0("k",Z_old[i])]]=list(Phi=G0_sample_log(),
                                                      n_k=1)
      K_n_=K_n_+1
      newlabels=c(newlabels,Z_New_i)
      
      #   print("here")
      #  print(i)
    }else{#
      Z_old[i]=newlabels[Z_New_i] 
      
      Internal_state_log[[paste0("k",Z_old[i])]]$n_k=Internal_state_log[[paste0("k",Z_old[i])]]$n_k+1
    }
    #check
    if(any(newlabels!=as.numeric(gsub("\\D", "", 
                                      names(Internal_state_log))))){
      cat("\n PROB \n")
      stop()
    }
    
  }
  
  return(  list(Z_updated=Z_old, Phis_inviewV_list=Internal_state_log, 
                Av=newlabels, K_v=K_n_,
                N_k=sapply(Internal_state_log, function(x) x$n_k)))
} 





############################ ####################################################
#6. Update gamma_d's in case of variable selection (only view 0 and 1)
##############################################################################
View_allocations_update=function(X, nvariables, 
                                 Z_=STATE$Clusterall[,1],# Z view 1 STATE view 1
                                 STATE_clusterspecpar=STATE[[paste0("view",view)]],
                                 PHI_NULL=STATE$Phi0,
                                 Succ_d0=S0_all_dmat,
                                 Hyperpar_log=log(c(nu1=0.5, nu0=0.5))){
  
  
  
  clusterspec_logs=sapply(STATE_clusterspecpar, function(x) log(x$Phi), 
                          simplify = "array"  )
  
  Matrix=apply(cbind(X,Z_),1, function(x, PAR=clusterspec_logs){
    PAR[ , , paste0("k",x[nvariables+1])][cbind(1:nvariables, x[1:nvariables]), drop=FALSE]
    
  })
  
  log_view1probs=rowSums(Matrix)+Hyperpar_log["nu1"]
  
  log_view0prob=rowSums(Succ_d0*log(PHI_NULL))+Hyperpar_log["nu0"]
  
  Matrix2=cbind(log_view0prob, log_view1probs)
  Gammas=apply(Matrix2, 1, function(x){
    if(x[1]>x[2]){ 
      y=x-x[1]
      SampleMult1(exp(y)/sum(exp(y)))-1
    }else{
      y=x-x[2]
      SampleMult1(exp(y)/sum(exp(y)))-1
    }}
  )
  
  names(Gammas)=NULL
  return(Gammas)
}


############################ ####################################################
#6. Update gamma_d's in case of Multi-view BPR (only 3 view: 0(null) 1(relevant) 2(irrelevant))
##############################################################################
View_allocations_update_3view=function(X, nvariables=20, 
                                       Z_1=STATE$Clusterall[,1],
                                       Z_2=STATE$Clusterall[,2],
                                       
                                       STATE_clusterspecpar1=STATE[[paste0("view",1)]],
                                       STATE_clusterspecpar2=STATE[[paste0("view",2)]],
                                       Succ_d0=S0_all_dmat,
                                       PHI_NULL=STATE$Phi0,
                                       Hyperpar_log=log(c(nu1=1/3, nu0=1/3, nu2=1/3))){
  
  
  
  clusterspec_logs1=sapply(STATE_clusterspecpar1, function(x) log(x$Phi), simplify = "array")
  clusterspec_logs2=sapply(STATE_clusterspecpar2, function(x) log(x$Phi), simplify = "array")
  
  
  Matrix1=apply(cbind(X,Z_1),1, function(x, PAR=clusterspec_logs1){
    PAR[ , , paste0("k",x[nvariables+1])][cbind(1:nvariables, x[1:nvariables]), drop=FALSE]
    
  })
  Matrix2=apply(cbind(X,Z_2),1, function(x, PAR=clusterspec_logs2){
    PAR[ , , paste0("k",x[nvariables+1])][cbind(1:nvariables, x[1:nvariables]), drop=FALSE]
    
  })
  
  log_view2probs=rowSums(Matrix2)+Hyperpar_log["nu2"]
  
  
  log_view1probs=rowSums(Matrix1)+Hyperpar_log["nu1"]
  
  log_view0prob=rowSums(Succ_d0*log(PHI_NULL))+Hyperpar_log["nu0"]
  
  Matrix_PROBS=cbind(log_view0prob, log_view1probs, log_view2probs)
  Gammas=apply(Matrix_PROBS, 1, function(x){
    
    y=x-max(x)
    SampleMult1(exp(y)/sum(exp(y)))-1
  }
  )
  
  names(Gammas)=NULL
  return(Gammas)
}



############################ ####################################################
#6. sample from the predictive of Y: given  for a new predictive units GIVEN
#                                    their profiles (x) and prognostic variables (w)
#Note these values must be in a matrix of size Nnewpredictiveunits X numb. covariates 
##############################################################################
Predictive_Y=function(Matrix_pred_x=Pred_dataset_X ,
                      npred=seq_along( Pred_dataset_X[,1]),
                      Matrix_pred_W=Pred_dataset_W, 
                      ry=R_y,
                      Gammas_=STATE$Gamma,
                      Beta_m=Mat_betas,
                      Alpha_v1=STATE$Alpha[1],
                      HY_theta=Hyperparameters$theta,
                      Log_prior_PRED_x=LogPrior_pred_phi_FOR_predictive_units,
                      state_view1=STATE$view1){
  
  #all units
  Gamma_logic= Gammas_==1
  
  Lp_all=Matrix_pred_W%*%Beta_m
  
  
  state_view1[["New"]]=NA
  
  #matrix thetas tilde ry-1  x Npred  
  ThetaS_tilde=  sapply(npred, function(ii){ 
    x_tilde=Matrix_pred_x[ii,Gamma_logic]
    
    #existing
    Probs_existing_new=sapply(state_view1, Internal_pred, Gl=Gamma_logic, xtl=x_tilde)
    Probs_existing_new["New"]=sum(Log_prior_PRED_x[ii,Gamma_logic])+log(Alpha_v1)
    maxi=max(Probs_existing_new)
    
    z_tilde=sample( names(Probs_existing_new ) , size = 1, 
                    prob = exp(Probs_existing_new-maxi))
    #print(z_tilde)
    if(z_tilde!="New"){
      theta_tilde=state_view1[[z_tilde]]$Theta
    }else{
      theta_tilde=G_0theta(hyperparameters = HY_theta, RY = ry)
    }
    
  })
  #draw from Y
  etas=compute_log_etas_thetaMAT_cpp(Exp = TRUE, Linpred_mat = Lp_all,
                                     THETA_mat = t(ThetaS_tilde) )
  
  
  Y_tilde=apply(etas, 1, SampleMult1)
  
}



Internal_pred=function(k, Gl=Gamma_logic, xtl=x_tilde){
  if(is.na(k)[1] ){return(NA)}
  log(k$n_k)+sum(log(k$Phi[cbind(which(Gl),xtl )]))
}


##################################################

##################################
##SAMPLE FROM THETA PRIOR
##################################
G_0theta=function(hyperparameters, RY   ){
  R_Ym1=RY-1
  
  mu=hyperparameters[["mu"]]
  Tr=rt(R_Ym1,df=hyperparameters[["df"]])
  
  
  return(hyperparameters[["mu"]]+hyperparameters[["scale"]]*Tr)
  
}



##################################
# COMPUTE PROBS. RESPONSE  
##Compute etas, for each unit (rows Linpred), the categories are the col
##################################
#  compute_log_etas=function(Exp=TRUE, 
#                                   #n_P_W=n_cov_W
#                                  Linpred_mat, # Wi%*%BETA, 
#                                  THETA#=c(-5,3) #MATRIX OR VECTOR
#                                  ){
# N=dim(Linpred_mat)[1]
# Rym1=dim(Linpred_mat)[2]
# #Computer unit specific probabilities 
# #check
# #if(dim(Linpred_mat)!=c(N,length(THETA))){stop()}
# 
# if(is.matrix(THETA)){
#   thetamat=THETA
#   #print("matrix")
# }else{ 
#   #print("vector")
#   thetamat=matrix(THETA, nrow = dim(Linpred_mat)[1], 
#                 ncol =dim(Linpred_mat)[2], byrow = TRUE )
# }
# #Matrix of 
# #theta_+_r +w_i^T%*%Beta_r, NOTE:Linpred_mat=w_i^T%*%Beta_r
# #by row i, by col r=1,..., Ry-1
# Nums_1_to_rym1=thetamat+Linpred_mat
# #add last category lin pred with theta, ie add col of 0
# Nums_1_to_ry=cbind(Nums_1_to_rym1, 0)
# 
# #compute the probs. to improve the stability we subtract the max ie ach row
# etas_i=t(apply(Nums_1_to_ry, 1, FUN = function(x){
#                                               y=exp(x-max(x))
#                                               y/sum(y)}) )
#
# if(!Exp){ log(etas_i)}else{etas_i}
# 
# } 

########  
#no  
#compute_log_etas_thetamatrix=function(#N,
#    exp=TRUE, 
#    #n_P_W=n_cov_W
#    Linpred_mat, # Wi%*%BETA, 
#    THETA_mat#=c(-5,3)
#  ){
#    N=dim(Linpred_mat)[1]
#    Rym1=dim(Linpred_mat)[2]
#    #Computer unit specific probabilities 
#    #check
#    dim(Linpred_mat)==dim(Linpred_mat)
#    
#    thetamat=THETA_mat
#    #Matrix of 
#    #theta_+_r +w_i^T%*%Beta_r, NOTE:Linpred_mat=w_i^T%*%Beta_r
#    #by row i, by col r=1,..., Ry-1
#    Nums_1_to_rym1=thetamat+Linpred_mat
#    #add last category lin pred with theta, ie add col of 0
#    Nums_1_to_ry=cbind(Nums_1_to_rym1, 0)
#    
#    #compute the probs. to improve the stability we subtract the max ie ach row
#    etas_i=t(apply(Nums_1_to_ry, 1, FUN = function(x){
#      y=exp(x-max(x))
#      y/sum(y)}) )
#    
#    if(!exp){ log(etas_i)}else{etas_i}
#    
#  }   
# 

##########################################
# Prior/log prior independent LS t-student  
###########################################
Prior_LST_dens_equalHyper_log=function(x,  
                                       exp=FALSE, #
                                       mu=hyperparameters[["mu"]],
                                       scale=hyperparameters[["scale"]],
                                       df=hyperparameters[["df"]] 
){
  #NOTE scale=sigma*sqrt(df-2/df)
  #THIS FUNCTION OMIT THE NORMALIZING CONSTANT
  #K is the number of random variables
  Normkonst=length(x)*(log(gamma(df/2+0.5))-log(gamma(df/2))-0.5*log(pi*df)-log(scale))
  
  
  frac=(x^2-2*mu*x+mu^2)/scale^2 
  complogf=log(1+frac/df)
  LOGf=-(df/2+0.5)*sum(complogf) 
  #return
  if(!exp){ LOGf}else{exp(LOGf)}
}  

##################################
#updating theta
#################################
##compute internal state for multiview
#Thetas are not log
Compute_Internal_stateview1=function(Oldstate=STATE[["view1"]], Gamma_ind=Gamma_v){
  lapply(Oldstate, function(x){
    if(!is.null(x$Theta)){tet=x$Theta
    }else{tet=G_0theta(hyperparameters = Hyperparameters$theta,RY=R_y); print("Generate_new_theta")}
    
    list(Phi=log(x$Phi)[Gamma_ind, ,drop=FALSE ], 
         Theta=tet, n_k=x$n_k )})
}

###########################################
####UPDATE Z_i for  view: USING RCPP etas (modified algoritm 8 Neal)
############################################
Update_Relevant_Zi=function(Y ,
                            RY=R_y,
                            X_is_ds =X[,Gamma_v, drop=FALSE], 
                            n_ =N, 
                            alpha_v=1,#=1,
                            LINPREDBetas=matrix(0, nrow = N, ncol = R_y-1),#no covariate
                            
                            A_v_=unname(STATE$relab[[paste0("view",view)]][,"A_v"]),#= NULL,#unname(A_v_Nk[,"A_v"]),
                            
                            M=10,
                            
                            G0Y_sample=function(x) G_0theta(hyperparameters = Hyperparameters$theta, RY=R_y),
                            
                            X_Prior_pred_prod_philog_i=rowSums(LogPrior_pred_phi_i[,Gamma_v, drop=FALSE]), 
                            G0X_sample_log=function() {log(t(apply(Hyperparameters$clusteringvariable[Gamma_v, ,drop=FALSE ],1, SampleDir)) )},
                            
                            Internal_state_log =Compute_Internal_stateview1(Oldstate = STATE[["view1"]], Gamma_ind = Gamma_v),
                            Z_old=STATE$Clusterall[, view] 
){#=function() {log(apply(hyp_a_mat[Gamma_v, ,drop=FALSE ],1, SampleDir)) }
  
   
  
  indexN_var_view_v=seq_along(X_is_ds[1,])#1:ncol(X_is_ds)
  
  # print(nrow(X_is_ds)==n_)
  K_n_=  length(Internal_state_log)
  newlabels=A_v_#names(Z_old)
  
  
  
  
  for (i in 1:n_) {
    #    ttt=Sys.time()
    #Generate M additional parameters: only theta
    Matrix_additonal_thetas=sapply(1:M,  G0Y_sample) # Ry-1 x M, by row theta1 theta2
    #  print(Sys.time()-ttt)
    
    #   print(newlabels)
    #  cat("\n")
    #check if n_-i,k==1 if yes rm the k st k==z_^v[i] component or decrease n_k
    char_Zi=paste0("k",Z_old[i])
    n_k_i=Internal_state_log[[char_Zi]]$n_k
    
    #  cat("s=",s," i-",i, " zi",Z_old[i], " ni=", n_k_i, "\n")
    
    if(n_k_i==1){#
      #if i is the only unit associated to this  k:
      #then 1:, k^-= numbgroup -1; and rm from intstate 
      # 2: let theta_k be theta_{k^- + 1} in   List_additonal_thetas  
      theta_ki=Internal_state_log[[char_Zi]]$Theta
      Internal_state_log[[char_Zi]]=NULL
      K_n_=K_n_-1
      
      Matrix_additonal_thetas[,1]=theta_ki
      #relabeling the existing groupp
      newlabels=unique(Z_old[-i])
      Z_old=match(Z_old, newlabels)
      names(Internal_state_log)= paste0("k",
                                        match(as.numeric(gsub("\\D", "", 
                                                              names(Internal_state_log))),newlabels ))
      newlabels=as.numeric(gsub("\\D", "", names(Internal_state_log)))
      #  print("Relabel")
    }else{#DRAWN m new values from G0 theta
      Internal_state_log[[char_Zi]]$n_k=n_k_i-1
      
      
    }
    
 
    
    #existing component
    Probs_Existing_log=sapply(Internal_state_log, function(x){
      log(x$n_k)+sum(x$Phi[cbind(indexN_var_view_v, X_is_ds[i,])])+ #X and n
        compute_log_etas_cpp( Exp=FALSE,  #model Y
                              Linpred_mat = LINPREDBetas[i, , drop=FALSE], # vec 1xRy-1 of 0's
                              THETA =x$Theta)[, Y[i]]
    })
    
    #new component
    Prob_new_log= log(alpha_v)-log(M)+X_Prior_pred_prod_philog_i[i] +
      apply(Matrix_additonal_thetas, 2, function(x){
        compute_log_etas_cpp( Exp=FALSE,  #model Y
                              Linpred_mat = LINPREDBetas[i, , drop=FALSE], # vec 1xRy-1 of 0's
                              THETA =x)[, Y[i]]
      })
    # print(Sys.time()-ttt)
    
     
    Probs_unorm=c(Probs_Existing_log,Prob_new_log)
    maxi=max(Probs_unorm)
    
    #sample component ie z^v_i
    Z_New_i=SampleMult1(Probs = exp(Probs_unorm-maxi) )
    
    
    if(Z_New_i>=K_n_+1){#add
      # print("new")
      Z_old[i]=K_n_+1 #Z_New_i #new label K_n_
      Internal_state_log[[paste0("k",Z_old[i])]]=list(Phi=G0X_sample_log(),#<- ??? mettere oldPhi?
                                                      Theta=Matrix_additonal_thetas[, Z_New_i-K_n_],
                                                      n_k=1)
      K_n_=K_n_+1
      newlabels=c(newlabels,Z_old[i])
      
 
    }else{#
      Z_old[i]=newlabels[Z_New_i] 
      
      Internal_state_log[[paste0("k",Z_old[i])]]$n_k=Internal_state_log[[paste0("k",Z_old[i])]]$n_k+1
    }
    #check
    if(any(newlabels!=as.numeric(gsub("\\D", "", 
                                      names(Internal_state_log))))){
      cat("\n PROB \n")
      stop()
    }
    
  }
  
  return(  list( Z_updated=Z_old, Param_inviewV_list=Internal_state_log, 
                 Av=newlabels, K_v=K_n_,
                 N_k=sapply(Internal_state_log, function(x) x$n_k)))
} 




##################################
#Updating beta
##################################
#utilities for beta  
proposal_forbeta=function(xt_1, delta=0.1, 
                          V=NULL){
  xt_1=as.vector(xt_1)
  if(is.null(V)){V=diag(length(xt_1))}
  SIGMA= V*delta
  MASS::mvrnorm(1,mu=xt_1, Sigma = SIGMA)
  
}

#Not  nedded in sym MH
# density_proposal=function(x, xt_1, delta){
#    SIGMA= diag(length(xt_1))*delta
#   
#    mvtnorm::dmvnorm(x =x ,mean = xt_1, sigma =SIGMA, log = TRUE  )

#  }

#  
#sampler for beta
##################################
Update_beta=function(betaold=STATE$view1$beta, 
                     proposal_func=proposal_forbeta, 
                     Y_, #Y st z_i kth group
                     THETA_=c(0,0),##thetas
                     Wi_=Wi, 
                     Delta=1,
                     V=NULL,#, variance proposal
                     n.categ=3,
                     Hyperp # list "mu" "scale" "df", 
                     
){
  
  #this is vector
  Beta_prop=proposal_func(xt_1 = betaold, 
                          delta = Delta,V=V)
  Beta_prop_mat=matrix(Beta_prop, ncol=n.categ)
  
  
  LP_old=Wi_%*%betaold
  LP_prop=Wi_%*%Beta_prop_mat
  #compute ratio log proposal
  #Q(x_old|x_proposed)/Q(x_proposed|x_old)
  logprorratio=0   #<-IF Metropolis
  #compute logratio post(theta|thetaold)/post(thetaold|theta)
  logr=logprorratio+
    #Prior proposal
    Prior_LST_dens_equalHyper_log(x=Beta_prop, 
                                  exp=FALSE,  
                                  mu=Hyperp[["mu"]],
                                  scale=Hyperp[["scale"]],
                                  df=Hyperp[["df"]] ) +
    #lik proposal
    sum(compute_log_etas_thetaMAT_cpp(Linpred_mat = LP_prop, 
                                      THETA = THETA_, 
                                      Exp = FALSE)[cbind(seq_along(Y_), Y_)])-
    #prior old 
    Prior_LST_dens_equalHyper_log(x=as.vector(betaold), 
                                  exp=FALSE,  
                                  mu=Hyperp[["mu"]],
                                  scale=Hyperp[["scale"]],
                                  df=Hyperp[["df"]] )-
    #lik old
    # sum(Etas_old[cbind(1, Y_)])
    sum(compute_log_etas_thetaMAT_cpp(Linpred_mat =LP_old, 
                                      THETA = THETA_, 
                                      Exp = FALSE)[cbind(seq_along(Y_), Y_)])
  
  
  if (logr>log(runif(1))){
    list(beta=Beta_prop_mat, betavec=Beta_prop,
         acc=1, Lp=LP_prop) 
  }else{
    list(beta=betaold,  betavec=as.vector(betaold),
         acc=0, Lp=LP_old) 
    
  }
  
}





##################################
#Updating theta
##################################

#Proposal theta
proposal_fortheta=function(xt_1, delta){
  SIGMA= diag(length(xt_1))*delta#var(RES)#diag(2)*100
  MASS::mvrnorm(1,mu=xt_1, Sigma = SIGMA)
  
}


#updating theta
Update_theta=function(thetaold=STATE$view1$theta, 
                      proposal_func, 
                      Delta,
                      Y_, #Y st z_i kth group
                      Hyperp, # list "mu" "scale" "df", 
                      MATRIX_wBetas=LP){
  
  
  
  Theta_prop=proposal_func(xt_1 = thetaold, delta = Delta)
  #compute ratio log proposal
  #Q(x_old|x_proposed)/Q(x_proposed|x_old)
  logprorratio=0 #<-if metropolis
  
  #compute logratio post(theta|thetaold)/post(thetaold|theta)
  logr=logprorratio+
    #Prior proposal
    Prior_LST_dens_equalHyper_log(x=Theta_prop,  
                                  exp=FALSE,  
                                  mu=Hyperp[["mu"]],
                                  scale=Hyperp[["scale"]],
                                  df=Hyperp[["df"]] ) +
    #lik proposal
    sum(compute_log_etas_cpp(Linpred_mat =MATRIX_wBetas, 
                             THETA = Theta_prop, 
                             Exp = FALSE)[cbind(seq_along(Y_), Y_)])-
    #prior old 
    Prior_LST_dens_equalHyper_log(x=thetaold,  
                                  exp=FALSE,  
                                  mu=Hyperp[["mu"]],
                                  scale=Hyperp[["scale"]],
                                  df=Hyperp[["df"]] )-
    #lik old
    # sum(Etas_old[cbind(1, Y_)])
    sum(compute_log_etas_cpp(Linpred_mat =MATRIX_wBetas, 
                             THETA = thetaold, 
                             Exp = FALSE)[cbind(seq_along(Y_), Y_)])
  
  
  if (logr>log(runif(1))){
    list(theta=Theta_prop, acc=1) 
  }else{
    list(theta=thetaold, acc=0) 
    
  }
  
}
#-----end



























