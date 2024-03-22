############################
##MV-BPR 3 view; Y: categorical multinomial regression regression with prognostic variables
#                X: categorical variables                         
############################

############################
# Data and preliminaries processing
rm(list=ls())
  Filepath_data="/home/lorenzo/Dropbox/Testdataset/alternative/Sim3/Datasets/gData_3-10-7_rep10"#File path of csv containing dataframe with y, all x, all w 
  File_predictiveDS=NULL#"/home/lorenzo/Dropbox/Testdataset/alternative/Sim3/Datasets/gData_3-10-7Outsample"#File path of csv containing dataframe with  all predictive  x_tilda,   w_tilda 
 
  #^ #NOTE very importand the colomns containing the clustering variables must be nammed as 
  # d1 d2 d3 ...etc.  In the same wat the cols containing the the prognostic variable with
  # w1 w2 ....etc., the must be nammed as  y    
  
 #set.seed()
  
  
  library(Rcpp)
  
  gc()
 
  source("/home/lorenzo/Desktop/github/Functions.R")

  sourceCpp("/home/lorenzo/Desktop/github/etascpp.cpp")
  sourceCpp("/home/lorenzo/Desktop/github/etas_thetamat.cpp")
  
   
  
  
  
  #filepred=paste0("/home/lorenzo/Dropbox/Testdataset/alternative/Sim3/Datasets/",
   #               sub("_rep.*", "", fff), "Outsample")
  
  
   
  print(File_predictiveDS)
  print(Filepath_data)
  #####---
  # TRUE GAMMAS
  #READ_TRUE_GAMMAS=read.csv("/home/lorenzo/Dropbox/Testdataset/Sim3/Datasets/aRealgammas")
  #TRUE_GAMMAS=as.vector(READ_TRUE_GAMMAS[READ_TRUE_GAMMAS$X==gsub("_rep[0-9]+", "", fff),-1])==1
  
  
  print( c("Multiview \n" ))
  
  ### Predictive units 
  if(!is.null(File_predictiveDS)){
    
  Pred_dataset=read.csv(File_predictiveDS)
  
  Indeces_pred=1:nrow(Pred_dataset)  
  
  
  Pred_dataset_X=as.matrix(Pred_dataset[Indeces_pred, 
                                        grep("^d", colnames(Pred_dataset)),
                                        drop=FALSE]) 
  
  
  Pred_dataset_W=as.matrix(Pred_dataset[Indeces_pred , 
                                        grep("^w", colnames(Pred_dataset)),  
                                        drop=FALSE])  
  
  
  rownames(Pred_dataset_W)=rownames(Pred_dataset_X)=c(
                                                      1:nrow(Pred_dataset))
  
  rm(Pred_dataset)
  
  colnames(Pred_dataset_X)=paste0("d",1:D)
  
  }
  
  
  
  
  ##############-
  #MODEL DATA
  Data=read.csv(Filepath_data)
  
  
  
  #Data=nulldata
  #READ DATA 
  nview=3
  Y=Data[, "y"]
  
  X=Data[ , grep("^d", colnames(Data))] 
  #print(colnames(X))
  
  X=as.matrix(X)
  N=nrow(X)
  D=ncol(X)
  
  colnames(X)=paste0("d",1:D)
  
  
  #prognostic vriables
  
  W=as.matrix(Data[ , grep("^w", colnames(Data))])
  P_w=ncol(W)
  colnames(W)=paste0("w",seq_len(P_w))
  
  
  
  #number of category by clustering variables
  R_ds=apply(X, 2, function(x){max(x)})
  #number of category of outcome  
  R_y=length(unique(Y))
  
  R_betas=(R_y-1)*P_w
  
  #note:
  #Beta matrix 
  #BETA:=(   |       ...         |   )
  #      (beta_r1    ...    beta_ry-1)  
  #      (   |       ...         |   )
  
  #IF no W variables the linear predictor matrix is a null matrix n x R_y-1
  
  #LPnull=matrix(0, nrow = N, ncol = R_y-1)
  
  ####################################################################
  ##                     SET HYPERPARAMETERS                        ##
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  #---------------------CHANGE hyperparameters Here-------------------
  
  #SET HYPERPARAMETERS: In this formulation Phi_k,d=Phi_prior  (k=0,1,2,....,and d=1,...,D )
  
  Hyperparameters=list()
  hyp_a_mat=matrix(NA, ncol=max(R_ds), nrow = D)
  rownames(hyp_a_mat)=paste0("d",1:D)
  for (dd in 1:D){
    hyp_a_mat[dd,1:R_ds[dd]]= rep(1, R_ds[dd]) #Phi_prior: <- CHANGE HERE
    
  }
  
  ##Phis_k and phi_0 (mat containing all d) 
  Hyperparameters[["clusteringvariable"]]=hyp_a_mat
  
  #############hyperprior alpha~gamma(a,b) 
  #hyperhyperparameters gamma(a,b) of alphas FOR all view
  Hyperparameters[["hyperhyperprioralphas"]]["a"]=2  #<- CHANGE HERE
  Hyperparameters[["hyperhyperprioralphas"]]["b"]=1  #<- CHANGE HERE
  
  #-->>ADD hyperprior THETA BETAS GAMMAS  
  ##############Thetas (view 1)
  Hyperparameters[["theta"]][["mu"]]=0               #<- CHANGE HERE
  Hyperparameters[["theta"]][["scale"]]=2.5          #<- CHANGE HERE
  Hyperparameters[["theta"]][["df"]]=7              #<- CHANGE HERE
  
  ##############Betas
  Hyperparameters[["beta"]][["mu"]]=0               #<- CHANGE HERE
  Hyperparameters[["beta"]][["scale"]]=2.5          #<- CHANGE HERE
  Hyperparameters[["beta"]][["df"]]=7               #<- CHANGE HERE
  
  rm(hyp_a_mat)
  
  
  ###Gammas (NOTE for computational reason specify these hyperparameter on log scale)
  Hyperparameters[["hyperGamma_LOG"]]["nu0"]=log(1/3) #<- CHANGE HERE
  Hyperparameters[["hyperGamma_LOG"]]["nu1"]=log(1/3) #<- CHANGE HERE
  Hyperparameters[["hyperGamma_LOG"]]["nu2"]=log(1/3) #<- CHANGE HERE
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ####################################################################
  ##            COMPUTE SOME USEFUL (CONSTANT) QUANTITIES           ##
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  #1. Sufficient statistics array:
  #        D matrix of n x max(R_d), it indicate which category 
  #        of the d-th variable the is in i-th 
  S_d_comp=array(NA, dim = c(N,max(R_ds),   D))
  dimnames(S_d_comp)=list(NULL,NULL, paste0("d",1:D))
  for (dd in 1:D ){
    S_d_comp[1:N,1:R_ds[dd], dd ]=
      apply(matrix(1:R_ds[dd], ncol = R_ds[dd], nrow = N,
                   byrow = TRUE),
            2,
            function(x)   (x==X[,dd])*1)
  }
  
  #2. compute the matrix of all priorpredictive value for all units and d
  #HERE WE ASSUME HYPERPARAMETER EQUAL IN ALL VIEWs
  #this matrix has for each row the value of log pred for each unit
  #for the variable d (in col)
  LogPrior_pred_phi_i=t(apply(X,1, function(x, A_mat=Hyperparameters$clusteringvariable){
    logSumden=log(rowSums(A_mat, na.rm = TRUE))
    lognum_d=log(A_mat[ names(x), ][cbind(1:D, x)])
    return(lognum_d-logSumden)       
    
  }))
  
  #3. 
  #this matrix D x max(R_d) count the total success in each caegory for all D variables 
  S0_all_dmat=t(apply(S_d_comp, 3,   colSums ) ) 
  
  
  
  #4. same of 2. but for the predictive units
  if(!is.null(File_predictiveDS)){
    
  LogPrior_pred_phi_FOR_predictive_units=
    t(apply(Pred_dataset_X,1, 
            function(x, A_mat=Hyperparameters$clusteringvariable){
              logSumden=log(rowSums(A_mat, na.rm = TRUE))
              lognum_d=log(A_mat[ names(x), ][cbind(1:D, x)])
              return(lognum_d-logSumden)       
              
            }))
  }
  ####################################################################
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #----------------------#Initialize STATE MCMC----------------------#
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ####################################################################
  
  
  #INITIALIZE STATE LIST
  #Gamma# rep(1, D) and no updating in sampler if no variable selection
  Gamma=sample(c(0,1,2), size = D,replace = TRUE)
  
  #################
  #General state: Gammas Phi0, Z^(v)'s 
  STATE=list()
  STATE[["Gamma"]]=Gamma# c(rep(1,10),rep(0,10))
  STATE[["Phi0"]]=  t( apply(Hyperparameters$clusteringvariable,1, SampleDir)   )
  STATE[["Alpha"]]=rchisq(nview-1, df=2)
  
  STATE[["Betas"]]= c(0,0)
  
  #BETA:=(   |       ...         |   )
  #      (beta_r1    ...    beta_ry-1)  
  #      (   |       ...         |   ) X [P_W X RY-1]
  #AS.VECTOR(BETA)= (--beta_r1--, ..., --beta_ry--)
  
  
  
  #################
  ### view specific:
  #Allocation variables Z   #ncol=nview-1
  z_temp_view=list(); Z_mats=matrix(NA, ncol = nview-1, nrow = N)
  z_temp_view[[1]]=sample(rep(1:10, size = 1  ,length.out=N))
  z_temp_view[[2]]=sample(rep(1:10, size = 1  ,length.out=N))
  
  STATE[["Clusterall"]]=NA
  STATE[["relab"]]=list()
  
  #cluster specific parameters accordingly to Z^(v)
  for (view in 1:(nview-1)) {
    #Z: matrix containing the Z^(v)'s   
    Z_mats[,view]= z_temp_view[[view]]  
    #clust par  
    TAB=table( (Z_mats[,view])) 
    for (kkk in as.numeric(names(TAB))) {
      #given Z initialize state in each k
      ##phi_k
      STATE[[paste0("view",view)]][[paste0("k",kkk)]][["Phi"]]= 
        t( apply(Hyperparameters$clusteringvariable,1, SampleDir)   )   
      ##theta_k IF view=1
      if(view==1){
        STATE[[paste0("view",view)]][[paste0("k",kkk)]][["Theta"]]= 
          G_0theta(hyperparameters = Hyperparameters$theta, RY = R_y)} 
      #number of obs in k
      STATE[[paste0("view",view)]][[paste0("k",kkk)]][["n_k"]]=
        as.integer(TAB[as.character(kkk)])
      
      
    }
    #relabeling 
    STATE$relab[[paste0("view",view)]]=cbind(A_v=as.numeric(gsub("\\D", "", names(STATE[[paste0("view",view)]]))) ,N_k=NA )
    
    
  }
  
  #add Z matrix to state
  STATE[["Clusterall"]]=Z_mats
  #STATE$relab[[paste0("view",1)]]=cbind(A_v=as.numeric(gsub("\\D", "", names(STATE$view1))) ,N_k=NA )
  
  
  #unit specific theta/phi view 1 only
  STATE[["unit_specific_theta"]] =matrix(NA, nrow=R_y-1, ncol=N,
                                         dimnames =list(NULL, NULL  ) )
  # STATE[["unit_specificview1"]][["Phi"]] =array(NA, dim=c(D, max(R_ds), N),
  #                                         dimnames =list(colnames(X), NULL, NULL ) )
  
  
  #MATRIX FOR THETA and for BETA
  Mat_theta=matrix(NA, nrow =N, ncol = R_y-1 )
  
  for (kkk in unique(STATE$Clusterall[,1])) {
   # Mat_theta[STATE$Clusterall[,1]==kkk, 1]=STATE$view1[[paste0("k",kkk)]]$Theta[1]
    # Mat_theta[STATE$Clusterall[,1]==kkk, 2]=STATE$view1[[paste0("k",kkk)]]$Theta[2]
    
    for(nubcat in 1:R_y-1){
    Mat_theta[STATE$Clusterall[,1]==kkk, nubcat]=STATE$view1[[paste0("k",kkk)]]$Theta[nubcat]
    }
  }
  
  #Initialize beta matrix and Lp
  Mat_betas=matrix(STATE$Betas, ncol=R_y-1)
  
   
  if(P_w==0){
    Lp_state=Mat_theta*0 #LP of all zeros
  }else{
    Lp_state=W%*%Mat_betas
    
  }
  
  rm(view, Z_mats, z_temp_view, TAB, Gamma, kkk)
  
  ####################################################################
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #----------------------------#MCMC SETTINGS# ----------------------#
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ####################################################################
  #SAMPLER settings
  NBurnin = 0            #<- CHANGE HERE
  Final_size= 10^4          #<- CHANGE HERE
  thin_interval= 1          #<- CHANGE HERE
  
  
  
  
  
  
  
  
  
  
  
  
  thin_index=seq(1, length.ou=Final_size, by=thin_interval)+NBurnin  
  
  S_max =  max( thin_index)  
  progbar <- txtProgressBar(0, S_max, style = 3)
  
  index_saving=1
  
  ###variance proposal beta: CHANGE HERE here if nedded
  vbetas= NULL
  
  
  #Lsit to store results
  RESULT=list()
  
  #save general 
  RESULT[["alpha"]]=matrix(NA, ncol = Final_size, nrow = nview-1)#rep(NA, Final_size )
  RESULT[["Gamma"]]= matrix(NA, ncol = Final_size, nrow = D) #MUST BE ARRAY
  RESULT[["PHI_null"]]=array(NA, dim = c(D, max(R_ds),   Final_size ), 
                             dimnames = list(  colnames(X), NULL, NULL ))
  
  RESULT[["Betas"]]=matrix(NA, ncol =R_betas , nrow = Final_size, 
                           dimnames =list(NULL,
                                          paste0(rep(paste0("beta_r", 1:(R_y-1)), each=(P_w) ), 
                                                 rep(paste0("_w",1:P_w), P_w))
                           ))
  #views
  RESULT[["view1"]]=list()
  RESULT[["view1"]][["Z"]]= matrix(NA, ncol = Final_size, nrow = N)
  
  #view 1 thehas specific
  RESULT[["view1"]][["Theta_us"]]=array(NA, dim = c(R_y-1, N,   Final_size ), 
                                        dimnames = list(  NULL, NULL, NULL ))
  #STATE[["unit_specific"]][["Theta"]])=matrix(NA, nrow=R_y, ncol =  )
  #vire2
  RESULT[["view2"]]=list()
  RESULT[["view2"]][["Z"]]= matrix(NA, ncol = Final_size, nrow = N)
  
  
  ###SAVE predictive distribution
  if(!is.null(File_predictiveDS)){
    
  RESULT[["predictive"]]=matrix(NA, nrow = Final_size, ncol = nrow(Pred_dataset_X),
                                dimnames =list(NULL, rownames(Pred_dataset_X)) )
  
  }
  ##---------------------------------------------------------------------------
  tt=Sys.time()
  for (s in 1:S_max) {
    setTxtProgressBar(progbar, s)
    
    
    #################################
    ###  View 0 -- GIVEN gammas
    #################################
    
    #UPDATE Phi 0 for all d
    Phi0s=t( apply(Hyperparameters$clusteringvariable+S0_all_dmat*(STATE$Gamma==0),1, 
                   SampleDir)  )
    
    #PROVA[,,s]=Phi0s
    STATE$Phi0=Phi0s  
    
    
    
    #mean(TIME)
    #colMeans(t (apply(PROVA, 3, function(x) x[11,] ))[500:S_max,],na.rm = TRUE) 
    
    #reshape2::melt(t (apply(PROVA, 3, function(x) x[4,])))
    #ggplot2::ggplot(reshape2::melt(t (apply(PROVA, 3, function(x) x[4,]))),
    #               ggplot2::aes(x=value))+ggplot2::geom_density()+
    #  ggplot2::facet_wrap(Var2~.)
    
    
    
    #################################
    ###  View v-th -- GIVEN gammas
    #################################
    
    
    #################################  
    ##MCMC DIRICHLET Relevant view (1)
    # Neal 8th algorithm   
    #--------------------------------
    view=1
    Gamma_v=STATE$Gamma==view#view
    ################################  
    #UPDATE Z_i v-th Relevant view  
    Updated_Zi_REL_list=Update_Relevant_Zi(Y=Y,
                                           RY = R_y,
                                           X_is_ds = X[,Gamma_v,drop=FALSE] ,
                                           n_ =N ,
                                           M = 7, #alg 8 parameter
                                           
                                           alpha_v =STATE$Alpha[view],
                                           
                                           LINPREDBetas = Lp_state,
                                           
                                           A_v_ = unname(STATE$relab[[paste0("view",view)]][,"A_v"]),
                                           Internal_state_log =Compute_Internal_stateview1(Oldstate = STATE$view1,
                                                                                           Gamma_ind = Gamma_v) ,
                                           Z_old =STATE$Clusterall[,view] , 
                                           G0Y_sample =function(x){G_0theta(hyperparameters = Hyperparameters$theta, RY = R_y)} ,
                                           
                                           X_Prior_pred_prod_philog_i =rowSums(LogPrior_pred_phi_i[,Gamma_v, drop=FALSE]), 
                                           G0X_sample_log = function(){log(t(apply(Hyperparameters$clusteringvariable[Gamma_v, ,drop=FALSE ],1, SampleDir)) )}
    )
    
    ##Change STATE view relevant: Update Z_^v's AND save relabeled A_v_Nk 
    #Z_i^(view)  
    STATE$Clusterall[, view]=Z_k_v=Updated_Zi_REL_list$Z_updated
    
    #Active components
    STATE$relab[[paste0("view",view)]]=A_v_Nk=cbind(A_v=Updated_Zi_REL_list$Av, 
                                                    N_k=Updated_Zi_REL_list$N_k)
    # A_v_Nk=cbind(A_v=Updated_Zi_REL_list$Av, N_k=Updated_Zi_REL_list$N_k)
    
    # N_k= Updated_Zi_irr_list$N_k
    K_v=Updated_Zi_REL_list$K_v
    
    Parameters_postZ=Updated_Zi_REL_list$Param_inviewV_list
    
    
    ################################  
    #UPDATE PAREMETERS PHI^(v) and THETA -- only for k in active components
    
    S_d_comp_givenGAMMA=S_d_comp#[, , Gamma==1]
    S_d_comp_givenGAMMA[, , !Gamma_v]=S_d_comp_givenGAMMA[, , !Gamma_v]*0
    
    STATE[[paste0("view",view)]]=list()
    for (ik in  1:K_v ) {
      # Z_st_eq_to_k= Z
      #PHI: update PHI^L_K for all d 
      Sv_k_all_dmat=t(apply(S_d_comp_givenGAMMA, 3,  
                            function(x#, N=A_v_Nk[ik, "N_k"]
                            ){
                              colSums(x[Z_k_v==A_v_Nk[ik, "A_v"], ,drop=FALSE])}  ))
      
      Phi_v_k=t( apply(Hyperparameters$clusteringvariable+Sv_k_all_dmat,1, SampleDir)  ) 
      
      #print(Phi_v_k)
      #THETA:
      # ttt=Sys.time()
      
      theta_v_k=   Update_theta(thetaold=Parameters_postZ[[paste0("k",A_v_Nk[ik, "A_v"])]][["Theta"]], 
                                proposal_func=proposal_fortheta, 
                                Delta=0.5,
                                Hyperp=Hyperparameters$theta, # list "mu" "scale" "df", 
                                Y_=Y[Z_k_v==A_v_Nk[ik, "A_v"]], #Y st z_i kth group
                                MATRIX_wBetas=Lp_state[Z_k_v==A_v_Nk[ik, "A_v"], ,drop=FALSE])
      #print(Sys.time()-ttt)
      #save theta in a matrix 
      ##IMPORTANTE METTERE TUTTLE LE COLONNE IN ACCORDO TO RY-1
      Mat_theta[Z_k_v==A_v_Nk[ik, "A_v"], 1] =   theta_v_k$theta[1]
      Mat_theta[Z_k_v==A_v_Nk[ik, "A_v"], 2] =   theta_v_k$theta[2]
      
      
      
      #Update state
      STATE[[paste0("view",view)]][[paste0("k",A_v_Nk[ik, "A_v"])]][["Phi"]]=Phi_v_k
      STATE[[paste0("view",view)]][[paste0("k",A_v_Nk[ik, "A_v"])]][["Theta"]]=theta_v_k$theta
      STATE[[paste0("view",view)]][[paste0("k",A_v_Nk[ik, "A_v"])]][["n_k"]]=A_v_Nk[ik, "N_k"]
      
      #unit specific parameters only for view 1 
      STATE[["unit_specific_theta"]][ , Z_k_v==A_v_Nk[ik, "A_v"]]=theta_v_k$theta
      
    }#end for updating parameters phi and theta in active component 
    
    ################################  
    #UPDATE ALPHA^(v)
    #cat(K_v, STATE$Alpha[view], "\n")
    STATE$Alpha[view]=Update_alpha(K_nuberCluster=K_v, 
                                   nobs=N, 
                                   alphaold=STATE$Alpha[view], 
                                   Hyp2alpha=Hyperparameters$hyperhyperprioralphas)
    
    
    
    #################################  
    ##MCMC DIRICHLET Irrelevant views
    # Neal 2nd algorithm   
    #--------------------------------
    
    #L=2#<z-cambiare
    for (view in  seq_len(nview-2)+1) {
      #view=1
      Gamma_v= STATE$Gamma==view#view
      #print("view2 estimation")
      #--------------------------------    
      #view-th views
      #--------------------------------
      
      ################################  
      #UPDATE Z_i v-th irrelevant view  
      Updated_Zi_irr_list=Update_irrelevantview_Zi(X_is_ds=X[,Gamma_v, drop=FALSE],
                                                   n_=N,    
                                                   A_v_  =unname(STATE$relab[[paste0("view",view)]][,"A_v"]), #oldlables
                                                   alpha_v=STATE$Alpha[view],#STATE
                                                   Prior_pred_prod_philog_i=rowSums(LogPrior_pred_phi_i[,Gamma_v, drop=FALSE]),
                                                   Internal_state_log=lapply(STATE[[paste0("view",view)]], function(x) list(Phi=log(x$Phi)[Gamma_v, ,drop=FALSE ], n_k=x$n_k)) ,                                #Phi_ks_list=STATE[[paste0("view",view)]][["Phi"]],
                                                   Z_old=STATE$Clusterall[, view] ,
                                                   G0_sample_log=function() {log(t(apply(Hyperparameters$clusteringvariable[Gamma_v, ,drop=FALSE ],1, SampleDir)) )}
      ) 
      
      
      ##Change STATE view v: Update Z_^v's AND save relabeled A_v_Nk 
      #Z_i^(view)  
      STATE$Clusterall[, view]=Z_k_v=Updated_Zi_irr_list$Z_updated
      
      #Active components
      STATE$relab[[paste0("view",view)]]=cbind(A_v=Updated_Zi_irr_list$Av, 
                                               N_k=Updated_Zi_irr_list$N_k)
      A_v_Nk=cbind(A_v=Updated_Zi_irr_list$Av, N_k=Updated_Zi_irr_list$N_k)
      
      # N_k= Updated_Zi_irr_list$N_k
      K_v=Updated_Zi_irr_list$K_v
      
      ################################  
      #UPDATE PAREMETERS PHI^(v) -- only for k in active components
      
      S_d_comp_givenGAMMA=S_d_comp#[, , Gamma==1]
      S_d_comp_givenGAMMA[, , !Gamma_v]=S_d_comp_givenGAMMA[, , !Gamma_v]*0
      
      STATE[[paste0("view",view)]]=list()
      for (ik in  1:K_v ) {
        # Z_st_eq_to_k= Z
        #update PHI^L_K for all d 
        Sv_k_all_dmat=t(apply(S_d_comp_givenGAMMA, 3,  
                              function(x#, N=A_v_Nk[ik, "N_k"]
                              ){
                                colSums(x[Z_k_v==A_v_Nk[ik, "A_v"], ,drop=FALSE])}  ))
        
        Phi_v_k=t( apply(Hyperparameters$clusteringvariable+Sv_k_all_dmat,1, SampleDir)  ) 
        
        #print(Phi_v_k)
        STATE[[paste0("view",view)]][[paste0("k",A_v_Nk[ik, "A_v"])]][["Phi"]]=Phi_v_k
        STATE[[paste0("view",view)]][[paste0("k",A_v_Nk[ik, "A_v"])]][["n_k"]]=A_v_Nk[ik, "N_k"]
        
        #unit specific parameters only for view 1 METTERE UN IF
        # STATE$unit_specific$Phi[, , Z_k_v==A_v_Nk[ik, "A_v"]  ]=Phi_v_k
        
        
      }
      
      #Compute Phi/theta unit specific
      # STATE[["unit_specific"]][["Phi"]]= 
      
      
      ################################  
      #UPDATE ALPHA^(v)
      #cat(K_v, STATE$Alpha[view], "\n")
      STATE$Alpha[view]=Update_alpha(K_nuberCluster=K_v, 
                                     nobs=N, 
                                     alphaold=STATE$Alpha[view], 
                                     Hyp2alpha=Hyperparameters$hyperhyperprioralphas)
      # print(STATE$Alpha[view])
      
    }#end view for
    
    
    
    
    
    #################################
    ###  BETAS
    #################################
    if(P_w!=0){
    Betareturn=Update_beta(betaold=Mat_betas, #matrix 
                           proposal_func= proposal_forbeta , 
                           Y_=Y, #Y st z_i kth group
                           Wi_=W,V=vbetas,
                           Delta =1,
                           n.categ=R_y-1, #Ry-1
                           THETA_ = Mat_theta, #THETA MATRIX
                           Hyperp = Hyperparameters$beta
                           #Etas_old=Etasold,# etas z_i in kth group, 
                           #old_theta=
    )
    #update state
    Mat_betas=Betareturn$beta
    Lp_state=Betareturn$Lp
    
    STATE[["Betas"]]=Betareturn$betavec
    
    }
    #################################
    ###  VIEW ALLOCATION VARIABLES
    #################################
    #THIS FUNCTION ONLY MULTIVIEW
    
    
    STATE$Gamma=View_allocations_update_3view(X, 
                                              nvariables=D, 
                                              Z_1=STATE$Clusterall[,1],
                                              Z_2=STATE$Clusterall[,2],
                                              
                                              STATE_clusterspecpar1=STATE[[paste0("view",1)]],
                                              STATE_clusterspecpar2=STATE[[paste0("view",2)]],
                                              Succ_d0=S0_all_dmat,
                                              PHI_NULL=STATE$Phi0,
                                              Hyperpar_log=Hyperparameters[["hyperGamma_LOG"]])
    
    
    
    
    
    
    
    #-----------------------------------#
    # SAVE RESULTS #   +   # PREDICTIVE #   
    #-----------------------------------#
    if(s==thin_index[index_saving]){
      #global
      RESULT$Gamma[, index_saving]=STATE$Gamma
      RESULT$PHI_null[,,index_saving]=STATE$Phi0
      RESULT$Betas[index_saving, ]=STATE$Betas 
      #views:
      #view1
      RESULT$alpha[, index_saving ]=STATE$Alpha[1]
      RESULT$view1$Z[ , index_saving ]=STATE$Clusterall[, 1]
      RESULT$view1$Theta_us[ , ,   index_saving]=STATE$unit_specific_theta 
      #view2
      RESULT$alpha[, index_saving ]=STATE$Alpha[2]
      RESULT$view2$Z[ , index_saving ]=STATE$Clusterall[, 2]
      
      
      
      
      ################################
      #PREDICTIVE
      if(!is.null(File_predictiveDS)){
      Y_tildeS=Predictive_Y(Matrix_pred_x=Pred_dataset_X ,
                            npred=seq_along( Pred_dataset_X[,1]),
                            Matrix_pred_W=Pred_dataset_W, 
                            ry=R_y,
                            Gammas_=STATE$Gamma,
                            Beta_m=Mat_betas,
                            Alpha_v1=STATE$Alpha[1],
                            HY_theta=Hyperparameters$theta,
                            Log_prior_PRED_x=LogPrior_pred_phi_FOR_predictive_units,
                            state_view1=STATE$view1)
      
      RESULT[["predictive"]][index_saving,]=Y_tildeS
      
      }
      
      
      ##~~~
      index_saving=1+index_saving
      ##~~~
      
    }
    
  }#end SAMPLER  
  
  
  
  
  print(tt- Sys.time())
###################################################################################################
#clear space but retain the rusults
rm(setdiff(ls(),"RESULTS"))  
#######################

#->RESULT: this is a list     
  
   
  