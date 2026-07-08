# This script prepare the function that will be used to compute the scoring
# rules to compare and validate point process models that is source in the 
# first test

# load packages ####
library(spatstat)
library(data.table)

set.seed(2026)

### set parameter ###
models = c('ihP','LGCP','lwppp')

# temporary
fit_ihP <- goodfit
fit_ihP <- fit1
fit_LGCP <- fitox1
fit_lwpp <- lwppp
fit <- list("fit_ihP"=fit_ihP,"fit_LGCP"=fit_LGCP,"fit_lwpp"=fit_lwpp)
#

##### functions #####

# function to simulate point patterns with these parameters

sim = function(type)
{
  if(type == 'ihP') dat = simulate(fit_ihP)$`Simulation 1`
  if(type == 'LGCP') dat = simulate(fit_LGCP)$`Simulation 1`
  if(type == 'lwpp') dat = simulate(fit_lwpp)$`Simulation 1`
  return(dat)
}

#' Auxiliary function, returns for several models multiple realizations of the estimator
#' 
#' @param mod Vector of models (list of model's type)
#' @param est {The considered estimator. Needs to be a function that takes a ppp as input and returns a vector vec, the sum over this vector is then considered to be the integral over mathcal R
#'            i.e. each entry of the vector resembles \wh T(dat,r) for a fixed r. The length of vec may not depend on the input data.}
#' @param N The number of realizations of the estimator for each model
#' 
#' @return data table with length(mod) x N columns, each one labelled paste0(mod,i) where i = 1:N. Each column has n rows where n is the output length of est

computeNEst = function(mod,est,N,silence = TRUE)
{
  # get the length of the output of est:
  l_est = length(est(sim(mod[1]))) #number of distance r tested
  
  #initialize data table
  dt = data.table()
  
  for(mm in mod){
    
    if(! silence) print(mm)
    
    # initialize matrix:
    dt_temp = matrix(0,nrow = l_est,ncol = N)
    
    #simulate and fill:
    for(i in 1:N)
    {
      #if(!silence & (i %% 10 == 0)) print(paste0(i,'/',N))
      
      dat = sim(mm)
      dt_temp[,i] = est(dat)
    }
    
    dt_temp = as.data.table(dt_temp)
    setnames(dt_temp,c(paste0(mm,1:N)))
    
    dt =c(dt,dt_temp)
    
    dt = as.data.table(dt)
  }
  return(dt)
}


# computes the CRPS from the output of computeNEst  
get_crps = function(dt,models,ret_abserrmat = FALSE)
{
  # get the integrated absolute error for all column combinations in dt:
  nn = ncol(dt)
  abs_err_matrix = matrix(0,nrow = nn,ncol = nn)
  
  # fill upper triangle matrix, thereafter symmetrize
  for(i in 1:(nn-1))
  {
    if(i %% 100 == 0) print(paste0(i,'/',nn-1))
    for(j in (i+1):nn)
    {
      abs_err_matrix[i,j] = mean(abs(dt[[i]]-dt[[j]]))
    }
  }
  
  # symmetrize:
  abs_err_matrix = abs_err_matrix + t(abs_err_matrix)
  
  
  crps = matrix(0,nrow = length(models),ncol = length(models))
  
  for(dat_mod_ind  in 1:length(models))
  {
    for(fc_mod_ind  in 1:length(models))
    {
      data_mod = models[dat_mod_ind]
      fc_mod = models[fc_mod_ind]
      
      # get relevant sections of the matrix containing the mean absolute errors
      N = nn/length(models)
      
      dt_fc_mat = abs_err_matrix[(dat_mod_ind-1) * N + 1:N, (fc_mod_ind-1) * N + 1:N]
      fc_fc_mat= abs_err_matrix[(fc_mod_ind-1) * N + 1:N,(fc_mod_ind-1) * N + 1:N]
      
      crps[dat_mod_ind,fc_mod_ind] = mean(dt_fc_mat[dt_fc_mat>0]) - 1/2 * mean(fc_fc_mat[fc_fc_mat>0])
      
    }
  }
  
  CRPS = data.table(crps)
  setnames(CRPS,models)
  
  if(ret_abserrmat) return(abs_err_matrix) else return(CRPS)
}

# computes the CRPS from the output of computeNEst with comparison from the
# real point pattern calculted estimator
get_crps_PPP = function(dt,models,est_PPP,ret_abserrmat = FALSE)
{
  # get the integrated absolute error for all column in dt against the true 
  # point pattern :
  nn = ncol(dt)
  abs_err_vector = rep(0, nn)
  
  # fill upper triangle matrix, thereafter symmetrize
  for(i in 1:(nn))
  {
    if(i %% 100 == 0) print(paste0(i,'/',nn-1))
    abs_err_vector[i] = mean(abs(dt[[i]]-est_PPP))
  }
  # get the integrated absolute error for all column combinations in dt:
  nn = ncol(dt)
  abs_err_matrix = matrix(0,nrow = nn,ncol = nn)
  
  # fill upper triangle matrix, thereafter symmetrize
  for(i in 1:(nn-1))
  {
    if(i %% 100 == 0) print(paste0(i,'/',nn-1))
    for(j in (i+1):nn)
    {
      abs_err_matrix[i,j] = mean(abs(dt[[i]]-dt[[j]]))
    }
  }
  # symmetrize:
  abs_err_matrix = abs_err_matrix + t(abs_err_matrix)
  
  
  # proper continuous ranked probability score
  crps = matrix(0,nrow = 1,ncol = length(models))
  
  for(fc_mod_ind  in 1:length(models))
    {
      fc_mod = models[fc_mod_ind]
      
      # get relevant sections of the matrix containing the mean absolute errors
      N = nn/length(models)
      
      dt_fc_mat = abs_err_vector[(fc_mod_ind-1) * N + 1:N]
      fc_fc_mat= abs_err_matrix[(fc_mod_ind-1) * N + 1:N,(fc_mod_ind-1) * N + 1:N]
      
      crps[1,fc_mod_ind] = mean(dt_fc_mat[dt_fc_mat>0]) - 1/2 * mean(fc_fc_mat[fc_fc_mat>0])
      
    }
  
  CRPS = data.table(crps)
  setnames(CRPS,models)
  row.names(CRPS) = "PPP"
  
  if(ret_abserrmat) return(abs_err_matrix) else return(CRPS)
}

get_crps_errbars = function(aem,models)
{
  n = length(models)
  N = dim(aem)[1]/n
  
  mod_2_vec = paste0('obs_',rep(models,each = n),'_fc_',rep(models,n))
  
  ret_dt = data.table()
  
  for(dat_mod_ind  in 1:length(models))
  {
    for(fc_mod_ind  in 1:length(models))
    {
      data_mod = models[dat_mod_ind]
      fc_mod = models[fc_mod_ind]
      
      # get relevant sections of the matrix containing the mean absolute errors
      dt_fc_mat = aem[(dat_mod_ind-1) * N + 1:N, (fc_mod_ind-1) * N + 1:N]
      fc_fc_mat= aem[(fc_mod_ind-1) * N + 1:N,(fc_mod_ind-1) * N + 1:N]
      
      
      crps1 = apply(X = dt_fc_mat,MARGIN = 2,FUN = mean)
      crps2 = mean(fc_fc_mat[fc_fc_mat>0])
      
      crps = crps1 - 1/2 * crps2
      
      name = paste0('obs_',models[dat_mod_ind],'_fc_',models[fc_mod_ind])
      dt_temp = data.table(crps)
      setnames(dt_temp,name)
      
      ret_dt = data.table(ret_dt,dt_temp)
      
    }
  }
  return(ret_dt) 
}


#### permutation tests ####

#' Run a permutation test of the pairwise difference between two vectors of numbers
#' @param a Vector, the scores from one method
#' @param b Vector, the scores from some other method
#' @param N Integer, the size of the permutation distribution
#' @return A list with the mean of the difference and the permutation distribution of that difference
#' @examples
#' N = 1e2
#' trend  = 1:N
#' a = trend + .01 + rnorm(N, .001)
#' b = trend - .01 + rnorm(N, .001)
#' l = permutation_test_difference(a,b)
#' q = sum(l$D <= l$d_bar) / length(l$D)
#' @author Alex
#' 
permutation_test_difference = function(a,
                                       b,
                                       pval = TRUE,
                                       N = 5e3){
  n = length(a)
  d = a - b
  d_bar = mean(d)
  D = NULL
  for(i in 1:N){
    swap = rbinom(n,1,0.5)
    w_swap = which(swap == 1)
    d_i = d
    d_i[w_swap] = -d_i[w_swap]
    D[i] = mean(d_i)
  }
  
  p_val = NULL
  if(pval)
  {
    p_val  = rank(x = c(d_bar,D))[1]/(N+1)
    
  }
  
  return(list(d_bar = d_bar, D = D,p_val = p_val))
}
