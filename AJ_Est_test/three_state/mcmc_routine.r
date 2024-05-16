library(mvtnorm, quietly=T);library(foreach, quietly=T);library(msm, quietly=T)
library(doParallel, quietly=T);library(deSolve, quietly=T)

# Construct the transition rate matrix
Q <- function(time,sex,betaMat){
    
    q1  = exp( c(1,time,sex) %*% betaMat[1,] )  
    q2  = exp( c(1,time,sex) %*% betaMat[2,] )  
    q3  = exp( c(1,time,sex) %*% betaMat[3,] )  
    q4  = exp( c(1,time,sex) %*% betaMat[4,] )  
    q5  = exp( c(1,time,sex) %*% betaMat[5,] )
    q6  = exp( c(1,time,sex) %*% betaMat[6,] )
    
    qmat = matrix(c( 0, q1, q2,
                     q3,  0, q4,
                     q5, q6,  0) ,nrow=3,byrow=TRUE)
    diag(qmat) = -rowSums(qmat)
    
    return(qmat)
}

model_t <- function(t,p,parms) {
    betaMat <- matrix(parms$b, ncol = 3, byrow = F)
    
    qmat = Q(t, parms$x_ik, betaMat)
    pmat = matrix(c(  p[1],  p[2],  p[3],
                      p[4],  p[5],  p[6],
                      p[7],  p[8],  p[9]),
                  nrow = 3, byrow = T)
    
    # Vectorizing the matrix multiplication row-wise
    dP = c(t(pmat %*% qmat))
    
    return(list(dP))
}

# Evaluating the log posterior
fn_log_post <- function(pars, prior_par, par_index, x, y, t, id, exact_time) {

    # Initial state probabilities
    init = c(1/3,1/3,1/3)

    resp_fnc = diag(3)

    beta <- matrix(pars[par_index$beta], ncol = 3)
    p_ic <- c( p1=1, p2=0, p3=0,
               p4=0, p5=1, p6=0,
               p7=0, p8=0, p9=1) # initial condition for deSolve

    # Parallelized computation of the log-likelihood
    log_total_val = foreach(i=unique(id), .combine='+', 
                            .export = c("model_t", "Q"), 
                            .packages = c("deSolve")) %dopar% {

        val = 1
        y_i = y[id == i]                # the observed state
        x_i = x[id == i,"sex",drop = F] # only the sex covariate
        t_i = t[id == i]                # continuous time

        f_i = init %*% diag(resp_fnc[, y_i[1]])
        log_norm = 0

        for(k in 2:length(t_i)) {

            out <- deSolve::ode(p_ic, times = t_i[(k-1):k], func = model_t,
                                parms = list(b=beta, x_ik = x_i[k,]))
            P <- matrix(c( out[2,"p1"],  out[2,"p2"],  out[2,"p3"],  
                           out[2,"p4"],  out[2,"p5"],  out[2,"p6"],  
                           out[2,"p7"],  out[2,"p8"],  out[2,"p9"]),
                        nrow = 3, byrow = T)

            # exact_time = T -> we need to account for exact transition time
            if(exact_time) {
                if(y_i[k] != y_i[k-1]) {
                    # Transition occurred and we have the exact transition time
                    val = f_i %*% P %*% Q(t_i[k], x_i[k,], beta) %*% diag(resp_fnc[, y_i[k]])
                } else {
                    val = f_i %*% P %*% diag(resp_fnc[, y_i[k]])
                }
            } else { 
                val = f_i %*% P %*% diag(resp_fnc[, y_i[k]])
            }

            norm_val = sqrt(sum(val^2))
            f_i = val / norm_val
            log_norm = log_norm + log(norm_val)
        }

        return(log(sum(f_i)) + log_norm)
    }

    mean = prior_par$prior_mean
    sd = diag(prior_par$prior_sd)
    log_prior_dens = dmvnorm( x=pars, mean=mean, sigma=sd, log=T)

    return(log_total_val + log_prior_dens)

}

# -----------------------------------------------------------------------------
# The mcmc routine for samping the parameters
# -----------------------------------------------------------------------------
mcmc_routine = function( y, x, t, id, init_par, prior_par, par_index,
                         steps, burnin, n_cores, exact_time){

  cl <- makeCluster(n_cores, outfile="")
  registerDoParallel(cl)

  pars = init_par
  n = length(y)
  n_par = length(pars)
  chain = matrix( 0, steps, n_par)

  group = list(c(par_index$beta[1:6]), 
               c(par_index$beta[7:12]),
               c(par_index$beta[13:18]))
  n_group = length(group)

  pcov = list();	for(j in 1:n_group)  pcov[[j]] = diag(length(group[[j]]))
  pscale = rep( .0001, n_group)
  accept = rep( 0, n_group)

  # Evaluate the log posterior of the initial parameters
  log_post_prev = fn_log_post( pars, prior_par, par_index, x, y, t, id, exact_time)

  if(!is.finite(log_post_prev)){
    print("Infinite log-posterior; choose better initial parameters")
    break
  }

  # Begin the MCMC algorithm --------------------------------------------------
  chain[1,] = pars
  for(ttt in 2:steps){
    for(j in 1:n_group){

      # Propose an update
      ind_j = group[[j]]
      proposal = pars
      proposal[ind_j] = rmvnorm( n=1, mean=pars[ind_j],sigma=pcov[[j]]*pscale[j])

      # Compute the log density for the proposal
      log_post = fn_log_post(proposal, prior_par, par_index, x, y, t, id, exact_time)

      # Only propose valid parameters during the burnin period
      if(ttt < burnin){
        while(!is.finite(log_post)){
          print('bad proposal')
          proposal = pars
          proposal[ind_j] = rmvnorm( n=1, mean=pars[ind_j],
                                     sigma=pcov[[j]]*pscale[j])
          log_post = fn_log_post(proposal, prior_par, par_index, x, y, t, id, exact_time)
        }
      }

      # Evaluate the Metropolis-Hastings ratio
      if( log_post - log_post_prev > log(runif(1,0,1)) ){
        log_post_prev = log_post
        pars[ind_j] = proposal[ind_j]
        accept[j] = accept[j] +1
      }
      chain[ttt,ind_j] = pars[ind_j]

      # Proposal tuning scheme ------------------------------------------------
      if(ttt < burnin){
        # During the burnin period, update the proposal covariance in each step
        # to capture the relationships within the parameters vectors for each
        # transition.  This helps with mixing.
        if(ttt == 100)  pscale[j] = 1

        if(100 <= ttt & ttt <= 2000){
          temp_chain = chain[1:ttt,ind_j]
          pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])

        } else if(2000 < ttt){
          temp_chain = chain[(ttt-2000):ttt,ind_j]
          pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])
        }
        if( sum( is.na(pcov[[j]]) ) > 0)  pcov[[j]] = diag( length(ind_j) )

        # Tune the proposal covariance for each transition to achieve
        # reasonable acceptance ratios.
        if(ttt %% 30 == 0){
          if(ttt %% 480 == 0){
            accept[j] = 0

          } else if( accept[j] / (ttt %% 480) < .4 ){ 
            pscale[j] = (.75^2)*pscale[j]

          } else if( accept[j] / (ttt %% 480) > .5 ){ 
            pscale[j] = (1.25^2)*pscale[j]
          }
        }
      }
      # -----------------------------------------------------------------------
    }
    # Restart the acceptance ratio at burnin.
    if(ttt == burnin)  accept = rep( 0, n_group)

    if(ttt%%1==0)  cat('--->',ttt,'\n')
  }
  # ---------------------------------------------------------------------------

  stopCluster(cl)
  print(accept/(steps-burnin))
  return(list( chain=chain[burnin:steps,], accept=accept/(steps-burnin),
               pscale=pscale))
}
# -----------------------------------------------------------------------------
