library(mvtnorm, quietly=T);library(foreach, quietly=T);library(msm, quietly=T)
library(doParallel, quietly=T);library(deSolve, quietly=T);library(expm, quietly=T)

# Construct the transition rate matrix
Q <- function(t,x_ik,beta){
    
    betaMat = matrix(beta, ncol = 3, byrow = F) 
    q1  = exp( c(1,t,x_ik) %*% betaMat[1,] )  # Transition from state 1 to state 2.
    q2  = exp( c(1,t,x_ik) %*% betaMat[2,] )  # Transition from state 1 to death.
    q3  = exp( c(1,t,x_ik) %*% betaMat[3,] )  # Transition from state 2 to state 3.
    q4  = exp( c(1,t,x_ik) %*% betaMat[4,] )  # Transition from state 2 to death.
    q5  = exp( c(1,t,x_ik) %*% betaMat[5,] )  # Transition from state 3 to death.
    
    qmat = matrix(c( 0,q1, 0,q2,
                     0, 0,q3,q4,
                     0, 0, 0,q5,
                     0, 0, 0, 0),nrow=4,byrow=TRUE)
    diag(qmat) = -rowSums(qmat)
    
    return(qmat)
}

# Setting up the differential equations for deSolve
model_t <- function(t,p,parms) {
    
    betaMat <- matrix(parms$b, ncol = 3, byrow = F)
    
    q1  = exp( c(1,t,parms$x_ik) %*% betaMat[1,] )  # Transition from state 1 to state 2.   
    q2  = exp( c(1,t,parms$x_ik) %*% betaMat[2,] )  # Transition from state 1 to death.     
    q3  = exp( c(1,t,parms$x_ik) %*% betaMat[3,] )  # Transition from state 2 to state 3.   
    q4  = exp( c(1,t,parms$x_ik) %*% betaMat[4,] )  # Transition from state 2 to death.     
    q5  = exp( c(1,t,parms$x_ik) %*% betaMat[5,] )  # Transition from state 3 to death.     
    
    dP = rep(1,9) # this is the vector with all differential equations
    
    dP[1] = p[1]*(-q1-q2)
    dP[2] = p[1]*q1 + p[2]*(-q3-q4)
    dP[3] = p[2]*q3 - p[3]*q5
    dP[4] = p[1]*q2 + p[2]*q4 + p[3]*q5
    dP[5] = p[5]*(-q3-q4)
    dP[6] = p[5]*q3 - p[6]*q5
    dP[7] = p[5]*q4 + p[6]*q5
    dP[8] = -p[8]*q5
    dP[9] = p[8]*q5
    
    return(list(dP))
    
}

# Evaluating the log posterior
fn_log_post <- function(pars, prior_par, par_index, x, y, t, id, disc, exact_time) {
    
    # Initial state probabilities
    init_logit = c( 1, exp(pars[par_index$pi_logit][1]), exp(pars[par_index$pi_logit][2]), 0)
    init = init_logit / sum(init_logit)
    
    # Misclassification response matrix
    resp_fnc = diag(4)
    
    beta <- pars[par_index$beta]
    p_ic <- c(p1=1,p2=0,p3=0,p4=0,p5=1,p6=0,p7=0,p8=1,p9=0) # initial condition for deSolve
    
    # Parallelized computation of the log-likelihood
    log_total_val = foreach(i=unique(id), .combine='+', 
                            .export = c("model_t", "Q"), 
                            .packages = c("deSolve", "expm")) %dopar% {
                                
        val = 1; disc_t_i = NULL
        y_i = y[id == i]                # the observed state
        x_i = x[id == i,"sex",drop = F] # only the sex covariate
        t_i = t[id == i]                # continuous time
        
        # If time is discretized for matrix exponential (disc = TRUE)
        if(disc) { disc_t_i = x[id == i,"disc_time",drop = F] }
        
        f_i = init %*% diag(resp_fnc[, y_i[1]])
        
        log_norm = 0
        
        for(k in 2:length(t_i)) {
            
            if(disc) {
                P = expm((t_i[k] - t_i[k-1]) * Q(disc_t_i[k-1], x_i[k-1,], beta))
            } else {
                out <- deSolve::ode(p_ic, times = t_i[(k-1):k], func = model_t,
                                    parms = list(b=beta, x_ik = x_i[k,]))
                
                P <- matrix(c(out[2,"p1"], out[2,"p2"], out[2,"p3"], out[2,"p4"],
                              0, out[2,"p5"], out[2,"p6"], out[2,"p7"],
                              0,  0, out[2,"p8"], out[2,"p9"],
                              0,  0,  0,  1), nrow = 4, byrow = T)
            }
            
            if(exact_time) {
                if(y_i[k-1] != y_i[k]) {
                    q_temp = Q(t_i[k], x_i[k,], beta)
                    diag(q_temp) = 0
                    val = f_i %*% P %*% q_temp %*% diag(resp_fnc[, y_i[k]])
                    # val = f_i * P[y_i[k-1], y_i[k-1]] * Q(t_i[k], x_i[k,], beta)[y_i[k-1], y_i[k]]   
                } else {
                    # val = f_i * P[y_i[k-1], y_i[k]]
                    val = f_i %*% P %*% diag(resp_fnc[, y_i[k]])
                }
            } else {
                # Checking if death (state = 4) has occurred
                if(y_i[k] != 4) {
                    # Checking if the row is observed or a censored row
                    if(y_i[k] != 99) {
                        val = f_i %*% P %*% diag(resp_fnc[, y_i[k]])
                    } else { # censor row (only happens when disc == TRUE)
                        val = f_i %*% P %*% diag(rowSums(resp_fnc[, 1:3]))
                    }
                } else { # death is observed
                    if(disc){
                        val = f_i %*% P %*% Q(disc_t_i[k], x_i[k,], beta) %*% diag(resp_fnc[, y_i[k]])
                    } else{
                        val = f_i %*% P %*% Q(t_i[k], x_i[k,], beta) %*% diag(resp_fnc[, y_i[k]])
                    }
                }   
            }
            
            norm_val = sqrt(sum(val^2))
            f_i = val / norm_val
            log_norm = log_norm + log(norm_val)
        }
        
        if(sum(f_i) < 0) f_i = 0
        
        return(log(sum(f_i)) + log_norm)
    }
    
    mean = prior_par$prior_mean
    sd = diag(prior_par$prior_sd)
    log_prior_dens = dmvnorm( x=pars[c(par_index$beta, par_index$pi_logit)], mean=mean, sigma=sd, log=T)
    
    return(log_total_val + log_prior_dens)
    
}

# -----------------------------------------------------------------------------
# The mcmc routine for samping the parameters
# -----------------------------------------------------------------------------
mcmc_routine = function( y, x, t, id, init_par, prior_par, par_index,
                         steps, burnin, n_cores, disc, exact_time){
    
    cl <- makeCluster(n_cores, outfile="")
    registerDoParallel(cl)
    
    pars = init_par
    n = length(y)
    n_par = length(pars)
    chain = matrix( 0, steps, n_par)
    
    group = list(c(par_index$beta))
    n_group = length(group)
    
    pcov = list();	for(j in 1:n_group)  pcov[[j]] = diag(length(group[[j]]))
    pscale = rep( .0001, n_group)
    accept = rep( 0, n_group)
    
    # Evaluate the log posterior of the initial parameters
    log_post_prev = fn_log_post( pars, prior_par, par_index, x, y, t, id, disc, exact_time)
    
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
            log_post = fn_log_post(proposal, prior_par, par_index, x, y, t, id, disc, exact_time)
            
            # Only propose valid parameters during the burnin period
            if(ttt < burnin){
                prop_count = 0
                while(!is.finite(log_post)){
                    print('bad proposal')
                    proposal = pars
                    proposal[ind_j] = rmvnorm( n=1, mean=pars[ind_j],
                                            sigma=pcov[[j]]*pscale[j])
                    log_post = fn_log_post(proposal, prior_par, par_index, x, y, t, id, disc, exact_time)

                    prop_count = prop_count + 1

                    if(prop_count > 20) { break }
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
                if(ttt == 100)  pscale[j] = .0001
                
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
        
        if(ttt %% 100 == 0) {
            print(matrix(pars[par_index$beta], ncol = 3))
            print(pscale)
            print(pcov)
            print(accept)
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
