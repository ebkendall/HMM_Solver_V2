library(msm)

# args = commandArgs(TRUE)
# num_iter = as.numeric(args[1]) 
# case_num = as.numeric(args[2])
num_iter = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
case_num = 3 # 1, 2 or 3

# case_num = 1 --> not exact transition time nor all transitions observed
# case_num = 2 --> exact transition times
# case_num = 3 --> exact transition times and all transitions observed

set.seed(num_iter)
print(num_iter)

# Set the sample size.  Note that the true cav data set has 622 subjects.
N <- 3000
# Choose the discretization of time.
dt <- 1/365

# The true values are set as the posterior means of the thinned last 15,000 steps
# from running the MCMC routine using the numerical ODE solution (seed 10).
# load('real_cav_analysis/Model_out/deSolve/mcmc_out_10.rda')
load('mcmc_out_10.rda')
chain = mcmc_out$chain[10000:25001, ]
ind_keep = seq(1, nrow(chain), by=10)
chain = chain[ind_keep, ]

trueValues = colMeans(chain)
trueValues[6:10] = 3 * trueValues[6:10]
trueValues[7] = trueValues[8]
trueValues[8] = trueValues[6]
par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)

betaMat <- matrix(trueValues[par_index$beta], ncol = 3, byrow = F)

# Misclassification error matrix
errorMat = diag(4)

# Initial state probabilities
initProbs_temp = c( 1, exp(trueValues[par_index$pi_logit][1]), exp(trueValues[par_index$pi_logit][2]), 0)
initProbs = initProbs_temp / sum(initProbs_temp)


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Collect information about the real data set.
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

ptnum <- unique(cav$PTNUM)
N_cav <- length(ptnum)
interObsTime <- NULL
propMale <- 0
propDeaths <- 0
NumObs <- rep(0,N_cav)
for(i in 1:N_cav){
    subject <- cav[cav$PTNUM==ptnum[i],,drop=FALSE]
    
    # The number of observations for each subject.
    NumObs[i] <- nrow(subject)
    
    # The times between observations.
    if(!(4 %in% subject$state)){  interObsTime <- c( interObsTime, round( diff(subject$years), 3))  }
    
    # Determine whether the subject is male.
    propMale <- propMale + as.integer(subject$sex[1]==1)
    
    # Determine whether the subject's death was observed.
    propDeaths <- propDeaths + as.integer(4 %in% subject$state)
}
propMale <- propMale / N_cav
propDeaths <- propDeaths / N_cav



#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Fill in the states using the transition rate matrix and error matrix estimated on the real data set.
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

Q <- function(time,sex,betaMat){
    
    q1  = exp( c(1,time,sex) %*% betaMat[1,] )  # Transition from state 1 to state 2.
    q2  = exp( c(1,time,sex) %*% betaMat[2,] )  # Transition from state 1 to death.
    q3  = exp( c(1,time,sex) %*% betaMat[3,] )  # Transition from state 2 to state 3.
    q4  = exp( c(1,time,sex) %*% betaMat[4,] )  # Transition from state 2 to death.
    q5  = exp( c(1,time,sex) %*% betaMat[5,] )  # Transition from state 3 to death.
    
    qmat = matrix(c( 0,q1, 0,q2,
                     0, 0,q3,q4,
                     0, 0, 0,q5,
                     0, 0, 0, 0),nrow=4,byrow=TRUE)
    diag(qmat) = -rowSums(qmat)
    
    return(qmat)
}


rawData <- NULL
propDeaths_sim <- 0
NumObs_sim <- NULL
for(i in 1:N){
    
    print(i)
    
    # Sample the gender, as proportional to the cav data set.
    sex <- as.integer(runif(1,0,1) < propMale)
    
    # Sample for an initial state.
    trueState <- sample(1:4, size=1, prob=initProbs)
    
    # Sample the remaining states until death.
    years <- 0
    time1 <- 0
    s <- trueState
    while(s < 4){
        # Infinitesimal transition rates.
        qmat <- Q(time1,sex,betaMat)
        
        # Possible next states.
        moveToStates <- which(qmat[s,] > 0)
        
        # Sample the wait times before transition to each of the next possible states.
        waitTimes <- rexp( n=length(moveToStates), rate= qmat[s,moveToStates])
        
        # If any of the wait times are smaller than dt, then transition to the state with the minimum wait time.
        min_waitTime <- min(waitTimes)
        if(min_waitTime < dt){  s <- moveToStates[ which(waitTimes == min_waitTime) ]  }
        
        time1 <- time1 + dt
        years <- c( years, time1)
        
        trueState <- c( trueState, s)
    }
    
    timeOfDeath <- tail(years,1)
    
    if(case_num == 1) {
        
        visitTimes <- NULL
        time2 <- 0
        
        while(time2 < min( 20, timeOfDeath)){
            
            visitTimes <- c( visitTimes, time2)
            time2 <- time2 + sample( interObsTime, size=1)
        }
        
        # If first visit time occurs after death, then subject is NOT entered into the study.
        if( !is.null(visitTimes) ){
            
            # If death occured before the study ended, then record the time of death.
            if( timeOfDeath < 20 ){  visitTimes <- c( visitTimes, timeOfDeath) }
            
            n_i <- length(visitTimes)
            state <- NULL
            for(k in 1:n_i){  state <- c( state, tail( trueState[ years <= visitTimes[k] ], 1))  }
            
            ptnum <- rep(i,n_i)
            years <- visitTimes
            rawData <- rbind( rawData, data.frame(ptnum,years,sex,state) )
            
            if(4 %in% state){  propDeaths_sim <- propDeaths_sim + 1  }
            NumObs_sim <- c( NumObs_sim, n_i)
        }  
    } else if(case_num == 2) {
        
        # Using a similar approach to case 1, we can use the inter-observation 
        # times to define which transitions we observe.
        visitTimes <- NULL
        time2 <- 0
        while(time2 < timeOfDeath){
            visitTimes <- c( visitTimes, time2)
            time2 <- time2 + sample( interObsTime, size=1)
        }
        visitTimes <- c( visitTimes, timeOfDeath) 
        n_i <- length(visitTimes)
        state <- NULL
        for(k in 1:n_i){  state <- c( state, tail( trueState[ years <= visitTimes[k] ], 1))  }
        
        # Get the exact transition times ---------------------------------------
        # if length(unique(state)) == 2, then observe initial state & state 4
        if(length(unique(state)) > 2) {
            # Which state transitions do we observe?
            obs_unique_states = unique(state)
            
            # What are the exact transition times?
            transition_times_pos = which(diff(trueState) != 0) + 1
            transition_times = c(0, years[transition_times_pos])
            transition_times_state = c(trueState[1], trueState[transition_times_pos])
            
            # Exact transition times for the observed transitions
            obs_exact_trans_times = transition_times[transition_times_state %in% obs_unique_states]
            # First time point is always 0
            for(jjj in 2:length(obs_unique_states)) {
                trans_time_jjj = obs_exact_trans_times[jjj]
                obs_state_min = min(which(state == obs_unique_states[jjj]))
                
                # Double check states match up
                if(state[obs_state_min] != obs_unique_states[jjj]) stop("sim error mismatch")
                
                visitTimes[obs_state_min] = trans_time_jjj
            }
        }
        
        if(sum(diff(visitTimes) < 0) > 0) stop("something went wrong with visitTimes")
        
        n_i = length(visitTimes)
        ptnum <- rep(i,n_i)
        years <- visitTimes
        rawData <- rbind( rawData, data.frame(ptnum,years,sex,state) )
        
        if(4 %in% state){  propDeaths_sim <- propDeaths_sim + 1  }
        NumObs_sim <- c( NumObs_sim, n_i)
        
    } else if(case_num == 3) {
        # Get ALL of the exact transition times 
        transition_times_pos = which(diff(trueState) != 0) + 1
        transition_times = years[transition_times_pos]
        transition_times_state = trueState[transition_times_pos]
        
        visitTimes = c(0, transition_times)
        state = c(trueState[1], transition_times_state)
        
        n_i = length(visitTimes)
        ptnum <- rep(i,n_i)
        years <- visitTimes
        rawData <- rbind( rawData, data.frame(ptnum,years,sex,state) )
        
        if(4 %in% state){  propDeaths_sim <- propDeaths_sim + 1  }
        NumObs_sim <- c( NumObs_sim, n_i)
    }
    
}

colnames(rawData) <- c('ptnum','years','sex','state')
N <- length(unique(rawData$ptnum))
propDeaths_sim <- propDeaths_sim / N

disc_time <- rawData$years

obstrue <- rep(0,nrow(rawData))

hold <- cbind(rawData,obstrue,disc_time)
hold <- hold[,c('ptnum','years','disc_time','sex','state','obstrue')]

tempRow <- rep(0,ncol(hold))
names(tempRow) <- c('ptnum','years','disc_time','sex','state','obstrue')

num <- 1
cavData <- NULL

cavData = hold

colnames(cavData) <- c('ptnum','years','disc_time','sex','state','obstrue')
rownames(cavData) <- NULL



save(cavData, file=paste("DataOut/cavData_case", case_num, "_it", num_iter, ".rda", sep=''))

obs_trans <- function(df) {
    count_transitions = matrix(0, nrow=4, ncol=4)
    total_trans = 0
    for(i in unique(df[,"ptnum"])){
        
        b_i_mle = c(df[df[,"ptnum"] == i, 'state'])
        
        for(t in 1:(length(b_i_mle) - 1)) {
            count_transitions[b_i_mle[t], b_i_mle[t+1]] = 
                count_transitions[b_i_mle[t], b_i_mle[t+1]] + 1
            total_trans = total_trans + 1
        }
    }
    print(count_transitions)   
}

print("sex = 0")
obs_trans(cavData[cavData$sex == 0, ])
print("sex = 1")
obs_trans(cavData[cavData$sex == 1, ])

print("proportion of state 4")
print(sum(cavData$state == 4) / length(unique(cavData$ptnum)))
