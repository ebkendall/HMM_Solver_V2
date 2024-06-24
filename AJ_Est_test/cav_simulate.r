floor_new <- function(t,p) {
    new_time = NULL
    if(p == 1) {
        new_time = t
    } else if (p==2) {
        monthSeq = seq(0,1,1/6)
        yearNum = floor(t)
        timeInd = max(which(monthSeq <= (t - yearNum))) # selects which month to go to
        new_time = yearNum + monthSeq[timeInd]
    } else if (p==3) {
        new_time = floor(t)
    } else if (p==4) {
        if(floor(t) %% 2 == 0) { # divisible by 2
            new_time = floor(t)
        } else {
            temp = t - 1
            new_time = floor(temp)
        }
    } else {
        print("Invalid input for p")
    }
    return(new_time)
}

censor_times <- function(t, p) {
    min_t = 0
    max_t = floor_new(max(t), p)
    new_time = c()
    if(p == 1) {
        new_time = t
    } else if (p==2) {
        new_time = seq(min_t, max_t, by = 1/6)
    } else if (p==3) {
        new_time = seq(min_t, max_t, by = 1)
    } else {
        new_time = seq(min_t, max_t, by = 2)
    }
    return(new_time)
}

library(msm)

# args = commandArgs(TRUE)
# num_iter = as.numeric(args[1]) 
# exact_time = as.logical(as.numeric(args[2]))
num_iter = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
exact_time = T

set.seed(num_iter)
print(num_iter)

# p = 1 --> update continuous
# p = 2 --> update every other month
# p = 3 --> update every year
# p = 4 --> update every other year

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
# The true values corresponding to the slope coefficient on time are scaled by
# 3 in order to magnify the effect of the different approaches to modelling 
# continuous time HMMs.
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
        
        # if(tail(trueState,1) != s) {
        #     years <- c( years, time1 + min_waitTime)
        # 
        # } else {
        #     years <- c( years, time1 + dt)
        # }
        # 
        # time1 <- time1 + dt
        
        time1 <- time1 + dt
        years <- c( years, time1)
        
        trueState <- c( trueState, s)
    }
    timeOfDeath <- tail(years,1)
    
    if(exact_time) {
        # We observe the exact transition time ---------------------------------
        
        # Get the exact observation times
        transition_times_pos = which(diff(trueState) != 0) + 1
        transition_times = years[transition_times_pos]
        transition_times_state = trueState[transition_times_pos]
        
        visitTimes = c(0, transition_times)
        state = c(trueState[1], transition_times_state)
        
        # if( timeOfDeath >= 20 ){
        #     print('removing state 4')
        #     print(rbind(visitTimes, state))
        #     visitTimes = visitTimes[-length(visitTimes)]
        #     state = state[-length(state)]
        # }
        # 
        # if(sum(abs(diff(visitTimes)) > 7) > 0) {
        #     new_vt = visitTimes[1]
        #     new_s  = state[1]
        #     diff_t = abs(diff(visitTimes))
        #     for(ttt in 2:length(visitTimes)) {
        #         if(diff_t[ttt-1] > 7) {
        #             new_vt = c(new_vt, c(mean(visitTimes[c(ttt-1, ttt)]), visitTimes[ttt]))
        #             new_s = c(new_s, c(state[ttt-1], state[ttt]))
        #         } else {
        #             new_vt = c(new_vt, visitTimes[ttt])
        #             new_s = c(new_s, state[ttt])
        #         }
        #     }
        #     
        #     visitTimes = new_vt
        #     state = new_s
        # }
        
        n_i = length(visitTimes)
        
        ptnum <- rep(i,n_i)
        years <- visitTimes
        rawData <- rbind( rawData, data.frame(ptnum,years,sex,state) )
        
        if(4 %in% state){  propDeaths_sim <- propDeaths_sim + 1  }
        NumObs_sim <- c( NumObs_sim, n_i)
        
    } else {
        # Sample inter-observation times from the cav data set.  Maximum of 20 years in study.
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
    }
    
}

colnames(rawData) <- c('ptnum','years','sex','state')
N <- length(unique(rawData$ptnum))
propDeaths_sim <- propDeaths_sim / N


# Key
# p = 1: continuous (no censor)
# p = 2: once every two months censor
# p = 3: year censor
# p = 4: once every two years censor
p = 1
disc_time <- sapply(rawData$years, floor_new, p = p)

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


if(exact_time) {
    # save(cavData, file=paste("supplement_code/DataOut/exactTime/cavData", num_iter, ".rda", sep=''))
    save(cavData, file=paste("DataOut/exactTime/cavData", num_iter, ".rda", sep=''))
} else {
    # save(cavData, file=paste("supplement_code/DataOut/interTime/cavData", num_iter, ".rda", sep=''))
    save(cavData, file=paste("DataOut/interTime/cavData", num_iter, ".rda", sep=''))
}

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
