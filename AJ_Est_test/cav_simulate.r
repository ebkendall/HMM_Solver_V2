library(msm)
# args = commandArgs(TRUE)
# num_iter = as.numeric(args[1])
num_iter = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

print(paste0("iteration ", num_iter))

set.seed(num_iter)
print(num_iter)

# p = 1 --> update continuous
# p = 2 --> update every other month
# p = 3 --> update every year
# p = 4 --> update every other year

# Set the sample size.  Note that the true cav data set has 622 subjects.
N <- 6000
# Choose the discretization of time.
dt <- 1/365

par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)

# The true values are set as the posterior means of the thinned last 15,000 steps
# from running the MCMC routine using the numerical ODE solution (seed 10).
# load('../Time_Capsuled_Code_JCGS/real_cav_analysis/Model_out/deSolve/mcmc_out_10.rda')
# chain = mcmc_out$chain[10000:25001, ]
# ind_keep = seq(1, nrow(chain), by=10)
# chain = chain[ind_keep, ]
# 
# trueValues = colMeans(chain)
# # The true values corresponding to the slope coefficient on time are scaled by
# # 3 in order to magnify the effect of the different approaches to modelling
# # continuous time HMMs.
# 
# # trueValues[6:10] = 3 * trueValues[6:10]
# trueValues[7] = -0.5
# trueValues[par_index$pi_logit] = c(0, 0)
# 
# save(trueValues, file = 'DataOut/trueValues.rda')
load('DataOut/trueValues.rda')


betaMat <- matrix(trueValues[par_index$beta], ncol = 3, byrow = F)

# Misclassification error matrix
# errorMat_temp = matrix(c(1, exp(trueValues[par_index$misclass][1]), 0, 0,
#                          exp(trueValues[par_index$misclass][2]), 1, exp(trueValues[par_index$misclass][3]), 0,
#                          0, exp(trueValues[par_index$misclass][4]), 1, 0,
#                          0,   0,   0, 1), ncol=4, byrow=TRUE)
# 
# errorMat = errorMat_temp / rowSums(errorMat_temp)
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
    # sex <- as.integer(runif(1,0,1) < propMale)
    sex <- sample(c(0,1), size = 1)
    
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
    
    # Sample inter-observation times from the cav data set.  Maximum of 20 years in study.
    visitTimes <- NULL
    time2 <- 0
    
    while(time2 < min( 20, timeOfDeath)){
        visitTimes <- c( visitTimes, time2)
        time2 <- time2 + sample( interObsTime, size=1) + runif(1, min = 0, max = 0.1)
    }
    
    # Get the exact observation times
    transition_times_pos = which(diff(trueState) != 0) + 1
    transition_times = years[transition_times_pos]
    transition_times_state = trueState[transition_times_pos]
    
    # # If the only transition is to death, then we don't need to do anything extra
    # if(length(transition_times_pos) == 1) {
    #     if(transition_times_state == 4) {
    #         transition_times_pos = NULL   
    #     }
    # }
    
    # If first visit time occurs after death, then subject is NOT entered into the study.
    if( !is.null(visitTimes) ){
    
        # visitTimes = c(0, transition_times)
        # state = c(trueState[1], transition_times_state)
        # n_i = length(visitTimes)
        # 
        # # Adding to get exact transition times ------------------------------------
        # if(length(transition_times_pos) > 0) {
        #     for(ttt in 1:length(transition_times_pos)) {
        #         if(transition_times_state[ttt] != 4) {
        #             visitTimes = c(visitTimes, transition_times[ttt])
        #         }
        #     }
        # }
        # visitTimes = sort(visitTimes)
        # # -------------------------------------------------------------------------
        
        # If death occured before the study ended, then record the time of death.
        if( timeOfDeath < 20 ){  visitTimes <- c( visitTimes, timeOfDeath) }
        
        n_i <- length(visitTimes)
        state <- NULL
        for(k in 1:n_i){  state <- c( state, tail( trueState[ years <= visitTimes[k] ], 1))  }
        
        ptnum <- rep(i,n_i)
        years <- visitTimes
        
        rawData <- rbind( rawData, data.frame(ptnum,years,sex,state) )
        # rawData <- rbind( rawData, data.frame(ptnum,years,state) )
        if(4 %in% state){  propDeaths_sim <- propDeaths_sim + 1  }
        NumObs_sim <- c( NumObs_sim, n_i)
    }
    
}

colnames(rawData) <- c('ptnum','years','sex','state')
# colnames(rawData) <- c('ptnum','years','state')
N <- length(unique(rawData$ptnum))
propDeaths_sim <- propDeaths_sim / N


# Add noise to the states.
# for(i in 1:nrow(rawData)){	rawData$state[i] <- sample(1:4, size=1, prob=errorMat[rawData$state[i],])  }

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Add censored rows.
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
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
cavData = hold

colnames(cavData) <- c('ptnum','years','disc_time','sex','state','obstrue')
rownames(cavData) <- NULL

save(cavData, file=paste("DataOut/cavData", num_iter, ".rda", sep=''))

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
