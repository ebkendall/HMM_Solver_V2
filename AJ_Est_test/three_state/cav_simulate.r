library(msm)

# args = commandArgs(TRUE)
# num_iter = as.numeric(args[1])
# exact_time = as.logical(as.numeric(args[2]))
num_iter = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
exact_time = T

print(paste0("iteration ", num_iter))

set.seed(num_iter)

# Set the sample size.  Note that the true cav data set has 622 subjects.
N <- 2000
# Choose the discretization for "instantaneous" time.
dt <- 1/1000

par_index = list( beta=1:18)



# trueValues=c(matrix(c(-3,  0.50, -0.4689827,
#                       -3,  0.50,  0.2557522,
#                       -3,  0.50, -0.1457067,
#                       -3,  0.50, -0.8164156,
#                       -3,  0.50,  0.5966361,
#                       -3,  0.50,  0.7967794), ncol = 3, byrow = T))
trueValues=c(matrix(c(-2.742924, 0.4553377, -0.03769474,
                      -3.242760, 0.5973436,  0.48973652,
                      -6.874222, 1.0878027, -0.50953974,
                      -6.919877, 1.2666980,  0.72539908,
                      -7.321201, 1.1865268, -0.77354619,
                      -7.471980, 1.3772730, -0.72377494), ncol = 3, byrow = T))

# init = c(1/3, 1/3, 1/3)
init = c(1,0,0)
 

betaMat <- matrix(trueValues[par_index$beta], ncol = 3, byrow = F)


#-------------------------------------------------------------------------------
# Collect information about the real data set.
#-------------------------------------------------------------------------------
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
    if(!(4 %in% subject$state)){
        interObsTime <- c( interObsTime, round( diff(subject$years), 3))  
    }
    
    # Determine whether the subject is male.
    propMale <- propMale + as.integer(subject$sex[1]==1)
    
    # Determine whether the subject's death was observed.
    propDeaths <- propDeaths + as.integer(4 %in% subject$state)
}
propMale <- propMale / N_cav
propDeaths <- propDeaths / N_cav
#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------
# Simulate the new data
#-------------------------------------------------------------------------------
rawData <- NULL
propDeaths_sim <- 0
NumObs_sim <- NULL
i=1
while(i <= N){
    
    print(i)
    
    sex <- sample(c(0,1), size = 1)
    
    # Sample for an initial state.
    trueState <- 1 # everyone starts in state 1
    # trueState <- sample(x = 1:3, size = 1, prob = init)
    
    # Sample the remaining states until death.
    years <- 0
    time1 <- 0
    s <- trueState
    max_time = 6
    while(time1 < max_time){
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

    if(exact_time) {
        # We observe the exact transition time ---------------------------------

        # Get the exact observation times
        transition_times_pos = which(diff(trueState) != 0) + 1
        transition_times = years[transition_times_pos]
        transition_times_state = trueState[transition_times_pos]

        visitTimes = c(0, transition_times)
        state = c(trueState[1], transition_times_state)
        
        # Add more observations per subject
        visitTimes2 <- NULL
        state2 <- NULL
        time2 <- 0
        
        while(time2 < timeOfDeath){
            visitTimes2 <- c( visitTimes2, time2)
            time2 <- time2 + sample( interObsTime, size=1) + runif(1, min = 0, max = 0.1)
        }

        if(length(visitTimes2) > 1) {
            visitTimes2 = visitTimes2[-1]
            for(k in 1:length(visitTimes2)){
                state2 <- c( state2, tail( trueState[ years <= visitTimes2[k] ], 1))
            }

            combo_visits = rbind(cbind(visitTimes, state), cbind(visitTimes2, state2))
            combo_visits = combo_visits[order(combo_visits[,1]),]

            visitTimes = combo_visits[,1]
            state = combo_visits[,2]
        }
        
        n_i = length(visitTimes)

    } else {
        # We do NOT observe the exact transition time --------------------------

        # Sample inter-observation times from the cav data set. Max of 20 years.
        visitTimes <- NULL
        time2 <- 0
        
        while(time2 < timeOfDeath){
            visitTimes <- c( visitTimes, time2)
            time2 <- time2 + sample( interObsTime, size=1) + runif(1, min = 0, max = 0.1)
        }

        n_i <- length(visitTimes)
        
        state <- NULL
        for(k in 1:n_i){  
            state <- c( state, tail( trueState[ years <= visitTimes[k] ], 1))  
        }
    }

    if(n_i > 1) {
        if(n_i > 6) {
            new_n_i = sample(4:6, 1)
            visitTimes = visitTimes[1:new_n_i]
            state = state[1:new_n_i]
            n_i = new_n_i
        }
        ptnum <- rep(i,n_i)
        years <- visitTimes
        
        rawData <- rbind( rawData, data.frame(ptnum,years,sex,state) )
        if(4 %in% state){  propDeaths_sim <- propDeaths_sim + 1  }
        NumObs_sim <- c( NumObs_sim, n_i)  
        i = i + 1
    } else {
        print(paste0(i, " not enough info"))
    }
}
#-------------------------------------------------------------------------------

colnames(rawData) <- c('ptnum','years','sex','state')
N <- length(unique(rawData$ptnum))
propDeaths_sim <- propDeaths_sim / N


p = 1

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

if(exact_time) {
    save(cavData, file=paste("DataOut/exactTime/cavData", num_iter, ".rda", sep=''))
} else {
    save(cavData, file=paste("DataOut/interTime/cavData", num_iter, ".rda", sep=''))
}

obs_trans <- function(df) {
    count_transitions = matrix(0, nrow=3, ncol=3)
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
