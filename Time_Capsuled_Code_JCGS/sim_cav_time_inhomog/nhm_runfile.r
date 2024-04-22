library(nhm)

args = commandArgs(TRUE)
ind = as.numeric(args[1]) 
set.seed(ind)
print(ind)

# Setting the initial values to the "true" values which generated the data
init_par= c(-2.31617310,  -1.28756312,  -1.10116400,  -2.52367543,  -2.10384797,
            0.27050001, -11.65470594,  -0.49306415,   0.28862090,   0.22731278,
            -0.39079609,  -0.05894252,  -0.32509646,   0.48631653,   0.99565810,
            -5.28923943,  -0.90870027,  -2.40751854,  -2.44696544,  -6.52252202,
            -6.24090500)

# ------------------------------------------------------------------------------
# Bespoke version of gompertz model, centered at 5 years: ----------------------
# ------------------------------------------------------------------------------
fourstate_illness_death<-function(t,z,x) {
    tc <- t-5
    q12<-exp(x[1])
    q14<-exp(x[2])
    q23<-exp(x[3])
    q24<-exp(x[4])
    q34<-exp(x[5])
    i12<-q12*exp(x[6]*tc + x[11]*z[1])
    i14<-q14*exp(x[7]*tc + x[12]*z[1])
    i23<-q23*exp(x[8]*tc + x[13]*z[1])
    i24<-q24*exp(x[9]*tc + x[14]*z[1])
    i34<-q34*exp(x[10]*tc + x[15]*z[1])
    q<-rbind(c(0,i12,0,i14),c(0,0,i23,i24),c(0,0,0,i34),c(0,0,0,0))
    diag(q) <- c(-i12-i14,-i23-i24,-i34,0)
    der<-array(0,c(4,4,15))
    der[1,1,1]<--i12
    der[1,2,1]<-i12
    der[1,1,2]<--i14
    der[1,4,2]<-i14
    der[2,2,3]<--i23
    der[2,3,3]<-i23
    der[2,2,4]<--i24
    der[2,4,4]<-i24
    der[3,3,5]<--i34
    der[3,4,5]<-i34
    der[,,6:10]<-tc*der[,,1:5]
    der[1,1,11]<--i12*z[1]
    der[1,2,11]<-i12*z[1]
    der[1,1,12]<--i14*z[1]
    der[1,4,12]<-i14*z[1]
    der[2,2,13]<--i23*z[1]
    der[2,3,13]<-i23*z[1]
    der[2,2,14]<--i24*z[1]
    der[2,4,14]<-i24*z[1]
    der[3,3,15]<--i34*z[1]
    der[3,4,15]<-i34*z[1]
    Q<-list(q=q,qp=der)
    return(Q)
}
attr(fourstate_illness_death,"npar")<-15
attr(fourstate_illness_death,"parnames")<-c("1->2 base:","1->4 base:","2->3 base:","2->4 base:","3->4 base:","1->2 NH","1->4 NH",
                                            "2->3 NH:","2->4 NH:","3->4 NH:","1->2 Cov1:","1->4 Cov1:","2->3 Cov1:","2->4 Cov1:","3->4 Cov1:")
attr(fourstate_illness_death,"parclass")<-rep(c("Trans","Nonhom","Cov"),c(5,5,5))

# ------------------------------------------------------------------------------
# SIMULATED CAV ----------------------------------------------------------------
# ------------------------------------------------------------------------------

# loading the simulated data set
load(paste0("sim_cav_time_inhomog/DataOut/Continuous/cavData", ind, ".rda"))

cavData = as.matrix(cavData)
nhm_data = cavData[, c("state", "years", "ptnum", "sex")]
colnames(nhm_data) = c("state", "time", "id", "cov1")
nhm_data = as.data.frame(nhm_data)
nhm_data$state = as.numeric(nhm_data$state)

# Now defining the terms necessary for nhm:
trans <- rbind(c(0, 1, 0, 2),
               c(0, 0, 3, 4),
               c(0, 0, 0, 5),
               rep(0,4))

nonh <- rbind(c(0, 1, 0, 2),
              c(0, 0, 3, 4),
              c(0, 0, 0, 5),
              rep(0,4))

covm <- list("cov1" = rbind(c(0, 1, 0, 2),
                            c(0, 0, 3, 4),
                            c(0, 0, 0, 5),
                            rep(0,4)))

emat = rbind(c(0, 1, 0, 0),
             c(2, 0, 3, 0),
             c(0, 4, 0, 0),
             rep(0,4))

initp <- c(0,1,2,0)


# Create a coarsened version of the dataset for numerical stability ------------
nhm_data2 <- nhm_data
nhm_data2$time <- floor(nhm_data2$time*1000)/1000 #Coarsen time
timegap<-c(0,diff(nhm_data2$time))
nhm_data2$time[timegap < 1e-10 & nhm_data2$time!=0] <- nhm_data2$time[timegap<1e-10 & nhm_data2$time!=0]+0.001
nhm_data2$time <- floor(nhm_data2$time*1000)/1000 
# ------------------------------------------------------------------------------

gomp_model_bespoke <- model.nhm(state~time, data=nhm_data2, subject = id, 
                                trans=trans, covariates="cov1",
                                type="bespoke", intens=fourstate_illness_death,
                                emat=emat, death = T, death.states = 4,
                                firstobs="misc",initp=initp)

# Translate the parameters so they correspond to a parametrization where time 
# is centered at 5 years.
init_par3 <- init_par
init_par3[1:5] <- init_par3[1:5] + 5*init_par3[6:10]

# set upper and lower bounds
ll<-rep(-Inf,21) # Lower bounds shouldn't be an issue
ul<-rep(5,21) # Default to 5 for upper bounds

gomp_fit  <- nhm(gomp_model_bespoke, initial=init_par3, 
                 control=nhm.control(obsinfo=FALSE,splits=c(seq(0.5,12,by=0.5)),
                                     constrained=TRUE,lower_lim=ll,upper_lim=ul,
                                     nlminb_control = list(trace=1)))

save(gomp_fit, file = paste0("sim_cav_time_inhomog/Model_out/nhm/nhm_out_", ind, ".rda"))
