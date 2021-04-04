source("./reacmat.R")
library(rje)
fname="./Farm_sbml.xml"

#loading sbml Farm reaction network
rn.proc(fname,F)

#Four organizations of the networ
CL <- c(1,2,3,4,5,6,7,8,11,12,13,14)
GL <- c(1,2,3,4,5,6,7,8,9,10,13,14)
GCL <-c(1,2,3,4,5,6,7,8,13,14)
TN <- 1:14

#Testing above proposed organization 
rn.linp_org(rn,rn$sp.id[CL])
rn.linp_org(rn,rn$sp.id[GL])
rn.linp_org(rn,rn$sp.id[GCL])
rn.linp_org(rn,rn$sp.id[TN])

#obtaining decompostion for each organization
CL_dcom <- rn.dcom(rn,CL)
GL_dcom <- rn.dcom(rn,GL)
GCL_dcom <- rn.dcom(rn,GCL)
TN_dcom <- rn.dcom(rn,TN)


#powerset of overproduced species for the reaction network
pset_op_TN <- powerSet(TN_dcom$opsp)

#decomposition for each combination of overproduced species
dcom_pset_TN <- unique(lapply(pset_op_TN, function(x) rn.dcom.opsp(rn,TN,x)))

# List of decomposition in text
dcom_pset_TN_sp <- lapply(dcom_pset_TN, function (x) rn.dcom2sp(rn,x))
#List of all fragile cycles (whiteout catalysts)
pset_TN_fcs <- unique(lapply(dcom_pset_TN, function(x) return(x$fcs)))


#cases in consideration
#first case, only money is overproduced
dcom_pset_TN[[2049]]
dcom_pset_TN_sp[[2049]]

#cases in consideration
#first case, only money is overproduced
dcom_pset_TN[[364]]
dcom_pset_TN_sp[[364]]


#Sable case
# obtaining feasible process for overproduced set of species Case #1
v_1 <- rn.linp_org_ov(rn,rn$sp.id,c(13))*11/23 # is multiply by 11/23 so both systems porduce same amount of money
# Needed species to trigger all reaction induced by the process v_1, i.e. T(v_1(0))
x_1 <- as.integer(round((rn$mr)%*%v_1))

# obtaining feasible process for overproduced set of species Case #2
v_2 <- rn.linp_org_ov(rn,rn$sp.id,c(2,4,6,8,11,13))
# Needed species to trigger all reaction induced by the process v_2, i.e. T(v_2(0))
x_2 <- as.integer(round((rn$mr)%*%v_2))
# initial concentration is the max so both networks can trigger all reactions
x <- pmax(x_1,x_2)
x
v_1
#Discrete dynamics for twenty steps Case1
rn.dina_step_stoi(rn,1:14,x,v_1,15)
round(v_2,2)
#Discrete dynamics for twenty steps Case2
rn.dina_step_stoi(rn,1:14,x,v_2,15)


# Cows perturbation
x[3] <- 3
x
round(v_2,2)
rn.dina_step_stoi(rn,1:14,x,v_1,15)
round(v_2,2)
rn.dina_step_stoi(rn,1:14,x,v_2,15)


#water production perturbation
v_2 <- rn.linp_org_ov(rn,rn$sp.id,c(2,4,6,8,11,13))
x_2 <- as.integer(round((rn$mr)%*%v_1))
v_1 <- rn.linp_org_ov(rn,rn$sp.id,c(13))*11/23
x_1 <- as.integer(round((rn$mr)%*%v_2))
x <- pmax(x_1,x_2)
v_1[1] <- (v_1[1]/2)
v_2[1] <- (v_2[1]/2)

x
round(v_1,2)
rn.dina_step_stoi(rn,1:14,x,v_1,23)
round(v_2,2)
rn.dina_step_stoi(rn,1:14,x,v_2,23)


#structural perturbation
v_2 <- rn.linp_org_ov(rn,rn$sp.id,c(2,4,6,8,11,13))
x_2 <- as.integer(round((rn$mr)%*%v_1))
v_1 <- rn.linp_org_ov(rn,rn$sp.id,c(13))*11/23
x_1 <- as.integer(round((rn$mr)%*%v_2))
x <- pmax(x_1,x_2)

rn_m <- rn
rn_m$sp.id[15] <- "mices"
rn_m$sp.idn <- rn_m$sp.id
rn_m$spc[15] <- 15
names(rn_m$spc) <-  rn_m$sp.id

rn_m$mr <- rbind(rn$mr,rep(0,ncol(rn_m$mr)))
rn_m$mp <- rbind(rn$mp,rep(0,ncol(rn_m$mp)))

rn_m$mr <- cbind(rn_m$mr,rep(0,nrow(rn_m$mr)))
rn_m$mp <- cbind(rn_m$mp,rep(0,nrow(rn_m$mp)))

rownames(rn_m$mr) <- rn_m$sp.id
rownames(rn_m$mp) <- rn_m$sp.id

rn_m$mr[c(11,15),18] <- 1
rn_m$mp[15,18] <- 2

v_1 <- c(v_1,10)
v_2 <- c(v_2,10*11/23)
x[15] <- 10

x
v_1
rn.dina_step_stoi(rn_m,1:15,x,v_1,20)
v_2
rn.dina_step_stoi(rn_m,1:15,x,v_2,20)


