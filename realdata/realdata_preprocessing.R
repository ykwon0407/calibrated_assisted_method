
#-----------------------------------------------------------------
# Load raw data
#-----------------------------------------------------------------

setwd("D:\\YC\\Research")
raw.data = read.csv("g4v2.csv", sep=",", header=T)
attach(raw.data)

# dim(raw.data) # 59566 279
# Y1=Cost; Y2=Cost_1stR;
# 100*sum(is.na(Y1))/length(Y1); #missing rate for the first cost # 13.69
# 100*sum(is.na(Y2))/length(Y2); #missing rate for the second cost # 13.71
# cbind(Y1,Y2) # missing for Y1 and Y2 does not occur in the same time

#-----------------------------------------------------------------
# Preprocess data
#-----------------------------------------------------------------

# Fraggrp5r2d2: 1 for different hospital
# SEXr1: 2 for female
#X=cbind(Fraggrp5r2d2, SEXr1=(SEXr1-1), BLACK, HISPANIC, RACEOTH, Agegrpd7_60s, Agegrpd7_70s, Agegrpd7_80sp, MedincHigh, MedincLow, MCAID, PRVT, PAY1OTH);
X=cbind(Fraggrp5r2d2, SEXr1, BLACK, HISPANIC, RACEOTH, Agegrpd7_60s, Agegrpd7_70s, Agegrpd7_80sp, MedincHigh, MedincLow, MCAID, PRVT, PAY1OTH);

Y=log(Cost_1stR); delta= 1-is.na(Y); #second cost
DSHOSPID = as.factor(DSHOSPID); id = as.numeric(DSHOSPID);
# cluster.size = as.vector(table(id)); #cluster size (n_i) is the number of repeated measures for each subjects.

data_raw = cbind(DSHOSPID, Y, delta, id, X); 
data = data.frame(data_raw[order(id),]);
head(data, n=10)
dim(data)
rm(list=setdiff(ls(), "data"))

detach(raw.data)
attach(data)

#-----------------------------------------------------------------
# Estimate Gamma
#-----------------------------------------------------------------

missing_index = which(delta==0) # missing.index indicates the people with missing cost
group_for_gamma = unique(id[missing_index]) # group.for.gamma indicates which groups are needed for estimating gamma
index_for_gamma = (id %in% group_for_gamma) # index.for.gamma indicates which indecies are needed for estimating gamma
data_for_gamma = data[index_for_gamma, ] # data.for.gamma indicates which data are needed for estimating gamma
length(unique(data_for_gamma[,4])) #198

p=dim(data_for_gamma)[2];
data_for_gamma = data_for_gamma[order(data_for_gamma[,4]),]; #ordering

isInclude = FALSE # include indicator variable or not
ifelse(isInclude, start_X<-5, start_X<-6)

X = data_for_gamma[,(start_X:p)];
X = matrix(unlist(X), nc=(p-(start_X-1)));
Y = data_for_gamma[,2]; id = data_for_gamma[,4]; delta = data_for_gamma[,3];
grp_index = c(0,cumsum(table(id)))
K_for_gamma = length(group_for_gamma)  # number of clusters

# initialize
fit_glm = glm(delta ~ -1+X, family = binomial(logit))

old = coef(fit_glm) #without intercept
Y.raw = Y; Y[is.na(Y)] <- 0;
while(TRUE){
	W = delta*exp(-X%*%old)
	U = matrix(0, nr=K_for_gamma, nc=length(old) ); 
	I = matrix(0, nr=length(old), nc=length(old) ); 

	for( i in 1: K_for_gamma){
		i.index = (grp_index[i]+1):(grp_index[i+1]);
		n.i = length(i.index); r.i = sum(delta[i.index])
		W.star.i = W[i.index]/sum(W[i.index]);
		W.dot.i =  W.star.i * ( matrix(rep( W.star.i%*%X[i.index,], n.i), nr=n.i, byrow=T) - X[i.index,] )
		
		#W1 = colSums(W[i.index]*X[i.index,])
		#W2 = t(X[i.index,])%*%(W[i.index]*X[i.index,])
		#I = I + (n.i-r.i)*( W1%*%t(W1)- W2*sum(W[i.index]) )/sum(W[i.index])^2;

		U[i,] = (n.i-r.i)*(t(W.star.i)%*%X[i.index,]) - (1-delta[i.index])%*%X[i.index,]
		I = I + (n.i-r.i)*( t(W.dot.i) %*% X[i.index,]);
	}

	new = old - solve(I)%*%colSums(U)
	error = sum((new-old)^2)
	old <- new
	print( error )
	if( error < 10^(-14) ){
		break;
	}	
}

# CALIBRATED PI ESTIMATION PARTS
est_pi = rep(1, dim(data)[1]);
est = rep(0, dim(data_for_gamma)[1]);
for( i in 1: K_for_gamma){
	i.index = (grp_index[i]+1):(grp_index[i+1])
	n.i = length(i.index); r.i = sum(delta[i.index])	
	Prob = 1/(1+1/(sum( (delta*exp(-X%*%old))[i.index] )/(n.i-r.i)) *(exp((X%*%old)[i.index])))
	est[i.index] = Prob
}
est_pi[index_for_gamma] = est

data_for_simul = data.frame(cbind(data,est_pi))
head(data_for_simul)

rm(list=setdiff(ls(), "data_for_simul"))
detach(data)
attach(data_for_simul)

#-----------------------------------------------------------------
# Analysis
#-----------------------------------------------------------------
library(lme4)
setwd("C:\\Users\\yc\\Dropbox\\Research\\Calibration-assisted inverse probability weighting method in clustered data\\simulation\\realdata")
source("realdata_function_code.R")

head(data_for_simul)
dim(data_for_simul)
sum( table(id) )

#data
data_obs = data_for_simul[data_for_simul$delta==1,]
dim(data_obs)

SEXr1 = SEXr1 -1;

#initialize
X = cbind(SEXr1, BLACK, HISPANIC, RACEOTH, Agegrpd7_60s, Agegrpd7_70s, Agegrpd7_80sp, MedincHigh, MedincLow, MCAID, PRVT, PAY1OTH);
fit = lmer(Y ~ SEXr1+BLACK+HISPANIC+RACEOTH+Agegrpd7_60s+Agegrpd7_70s+Agegrpd7_80sp+MedincHigh+MedincLow+MCAID+PRVT+PAY1OTH+(1|id), data=data_obs)

# Complete
old = c(fixef(fit),as.data.frame(VarCorr(fit))[,4]);
result_complete <- get_proposed_EE(X=(X[data_for_simul$delta==1,]), Y=(Y[data_for_simul$delta==1]), a=(id[data_for_simul$delta==1]), WEIGHT=rep(1, 51396), old = old)
com_est <- result_complete[[1]]
com_sd <- sqrt(diag(result_complete[[2]]))

# Proposed
Y[is.na(Y)] <- 0; WEIGHT <- delta/est_pi;
result_proposed_EE <- get_proposed_EE(X=X, Y=Y, a=id, WEIGHT=WEIGHT, old = old)
result_proposed_EM <- get_proposed_EM(X=X, Y=Y, a=id, WEIGHT=WEIGHT, old = old)
prop_EE_est <- result_proposed_EE[[1]]
prop_EE_sd <- sqrt(diag(result_proposed_EE[[2]]))
prop_EM_est <- result_proposed_EM[[1]]
prop_EM_sd <- sqrt(diag(result_proposed_EM[[2]]))

COEF = cbind( com_est, prop_EE_est, prop_EM_est ) #COEF
SD = cbind( com_sd, prop_EE_sd, prop_EM_sd ) #sd

colnames(COEF) <- c("complete case","proposed")
rownames(COEF) <- c("Intercept",colnames(X), "var.comp.random","var.error")
colnames(SD) <- c("complete case","EE","EM")
rownames(SD) <- rownames(COEF)

result_table<-cbind( cbind(COEF[,1], SD[,1], apply( cbind( pnorm( (COEF/SD)[,1] ), 1-pnorm( (COEF/SD)[,1] ) ), 1, min) ),
cbind(COEF[,2], SD[,2], apply( cbind( pnorm( (COEF/SD)[,2] ), 1-pnorm( (COEF/SD)[,2] ) ), 1, min) ),
cbind(COEF[,3], SD[,3], apply( cbind( pnorm( (COEF/SD)[,3] ), 1-pnorm( (COEF/SD)[,3] ) ), 1, min) ) )

a=pnorm( (COEF/SD)[,1:3] ); b=1-a;
p_value=a*(a<b)+b*(a>b)

result_table<-cbind( COEF[,1:3], SD[,1:3],  p_value)
colnames(result_table)<-c("COMP","EE","EM","COMP","EE","EM","COMP","EE","EM")
result_table
library(xtable)
xtable(result_table, digits=c(rep(3,4),rep(4,3),rep(3,3)) )

detach(data_for_simul)












