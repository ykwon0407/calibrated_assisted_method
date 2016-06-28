detach()
rm(list=ls())
gc()

library(lme4)
setwd("C:\\Users\\yc\\Dropbox\\Calibration-assisted inverse probability weighting method in clustered data\\simulation\\realdata")
source("realdata_function_code.R")
data.real = read.table("D:\\YC\\Research\\simul_2nd.txt")
head(data.real)
attach(data.real)
SEXr1 = SEXr1-1 # 2 stands for female (more data)
on.exit(detach(data.real))
dim(data.real)


sum( data.real$delta )

no.miss.cluster = unique(id[est.pi == 1])

min( table(id) )
max( table(id) )
length( table(id) )
sum( table(id) )

grp.index = c(0,cumsum(table(id)))
data.com = NULL
for( i in 1: length(unique(id)) ){
	i.index = (grp.index[i]+1):(grp.index[i+1]);
	if( i %in% no.miss.cluster){
		data.com = rbind(data.com, data.real[i.index,])
	}
}
length(no.miss.cluster)
head(data.com)
dim(data.com)

#data
data.com1 = data.real[data.real$delta==1,]
dim(data.com1)

#initialize
X = cbind(SEXr1, BLACK, HISPANIC, RACEOTH, Agegrpd7_60s, Agegrpd7_70s, Agegrpd7_80sp, MedincHigh, MedincLow, MCAID, PRVT, PAY1OTH);
fit = lmer(Y ~ SEXr1 + BLACK+ HISPANIC +RACEOTH+ Agegrpd7_60s+ Agegrpd7_70s+ Agegrpd7_80sp+ MedincHigh+ MedincLow + MCAID + PRVT + PAY1OTH + (1|id), data=data.com1)

# Complete
old = c(fixef(fit),as.data.frame(VarCorr(fit))[,4]);
result_complete <- get_proposed_EE(X=(X[data.real$delta==1,]), Y=(Y[data.real$delta==1]), a=(id[data.real$delta==1]), WEIGHT=rep(1, 51396), old = old)
com_est <- result_complete[[1]]
com_sd <- sqrt(diag(result_complete[[2]]))

# Proposed
Y[is.na(Y)] <- 0; WEIGHT <- delta/est.pi;
old = c(fixef(fit),as.data.frame(VarCorr(fit))[,4]);
result_proposed_EE <- get_proposed_EE(X=X, Y=Y, a=id, WEIGHT=WEIGHT, old = com_est)
result_proposed_EM <- get_proposed_EM(X=X, Y=Y, a=id, WEIGHT=WEIGHT, old = com_est)
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

cbind(COEF[,1], SD[,1], apply( cbind( pnorm( (COEF/SD)[,1] ), 1-pnorm( (COEF/SD)[,1] ) ), 1, min) )
cbind(COEF[,2], SD[,2], apply( cbind( pnorm( (COEF/SD)[,2] ), 1-pnorm( (COEF/SD)[,2] ) ), 1, min) )
cbind(COEF[,3], SD[,3], apply( cbind( pnorm( (COEF/SD)[,3] ), 1-pnorm( (COEF/SD)[,3] ) ), 1, min) )

