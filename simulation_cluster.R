library(boot);
library(lme4);
#----------------------------------------------------------------------------------------
# Default setting
#----------------------------------------------------------------------------------------

isCloglog<-FALSE;scenario="3";cluster<-1000;unit<-10;

set.seed(105);
switch(scenario, 
"1"={
isMissing <- TRUE;isVarying <- FALSE;isNormal <- TRUE;
},
"2"={
isMissing <- TRUE;isVarying <- TRUE;isNormal <- TRUE; 
},
"3"={
isMissing <- TRUE;isVarying <- TRUE;isNormal <- FALSE;   
},
{
   print('default')
})
gamma1<-1;beta0<-0.25;beta1<-0.5;
ifelse(isCloglog, gamma0<-0.4, gamma0<-1.0)
#ifelse(isNormal, sigma.a<-1, sigma.a<-(15/13))

sigma.a<-;sigma.e<-1;true<-0;M<-1000;
len<-4;res1=res2=res3=res4=matrix(0,nrow=M,ncol=5)
count.res1=count.res2=count.res3=count.res4=matrix(0,nrow=M,ncol=4)
var.res3=var.res4=matrix(0, nrow=M, ncol=len^2)
true_parameter<-c(beta0,beta1,sigma.a,sigma.e)
vari1=vari2=vari3=vari4=est1=est2=est3=est4=rep(0,4)

#----------------------------------------------------------------------------------------
# Function for estimating gamma
#----------------------------------------------------------------------------------------

calculate_nr <- function(old){
  w=delta*exp(-X*old); U=I=0; U.part=NULL;
  	for( i in 1: (cluster-full_cluster)){
      this_cluster<-(1:cluster)[-full_index][i]
      this_index<-which( long_cluster %in% this_cluster )

  		w.i = w[this_index]
  		U.part[i] = (length(this_index)-Obs[this_index])[1]*sum(w.i*(X-unObs_X)[this_index])/sum(w.i)
  		U = U + U.part[i];
  		I = I + (length(this_index)-Obs[this_index])[1]*(sum((-X[this_index])*w.i*(X-unObs_X)[this_index])*sum(w.i)-
        sum(w.i*(X-unObs_X)[this_index])*sum((-X[this_index])*w.i))/(sum(w.i))^2
  	}
  new = old - U/I; 
  return(c(new, U^2, I, U.part))
}

get_gamma <- function(){
    count=1; isBad=TRUE; temp=temp2=NULL; I=U.part=NULL;
    old = sum(X-unObs_X)/sum((X-unObs_X)*X);
    while( isBad ){
        count<-count+1;
        temp<-calculate_nr(old); 
        new<-calculate_nr(temp[1]);
        	if( new[2] < temp[2] ){
        		old<-temp[1];
        		if( new[2] < 10^(-15) ){
        			isBad<-FALSE;
        		}
        	}else{
        		Check<-TRUE;
        		while(Check){
        			print(":("); 
        			if(calculate_nr((old+temp[1])/2)[2] < calculate_nr(old)[2]){
        				Check<-FALSE;
        			}else{
        				old<-(old+temp[1])/2
        			}
        		}
        	}	
    }
    old = temp[1];
    I = -temp[3]; U.part = temp[4:(cluster-full_cluster+3)]; IF = sum(U.part)/I; Z=NULL;
    
    for( i in 1: (cluster-full_cluster)){
      this_cluster<-(1:cluster)[-full_index][i]
      this_index<-which( long_cluster %in% this_cluster )
      delta_i<-delta[this_index]; X_i<-X[this_index];
      alpha[i]<-(-log((length(this_index)-Obs[this_index])[1]/sum(delta_i*exp(-old*X_i))))
      est_pi[this_index]<-inv.logit(alpha[i]+old*X_i)
    }
    
    return(c(est_pi,old));
}

#----------------------------------------------------------------------------------------
# Function for generating database and result tables
#----------------------------------------------------------------------------------------

make_data_matrix <- function(cluster,unit){
  X = matrix(runif(cluster*unit, min=-1/2,max=1/2), nrow=unit, ncol=cluster)
  if(isNormal){
    a = matrix(rep(rnorm(cluster, mean=0, sd=sigma.a),each = unit), nrow=unit, ncol=cluster)
  }else{
    #a = matrix(rep(rt(cluster, df=15), each = unit), nrow=unit, ncol=cluster)
    a = matrix(rep((rexp(cluster, rate=1)-1), each = unit), nrow=unit, ncol=cluster)
  }

  if( isMissing ){ #True value for missing data analysis
    if( isCloglog ){
      proba = 1-exp(-exp(gamma0+0.6*a+gamma1*X)); #cloglog
    }else{
      proba = inv.logit(gamma0+0.6*a+gamma1*X); #logit probability
    }
  	
  	while(TRUE){
  		delta = matrix(rbinom(n=cluster*unit, size=1, prob = proba), nrow=unit, ncol=cluster)
  		if( (sum(colMeans(delta)==0)==0) ){
  			break;
  		}
  	}
  }else{
  	proba = matrix(1, nrow=unit, ncol=cluster)
  	delta = matrix(1, nrow=unit, ncol=cluster)
  }
	Y = matrix(beta0 + beta1*X + a + rnorm(cluster*unit, mean=0, sd=sigma.e), nrow=unit, ncol=cluster);
	return(list(X,Y,delta,a,proba));
}

trans <- function(vector, units, islong=FALSE){
  result_vec = c()
  cumsum_point <- c(0, cumsum(units))
  for(i in 1:length(units)){
    this_cluster <- (cumsum_point[i]+1):cumsum_point[(i+1)]
    if(islong){
      result_vec = c(result_vec, rep(sum(vector[this_cluster]), units[i]))
    }else{
      result_vec = c(result_vec, sum(vector[this_cluster]))
    }
  }
	return( result_vec )
}

old_trans <-function(vector, repe){
  result_value <- rep(colSums(matrix(vector, ncol = cluster)), each=repe)
  return(result_value)
}

res.cal <- function(result){
  BIAS = colMeans(result)[1:4]
  VAR=c(var(result[,1]),var(result[,2]),var(result[,3]),var(result[,4]) )
  M <- rbind(BIAS,VAR,sqrt(1000)*BIAS/sqrt(VAR))
  colnames(M) <- c('intercept','beta','D','sigma')
  return(M)
}

make_long_unit <- function(units){
  long_unit=c()
  for( i in 1:length(units)){
    long_unit = c(long_unit, rep(units[i], units[i]))   
  }
  return(long_unit)
}

make_long_cluster <- function(units){
  long_id=c()
  for( i in 1:length(units)){
    long_id = c(long_id, rep(i, units[i]))   
  }
  return(long_id)
}

#----------------------------------------------------------------------------------------
# Function for parameter
#----------------------------------------------------------------------------------------

get_proposed_EE <- function(X=raw_X, Y=raw_Y, a=raw_a, WEIGHT=rep(1,length(raw_Y)), old = rep(0.3, 4)){
  len = 4
  new = rep(0, 4)
  
  class_index= c()
  class_index[1] <- 0
  class_index[2:length(a)] <- sapply(2:length(a), function(x){if(a[(x-1)] != a[x]){return(1)}else{return(0)}})
  changeing_points = c(1, which(class_index==1), (length(a)+1))
  units = diff(changeing_points)

  long_unit<-make_long_unit(units)
  
  U=matrix(0, nr=length(units), nc=len); I=matrix(0, nr=len, nc =len);
  while(TRUE){
    const_a = old[3]/(old[4]+long_unit*old[3]); #I define const_a for convenience sake.
    short_const_a = old[3]/(old[4]+units*old[3]); 
    a.hat = trans(WEIGHT*(Y-old[1]-X*old[2]), units= units, islong=TRUE);
    zeta = trans(WEIGHT*(Y-old[1]-X*old[2]), units= units)^2 - trans((WEIGHT^2-1), units= units)*old[4]
    
    U.1 = trans( (WEIGHT*(Y-old[1]-X*old[2])-const_a*a.hat), units= units)
    U.2 = trans( X*(WEIGHT*(Y-old[1]-X*old[2])-const_a*a.hat), units= units) - trans( X*(WEIGHT-1)*a.hat, units= units)/units
    U.3 = zeta - trans(WEIGHT*(Y-old[1]-X*old[2])^2, units= units)  - units*(units-1)*old[3]
    U.4 = trans(WEIGHT*(Y-old[1]-X*old[2])^2, units= units) - short_const_a*zeta - units*old[4]
    
    U = cbind(U.1, U.2, U.3, U.4)
    
    I[1,1] = sum( units*(short_const_a*units-1) )
    I[1,2] = sum( trans(WEIGHT*X, units= units)*(units*short_const_a-1) )
    I[1,3] = sum( (trans(WEIGHT*(Y-old[1]-X*old[2]), units= units))*(-units*old[4]/((old[4]+units*old[3])^2)) )
    I[1,4] = sum( (trans(WEIGHT*(Y-old[1]-X*old[2]), units= units))*units*(short_const_a^2)/old[3] )
    
    I[2,1] = sum(trans(X, units= units)*(units*short_const_a-1))
    I[2,2] = sum( X*(const_a*trans(WEIGHT*X, units= units, islong=TRUE)-(WEIGHT*X)) ) + sum( trans( X*(WEIGHT-1)*trans(WEIGHT*X, units= units, islong=TRUE), units= units)/units )
    I[2,3] = sum( (trans(WEIGHT*(Y-old[1]-X*old[2]), units= units))*(trans(X, units= units))*(-old[4]/((old[4]+units*old[3])^2)) );
    I[2,4] = sum( (trans(WEIGHT*(Y-old[1]-X*old[2]), units= units))*(trans(X, units= units))*((short_const_a^2)/old[3]) );
    
    I[3,1] = 2*sum( trans(WEIGHT*(Y-old[1]-X*old[2]), units= units)*(1-units) );
    I[3,2] = 2*sum( -trans(WEIGHT*(Y-old[1]-X*old[2]), units= units)*trans(WEIGHT*X, units= units) + trans(WEIGHT*(Y-old[1]-X*old[2])*X, units= units));
    I[3,3] = sum( -units*(units-1) );
    I[3,4] = sum( -trans(WEIGHT^2-1, units= units) );
    
    I[4,1] = 2*sum( trans(WEIGHT*(Y-old[1]-X*old[2]), units= units)*(units*short_const_a-1) );
    I[4,2] = 2*sum( trans(WEIGHT*(Y-old[1]-X*old[2]), units= units)*trans(WEIGHT*X, units= units)*(short_const_a) - trans(WEIGHT*(Y-old[1]-X*old[2])*X, units= units) );
    I[4,3] = sum( zeta*(-old[4]/(old[4]+units*old[3])^2) );
    I[4,4] = sum( zeta*(short_const_a^2)/old[3] + short_const_a*trans((WEIGHT^2-1), units= units) - units );
    
    new = old - solve(I)%*%colSums(U);
    error=sum((new-old)^2); old=new;
    print(error)
    if( error<10^-8 ){
      vari = (solve(I)%*%var(U)%*%t(solve(I)))*cluster		
      break;
    }
  }
  
  return(list(old, vari))
}

get_proposed_EM <- function(X=raw.X, Y=raw.Y, a=raw.a, WEIGHT=rep(1,length(raw.Y)), old = rep(0.3, 4)){
  len = 4
  new = rep(0, 4)
  
  class_index= c()
  class_index[1] <- 0
  class_index[2:length(a)] <- sapply(2:length(a), function(x){if(a[(x-1)] != a[x]){return(1)}else{return(0)}})
  changeing_points = c(1, which(class_index==1), (length(a)+1))
  units = diff(changeing_points)
  
  long_unit=c()
  for( i in 1:length(units)){
    long_unit = c(long_unit, rep(units[i], units[i]))   
  }
  
  U=matrix(0, nr=length(units), nc=len); I=matrix(0, nr=len, nc =len);
  while(TRUE){
    q = old[3]/(old[4]+long_unit*old[3]); #I define q for convenience sake.
    short_q = old[3]/(old[4]+units*old[3]); 
    a.hat = trans(WEIGHT*(Y-old[1]-X*old[2]), units= units, islong=TRUE);
    zeta = trans(WEIGHT*(Y-old[1]-X*old[2]), units= units)^2 - trans((WEIGHT^2-1), units= units)*old[4]
    
    U.1 = trans( (WEIGHT*(Y-old[1]-X*old[2])-q*a.hat), units= units) - trans( (WEIGHT-1)*a.hat, units= units)/units
    U.2 = trans( X*(WEIGHT*(Y-old[1]-X*old[2])-q*a.hat), units= units) - trans( X*(WEIGHT-1)*a.hat, units= units)/units
    U.3 = old[4]*short_q+(short_q^2)*zeta - old[3]
    U.4 = trans(WEIGHT*(Y-old[1]-X*old[2])^2, units= units) + units*old[4]*short_q - (units + 2*old[4]/old[3])*(short_q^2)*zeta - units*old[4]
    
    U = cbind(U.1, U.2, U.3, U.4)
    
    I[1,1] = sum(units*(short_q*units-1))
    I[1,2] = sum(trans(WEIGHT*X, units= units)*(units*short_q-1))
    I[2,1] = sum(trans(X, units= units)*(units*short_q-1))
    I[2,2] = sum( X*(q*trans(WEIGHT*X, units= units, islong=TRUE)-(WEIGHT*X)) ) + sum( trans( X*(WEIGHT-1)*trans(WEIGHT*X, units= units, islong=TRUE), units= units)/units )
    I[1,3] = sum( (trans(WEIGHT*(Y-old[1]-X*old[2]), units= units))*(-units*old[4]/((old[4]+units*old[3])^2)) )
    I[1,4] = sum( (trans(WEIGHT*(Y-old[1]-X*old[2]), units= units))*units*(short_q^2)/old[3] )
    I[2,3] = sum( (trans(WEIGHT*(Y-old[1]-X*old[2]), units= units))*(trans(X, units= units))*(-old[4]/((old[4]+units*old[3])^2)) );
    I[2,4] = sum( (trans(WEIGHT*(Y-old[1]-X*old[2]), units= units))*(trans(X, units= units))*((short_q^2)/old[3]) );
    I[3,1] = sum( 2*(short_q^2)*trans(WEIGHT*(Y-old[1]-X*old[2]), units= units)*trans(WEIGHT*(-1), units= units) );
    I[3,2] = sum( 2*(short_q^2)*trans(WEIGHT*(Y-old[1]-X*old[2]), units= units)*trans(WEIGHT*(-X), units= units) );
    
    I[3,3] = sum( (old[4]/(old[4]+units*old[3]))^2 + 2*short_q*old[4]*zeta/(old[4]+units*old[3])^2 ) - length(units);
    I[3,4] = sum( units*short_q^2 - 2*(short_q^3)*zeta/(old[3]) - (short_q^2)*(trans(WEIGHT^2-1, units= units)) );
    I[4,1] = sum(2*WEIGHT*(Y-old[1]-X*old[2])*(-1)) + sum( (units*(short_q^2)-2*short_q)*2*(trans(WEIGHT*(Y-old[1]-X*old[2]), units= units)*trans(WEIGHT*(-1), units= units)) );
    I[4,2] = sum(2*WEIGHT*(Y-old[1]-X*old[2])*(-X)) + sum( (units*(short_q^2)-2*short_q)*2*(trans(WEIGHT*(Y-old[1]-X*old[2]), units= units)*trans(WEIGHT*(-X), units= units)) );
    I[4,3] = sum( units*(old[4]/(old[4]+units*old[3]))^2 + zeta*(units*2*short_q*old[4]/(old[4]+units*old[3])^2 - 2*old[4]/(old[4]+units*old[3])^2 ) );
    I[4,4] = sum( (units*short_q)^2 - trans(WEIGHT^2-1, units= units)*(units*(short_q^2)-2*short_q) + zeta*(2*units*(short_q^3)/(-old[3]) + 2*(short_q^2)/old[3] ) - units );
    
    if(FALSE){
      W = delta*exp(-gamma.fit*X)
      W.X = trans(W*X, 10); W.1 = trans(W*1, 10);
      del.pi = proba*(1-proba)*(X - W.X/W.1)
      U.5 = trans((WEIGHT-1)*X,1)
      U.star = cbind(U.1, U.2, U.3, U.4, U.5)
      
      I.star[1:4,1:4] = I[1:4,1:4]
      I.star[1,5] = sum( ((-WEIGHT^2)*del.pi*(Y-old[1]-X*old[2]-q*a.hat) + WEIGHT*(q*trans(WEIGHT^2*(del.pi)*(Y-old[1]-X*old[2]),1))))
      I.star[2,5] = sum( (X*(-WEIGHT^2)*del.pi*(Y-old[1]-X*old[2]-a.hat/unit) + X*WEIGHT*(trans(WEIGHT^2*(del.pi)*(Y-old[1]-X*old[2]),1)/unit))*unit*q)
      I.star[3,5] = sum(q^2*(2*trans(WEIGHT*(Y-old[1]-X*old[2]),1)*trans(-WEIGHT^2*(del.pi)*(Y-old[1]-X*old[2]),1)+2*old[4]*trans(WEIGHT^3*(del.pi),1) ))
      I.star[4,5] = sum( unit*(old[4]/(old[4]+unit*old[3]))^2 + zeta*(unit*2*q*old[4]/(old[4]+unit*old[3])^2 - 2*old[4]/(old[4]+unit*old[3])^2) ) 
      I.star[5,5] = sum(-WEIGHT^2*del.pi*X)
      new = old - solve(I.star)%*%colSums(U.star)
    }
    
    new = old - solve(I)%*%colSums(U);
    error=sum((new-old)^2); old=new;
    print(error)
    if( error<10^-8 ){
      vari = (solve(I)%*%var(U)%*%t(solve(I)))*cluster    
      break;
    }
  }
  
  return(list(old, vari))
}


#----------------------------------------------------------------------------------------
# Simulation
#----------------------------------------------------------------------------------------

for(k in 1:M){
  
  #-----------------------------------------------------------------------
  # generating data
  #-----------------------------------------------------------------------
  
  isPASS<-TRUE
  while(isPASS){
  isPASS<-FALSE
  raw = make_data_matrix(cluster,unit);
  raw_matrix<-sapply(raw,c);    
  if( isVarying ){
    units<-rbinom(cluster, unit, 0.9) #subset of raw matrix
    subset_index<-c();
    for(j in 1:cluster){
      cluster_index = ((j-1)*unit+1):((j-1)*unit+units[j])
      subset_index <- c(subset_index, cluster_index)
	if( sum(raw_matrix[cluster_index,3]) == 0 ){
		isPASS<-TRUE
		break;
	}
    }
    raw_matrix<-raw_matrix[subset_index,]
  }else{
    units<-rep(unit, cluster) #subset of raw matrix
  }
  long_cluster<-make_long_cluster(units)
  }

  #X=Y=delta=a=NULL;	
  X=raw_X=raw_matrix[,1]; Y=raw_Y=raw_matrix[,2];
  delta=raw_delta=raw_matrix[,3]; a=raw_a=raw_matrix[,4];
  proba=raw_proba=raw_matrix[,5]; raw_long_cluster=long_cluster
  hist(proba)
  
  #-----------------------------------------------------------------------
  # estimating gamma
  #-----------------------------------------------------------------------
  
  if( isMissing ){
      # raw_* stands for real data and variables without 'raw_' prefix are for gamma calculation
      if ( sum( tapply(raw_delta, long_cluster, FUN=sum) == units ) > 0 ) {
      	full_index = which( tapply(raw_delta, long_cluster, FUN=sum) == units )
      	full_cluster = length(full_index);
        del_full_cluster = which( long_cluster %in% full_index );
      	X<-raw_X[-del_full_cluster]; a<-raw_a[-del_full_cluster];
        delta<-raw_delta[-del_full_cluster]; Y<-raw_Y[-del_full_cluster];
        long_cluster<-long_cluster[-del_full_cluster]
      }else{
        full_cluster<-0;
	  full_index<-(cluster+1);
      }
    	Obs = unObs_X = item = Obs_X = rep(0,length(X)); alpha= NULL;
    	for( i in 1:(cluster-full_cluster) ){
        this_cluster<-(1:cluster)[-full_index][i]
        this_index<-which( long_cluster %in% this_cluster )
    		delta_i<-delta[this_index]; X_i<-X[this_index];
    		Obs[this_index]<-sum(delta_i); #r_i
    		unObs_X[this_index]<-rep(sum((1-delta_i)*X_i)/sum(1-delta_i), each=length(this_index));
    		Obs_X[this_index]<-rep(sum(delta_i*X_i)/sum(delta_i), each=length(this_index));
    	}
      est_pi = rep(0,length(X));
      fit_gamma = get_gamma()
      est_pi = fit_gamma[1:(length(fit_gamma)-1)]
      
      proba=c()
      if ( sum( tapply(raw_delta, raw_long_cluster, FUN=sum) == units ) > 0 ) {
      	proba=rep(1,length(raw_X));
      	proba[c(-del_full_cluster)] = as.vector(est_pi);
      }else{
      	proba = est_pi;
      }
  }
  
  X=raw_X; a=raw_a; delta=raw_delta; Y=raw_Y; WEIGHT=raw_delta/proba;
  X.obs=X[delta==1]; a.obs=a[delta==1]; delta.obs=delta[delta==1]; Y.obs=Y[delta==1];
  
  lmer_fit<-lmer(Y.obs ~ X.obs + (1|a.obs))
  init_old<-c(fixef(lmer_fit), unclass(as.data.frame(VarCorr(lmer_fit)))$vcov) 
  #-----------------------------------------------------------------------
  # full cases
  #-----------------------------------------------------------------------
  
  fit.full = get_proposed_EE(X = raw_X, Y=raw_Y, a=raw_a, WEIGHT = rep(1, length(raw_Y)), old = init_old)
  est1 = fit.full[[1]]
  vari1 = sqrt(diag(fit.full[[2]]))
  
  #-----------------------------------------------------------------------
  # complete cases
  #-----------------------------------------------------------------------
  
  fit.complete = get_proposed_EE(X = X.obs, Y=Y.obs, a=a.obs, WEIGHT = rep(1, length(Y.obs)), old = init_old)
  est2 = fit.complete[[1]]
  vari2 = sqrt(diag(fit.complete[[2]]))
  
  #-----------------------------------------------------------------------
  # proposed cases
  #-----------------------------------------------------------------------
  
  fit.proposed = get_proposed_EE(X=X, Y=Y, a=a, WEIGHT=WEIGHT, old = init_old)
  est3 = fit.proposed[[1]]
  vari3 = sqrt(diag(fit.proposed[[2]]))

  #-----------------------------------------------------------------------
  # proposed cases
  #-----------------------------------------------------------------------
  
  fit.proposed_EM = get_proposed_EM(X=X, Y=Y, a=a, WEIGHT=WEIGHT, old = init_old)
  est4 = fit.proposed_EM[[1]]
  vari4 = sqrt(diag(fit.proposed_EM[[2]]))

count.res1[k,] = as.numeric( (est1 + qnorm(0.025)*vari1-true_parameter)*(est1 - qnorm(0.025)*vari1-true_parameter) < 0 );
count.res2[k,] = as.numeric( (est2 + qnorm(0.025)*vari2-true_parameter)*(est2 - qnorm(0.025)*vari2-true_parameter) < 0 );
count.res3[k,] = as.numeric( (est3 + qnorm(0.025)*vari3-true_parameter)*(est3 - qnorm(0.025)*vari3-true_parameter) < 0 );
count.res4[k,] = as.numeric( (est4 + qnorm(0.025)*vari4-true_parameter)*(est4 - qnorm(0.025)*vari4-true_parameter) < 0 );


res1[k,] = c(est1-true_parameter, mean(delta))
res2[k,] = c(est2-true_parameter, mean(delta))
res3[k,] = c(est3-true_parameter, mean(delta))
res4[k,] = c(est4-true_parameter, mean(delta))

var.res3[k,] = as.vector(vari3)
var.res4[k,] = as.vector(vari3)
}

colMeans(count.res1)
colMeans(count.res2)
colMeans(count.res3)
colMeans(count.res4)

res.cal(res1);
res.cal(res2);
res.cal(res3);
res.cal(res4);

T = matrix(colMeans(var.res^2), nr=4)
cluster;unit;gamma1;

var(res3[,1:4])/T[1:4,1:4] #Check the ratio between simulation variance and estimates

#########################################################

library("xtable")
Bias = rbind(res.cal(res1)[1,],res.cal(res2)[1,],res.cal(res3)[1,],res.cal(res4)[1,])*1000
Variance = rbind(res.cal(res1)[2,],res.cal(res2)[2,],res.cal(res3)[2,],res.cal(res4)[2,])*1000
CP = rbind(colMeans(count.res1),colMeans(count.res2),colMeans(count.res3),colMeans(count.res4))*1000

DF<-data.frame(t(Bias),t(Variance),t(CP))
colnames(DF) =c("FULL","COMP","PROP_EE","PROP_EM","FULL","COMP","PROP_EE","PROP_EM","FULL","COMP","PROP_EE","PROP_EM") 
xtable(DF, digits=c(rep(1,9),rep(0,4)), display=c(rep("f",9),rep("f",4)), align=rep("c",13))

cluster;unit;isVarying;isNormal;isMissing;scenario;colMeans(res3)[5]
