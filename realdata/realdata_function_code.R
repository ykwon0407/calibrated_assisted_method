
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

trans_mat <- function(mat, units, islong=FALSE){
  result_mat = c()
  cumsum_point <- c(0, cumsum(units))
  for(i in 1:length(units)){
    this_cluster <- (cumsum_point[i]+1):cumsum_point[(i+1)]
    if(islong){
	if( length(this_cluster) == 1){
		result_mat = rbind(result_mat , matrix(mat[this_cluster,], nr=units[i], nc=12, byrow=TRUE) )
	}else{
		result_mat = rbind(result_mat , matrix(colSums(mat[this_cluster,]), nr=units[i], nc=12, byrow=TRUE) )
	}      
    }else{
	if( length(this_cluster) == 1){
		result_mat = rbind(result_mat , mat[this_cluster,])
	}else{
		result_mat = rbind(result_mat , colSums(mat[this_cluster,]))
	}      
    }
  }
	return( result_mat )
}


get_proposed_EE <- function(X=raw.X, Y=raw.Y, a=raw.a, WEIGHT=rep(1,length(raw.Y)), old = rep(0.3, 15), len = 15){
  new = rep(0, length(old))
  
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
    const_a = old[(len-1)]/(old[len]+long_unit*old[(len-1)]); #I define const_a for convenience sake.
    short_const_a = old[(len-1)]/(old[len]+units*old[(len-1)]); 
    a.hat = trans(WEIGHT*(Y-old[1]-X%*%old[2:(len-2)]), units= units, islong=TRUE);
    zeta = trans(WEIGHT*(Y-old[1]-X%*%old[2:(len-2)]), units= units)^2 - trans((WEIGHT^2-1), units= units)*old[len]
    
    U.1 = trans( (WEIGHT*(Y-old[1]-X%*%old[2:(len-2)])-const_a*a.hat), units= units)
    U.2 = trans_mat( (X*as.vector(WEIGHT*(Y-old[1]-X%*%old[2:(len-2)])-const_a*a.hat)), units= units) - trans_mat( X*(WEIGHT-1)*a.hat, units= units)/units
    U.3 = zeta - trans(WEIGHT*(Y-old[1]-X%*%old[2:(len-2)])^2, units= units)  - units*(units-1)*old[(len-1)]
    U.4 = trans(WEIGHT*(Y-old[1]-X%*%old[2:(len-2)])^2, units= units) - short_const_a*zeta - units*old[len]
    
    U = cbind(U.1, U.2, U.3, U.4)

    I[1,1] = sum( units*(short_const_a*units-1) )
    I[1,2:(len-2)] = colSums( trans_mat(WEIGHT*X, units= units)*(units*short_const_a-1) )
    I[1,(len-1)] = sum( (trans(WEIGHT*(Y-old[1]-X%*%old[2:(len-2)]), units= units))*(-units*old[len]/((old[len]+units*old[(len-1)])^2)) )
    I[1,len] = sum( (trans(WEIGHT*(Y-old[1]-X%*%old[2:(len-2)]), units= units))*units*(short_const_a^2)/old[(len-1)] )
    
    I[2:(len-2),1] = colSums(trans_mat(X, units= units)*(units*short_const_a-1))
    I[2:(len-2),2:(len-2)] = t(X)%*%(const_a*trans_mat(WEIGHT*X, units= units, islong=TRUE)-(WEIGHT*X)) + t(X)%*%((WEIGHT-1)*trans_mat(WEIGHT*X, units= units, islong=TRUE)/long_unit)
    I[2:(len-2),(len-1)] = colSums( (trans(WEIGHT*(Y-old[1]-X%*%old[2:(len-2)]), units= units))*(trans_mat(X, units= units))*(-old[len]/((old[len]+units*old[(len-1)])^2)) );
    I[2:(len-2),len] = colSums( (trans(WEIGHT*(Y-old[1]-X%*%old[2:(len-2)]), units= units))*(trans_mat(X, units= units))*((short_const_a^2)/old[(len-1)]) );
    
    I[(len-1),1] = 2*sum( trans(WEIGHT*(Y-old[1]-X%*%old[2:(len-2)]), units= units)*(1-units) );
    I[(len-1),2:(len-2)] = 2*colSums( -trans(WEIGHT*(Y-old[1]-X%*%old[2:(len-2)]), units= units)*trans_mat(WEIGHT*X, units= units) + trans_mat(WEIGHT*(Y-old[1]-X*old[2:(len-2)])*X, units= units) );
    I[(len-1),(len-1)] = sum( -units*(units-1) );
    I[(len-1),len] = sum( -trans(WEIGHT^2-1, units= units) );
    
    I[len,1] = 2*sum( trans(WEIGHT*(Y-old[1]-X%*%old[2:(len-2)]), units= units)*(units*short_const_a-1) );
    I[len,2:(len-2)] = 2*colSums( trans(WEIGHT*(Y-old[1]-X%*%old[2:(len-2)]), units= units)*trans_mat(WEIGHT*X, units= units)*(short_const_a) - trans_mat( (WEIGHT*as.vector(Y-old[1]-X%*%old[2:(len-2)])*X), units= units) );
    I[len,(len-1)] = sum( zeta*(-old[len]/(old[len]+units*old[(len-1)])^2) );
    I[len,len] = sum( zeta*(short_const_a^2)/old[(len-1)] + short_const_a*trans((WEIGHT^2-1), units= units) - units );
    
    new = old - solve(I)%*%colSums(U);
    error=sum((new-old)^2); old=new;
    print(error)
    if( error<10^-8 ){
      vari = (solve(I)%*%var(U)%*%t(solve(I)))*length(units)		
      break;
    }
  }
  
  return(list(old, vari))
}



get_proposed_EM <- function(X=raw.X, Y=raw.Y, a=raw.a, WEIGHT=rep(1,length(raw.Y)), old = rep(0.3, 15), len = 15){
  new = rep(0, length(old))
  
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
    const_a = old[(len-1)]/(old[len]+long_unit*old[(len-1)]); #I define const_a for convenience sake.
    short_const_a = old[(len-1)]/(old[len]+units*old[(len-1)]); 
    a.hat = trans(WEIGHT*(Y-old[1]-X%*%old[2:(len-2)]), units= units, islong=TRUE);
    zeta = trans(WEIGHT*(Y-old[1]-X%*%old[2:(len-2)]), units= units)^2 - trans((WEIGHT^2-1), units= units)*old[len]

    U.1 = trans( (WEIGHT*(Y-old[1]-X%*%old[2:(len-2)])-const_a*a.hat), units=units) - trans((WEIGHT-1)*a.hat, units= units)/units
    U.2 = trans( X*as.vector(WEIGHT*(Y-old[1]-X%*%old[2:(len-2)])-const_a*a.hat),units= units) - trans( X*(WEIGHT-1)*a.hat, units= units)/units
    U.3 = old[len]*short_const_a+(short_const_a^2)*zeta - old[(len-1)]
    U.4 = trans(WEIGHT*(Y-old[1]-X%*%old[2:(len-2)])^2, units= units) + units*old[len]*short_const_a - (units + 2*old[len]/old[(len-1)])*(short_const_a^2)*zeta - units*old[len]
    
    U = cbind(U.1, U.2, U.3, U.4)
    
    I[1,1] = sum(units*(short_const_a*units-1))
    I[1,2:(len-2)] = colSums(trans_mat(WEIGHT*X, units=units)*(units*short_const_a-1))
    I[1,(len-1)] = sum( (trans(WEIGHT*(Y-old[1]-X%*%old[2:(len-2)]), units= units))*(-units*old[len]/((old[len]+units*old[(len-1)])^2)) )
    I[1,len] = sum( (trans(WEIGHT*(Y-old[1]-X*old[2]), units= units))*units*(short_const_a^2)/old[(len-1)] )

    I[2:(len-2),1] = colSums(trans(X, units= units)*(units*short_const_a-1))
    I[2:(len-2),2:(len-2)] = t(X)%*%(const_a*trans_mat(WEIGHT*X, units= units, islong=TRUE)-(WEIGHT*X))  + t(X)%*%((WEIGHT-1)*trans_mat(WEIGHT*X, units= units, islong=TRUE)/long_unit)  
    I[2:(len-2),(len-1)] = colSums((trans(WEIGHT*(Y-old[1]-X%*%old[2:(len-2)]), units= units))*(trans_mat(X, units= units))*(-old[len]/((old[len]+units*old[(len-1)])^2)) );
    I[2:(len-2),len] = colSums((trans(WEIGHT*(Y-old[1]-X%*%old[2:(len-2)]), units= units))*(trans_mat(X, units= units))*((short_const_a^2)/old[(len-1)]) );
    
    I[(len-1),1] = 2*sum((short_const_a^2)*trans(WEIGHT*(Y-old[1]-X*old[2]), units= units)*trans(WEIGHT*(-1), units= units) );
    I[(len-1),2:(len-2)] = 2*colSums( (short_const_a^2)*trans(WEIGHT*(Y-old[1]-X*old[2]), units= units)*trans_mat(WEIGHT*(-X), units= units) );
    I[(len-1),(len-1)] = sum( (old[len]/(old[len]+units*old[(len-1)]))^2 + 2*short_const_a*old[len]*zeta/(old[len]+units*old[(len-1)])^2 ) - length(units);
    I[(len-1),len] = sum( units*short_const_a^2 - 2*(short_const_a^3)*zeta/(old[(len-1)]) - (short_const_a^2)*(trans(WEIGHT^2-1, units= units)) );

    I[len,1] = 2*sum(WEIGHT*(Y-old[1]-X%*%old[2:(len-2)])*(-1)) + sum( (units*(short_const_a^2)-2*short_const_a)*2*(trans(WEIGHT*(Y-old[1]-X%*%old[2:(len-2)]), units= units)*trans(WEIGHT*(-1), units= units)) );

    I[len,2:(len-2)] = 2*t(-X)%*%(WEIGHT*(Y-old[1]-X%*%old[2:(len-2)])) + colSums( (units*(short_const_a^2)-2*short_const_a)*2*(trans(WEIGHT*(Y-old[1]-X%*%old[2:(len-2)]), units= units)*trans_mat(WEIGHT*(-X), units= units)) );

    I[len,(len-1)] = sum( units*(old[len]/(old[len]+units*old[(len-1)]))^2 + zeta*(units*2*short_const_a*old[len]/(old[len]+units*old[(len-1)])^2 - 2*old[len]/(old[len]+units*old[(len-1)])^2 ) );
    I[len,len] = sum( (units*short_const_a)^2 - trans(WEIGHT^2-1, units= units)*(units*(short_const_a^2)-2*short_const_a) + zeta*(2*units*(short_const_a^3)/(-old[(len-1)]) + 2*(short_const_a^2)/old[(len-1)] ) - units );
    
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

