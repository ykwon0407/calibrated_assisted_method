
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


get_proposed_EE <- function(X=raw.X, Y=raw.Y, a=raw.a, WEIGHT=rep(1,length(raw.Y)), old = rep(0.3, 15)){
  len = 15
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
    const_a = old[14]/(old[15]+long_unit*old[14]); #I define const_a for convenience sake.
    short_const_a = old[14]/(old[15]+units*old[14]); 
    a.hat = trans(WEIGHT*(Y-old[1]-X%*%old[2:13]), units= units, islong=TRUE);
    zeta = trans(WEIGHT*(Y-old[1]-X%*%old[2:13]), units= units)^2 - trans((WEIGHT^2-1), units= units)*old[4]
    
    U.1 = trans( (WEIGHT*(Y-old[1]-X%*%old[2:13])-const_a*a.hat), units= units)
    U.2 = trans_mat( (X*as.vector(WEIGHT*(Y-old[1]-X%*%old[2:13])-const_a*a.hat)), units= units) - trans_mat( X*(WEIGHT-1)*a.hat, units= units)/units
    U.3 = zeta - trans(WEIGHT*(Y-old[1]-X%*%old[2:13])^2, units= units)  - units*(units-1)*old[14]
    U.4 = trans(WEIGHT*(Y-old[1]-X%*%old[2:13])^2, units= units) - short_const_a*zeta - units*old[15]
    
    U = cbind(U.1, U.2, U.3, U.4)

    I[1,1] = sum( units*(short_const_a*units-1) )
    I[1,2:13] = colSums( trans_mat(WEIGHT*X, units= units)*(units*short_const_a-1) )
    I[1,14] = sum( (trans(WEIGHT*(Y-old[1]-X%*%old[2:13]), units= units))*(-units*old[15]/((old[15]+units*old[14])^2)) )
    I[1,15] = sum( (trans(WEIGHT*(Y-old[1]-X%*%old[2:13]), units= units))*units*(short_const_a^2)/old[14] )
    
    I[2:13,1] = colSums(trans_mat(X, units= units)*(units*short_const_a-1))
    I[2:13,2:13] = t(X)%*%(const_a*trans_mat(WEIGHT*X, units= units, islong=TRUE)-(WEIGHT*X)) + t(X)%*%((WEIGHT-1)*trans_mat(WEIGHT*X, units= units, islong=TRUE)/long_unit)
    I[2:13,14] = colSums( (trans(WEIGHT*(Y-old[1]-X%*%old[2:13]), units= units))*(trans_mat(X, units= units))*(-old[15]/((old[15]+units*old[14])^2)) );
    I[2:13,15] = colSums( (trans(WEIGHT*(Y-old[1]-X%*%old[2:13]), units= units))*(trans_mat(X, units= units))*((short_const_a^2)/old[14]) );
    
    I[14,1] = 2*sum( trans(WEIGHT*(Y-old[1]-X%*%old[2:13]), units= units)*(1-units) );
    I[14,2:13] = 2*colSums( -trans(WEIGHT*(Y-old[1]-X%*%old[2:13]), units= units)*trans_mat(WEIGHT*X, units= units) + trans_mat(WEIGHT*(Y-old[1]-X*old[2:13])*X, units= units) );
    I[14,14] = sum( -units*(units-1) );
    I[14,15] = sum( -trans(WEIGHT^2-1, units= units) );
    
    I[15,1] = 2*sum( trans(WEIGHT*(Y-old[1]-X%*%old[2:13]), units= units)*(units*short_const_a-1) );
    I[15,2:13] = 2*colSums( trans(WEIGHT*(Y-old[1]-X%*%old[2:13]), units= units)*trans_mat(WEIGHT*X, units= units)*(short_const_a) - trans_mat( (WEIGHT*as.vector(Y-old[1]-X%*%old[2:13])*X), units= units) );
    I[15,14] = sum( zeta*(-old[15]/(old[15]+units*old[14])^2) );
    I[15,15] = sum( zeta*(short_const_a^2)/old[14] + short_const_a*trans((WEIGHT^2-1), units= units) - units );
    
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
    I[1,2:13] = sum(trans(WEIGHT*X, units= units)*(units*short_q-1))
    I[1,3] = sum( (trans(WEIGHT*(Y-old[1]-X*old[2]), units= units))*(-units*old[4]/((old[4]+units*old[3])^2)) )
    I[1,4] = sum( (trans(WEIGHT*(Y-old[1]-X*old[2]), units= units))*units*(short_q^2)/old[3] )

    I[2:13,1] = sum(trans(X, units= units)*(units*short_q-1))
    I[2:13,2:13] = sum( X*(q*trans(WEIGHT*X, units= units, islong=TRUE)-(WEIGHT*X)) ) + sum( trans( X*(WEIGHT-1)*trans(WEIGHT*X, units= units, islong=TRUE), units= units)/units )
    I[2:13,3] = sum( (trans(WEIGHT*(Y-old[1]-X*old[2]), units= units))*(trans(X, units= units))*(-old[4]/((old[4]+units*old[3])^2)) );
    I[2:13,4] = sum( (trans(WEIGHT*(Y-old[1]-X*old[2]), units= units))*(trans(X, units= units))*((short_q^2)/old[3]) );
    
    I[3,1] = sum( 2*(short_q^2)*trans(WEIGHT*(Y-old[1]-X*old[2]), units= units)*trans(WEIGHT*(-1), units= units) );
    I[3,2:13] = sum( 2*(short_q^2)*trans(WEIGHT*(Y-old[1]-X*old[2]), units= units)*trans(WEIGHT*(-X), units= units) );
    I[3,3] = sum( (old[4]/(old[4]+units*old[3]))^2 + 2*short_q*old[4]*zeta/(old[4]+units*old[3])^2 ) - length(units);
    I[3,4] = sum( units*short_q^2 - 2*(short_q^3)*zeta/(old[3]) - (short_q^2)*(trans(WEIGHT^2-1, units= units)) );

    I[4,1] = sum(2*WEIGHT*(Y-old[1]-X*old[2])*(-1)) + sum( (units*(short_q^2)-2*short_q)*2*(trans(WEIGHT*(Y-old[1]-X*old[2]), units= units)*trans(WEIGHT*(-1), units= units)) );
    I[4,2:13] = sum(2*WEIGHT*(Y-old[1]-X*old[2])*(-X)) + sum( (units*(short_q^2)-2*short_q)*2*(trans(WEIGHT*(Y-old[1]-X*old[2]), units= units)*trans(WEIGHT*(-X), units= units)) );
    I[4,3] = sum( units*(old[4]/(old[4]+units*old[3]))^2 + zeta*(units*2*short_q*old[4]/(old[4]+units*old[3])^2 - 2*old[4]/(old[4]+units*old[3])^2 ) );
    I[4,4] = sum( (units*short_q)^2 - trans(WEIGHT^2-1, units= units)*(units*(short_q^2)-2*short_q) + zeta*(2*units*(short_q^3)/(-old[3]) + 2*(short_q^2)/old[3] ) - units );
    
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

