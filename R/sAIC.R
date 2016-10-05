sAIC <- function(x, y=NULL, beta, family=c("binomial","poisson","ggm")){
	
    if( !is.matrix(x) ) stop("x must be a matrix.")
    if( mode(x)!="numeric" ) stop("x must be numeric.")
	
	family <- match.arg(family)
	n <- nrow(x)
	np <- ncol(x)
	
	if( family == "binomial" ){
		x.glm <- glm(y~x,family="binomial")
		predict.glm <- predict(x.glm)
		pr0 <- 1-1/(1+exp(predict.glm))
		pr <- 1-1/(1+exp(cbind(1,x)%*%beta))
		cc <- as.vector(beta)
		jc <- abs(cc) > 1e-5
		xj <- cbind(1,x)[ ,jc]
		Pi <- diag(as.vector(pr*(1-pr)))
		Pi0 <- diag(as.vector(pr0*(1-pr0)))
		j22 <- t(xj)%*%Pi%*%xj
		i22 <- t(xj)%*%Pi0%*%xj
		aic <- -2*sum(y*log(pr)+(1-y)*log(1-pr)) + 2*sum(diag(solve(j22)%*%i22))
	}
	if( family == "poisson" ){
		x.glm <- glm(y~x,family="poisson")
		predict.glm <- predict(x.glm)
		lambda0 <- exp(predict.glm)
		lambda <- exp(cbind(1,x)%*%beta)
		cc <- as.vector(beta)
		jc <- abs(cc) > 1e-5		
		xj <- cbind(1,x)[ ,jc]
		Lambda <- diag(as.vector(lambda))
		Lambda0 <- diag(as.vector(lambda0))
		j22 <- t(xj)%*%Lambda%*%xj
		i22 <- t(xj)%*%Lambda0%*%xj
		aic <- -2*sum(y*log(lambda)-lambda-lgamma(y+1)) + 2*sum(diag(solve(j22)%*%i22))
	}
	if( family == "ggm" ){
		invSigma_tic <- beta
		Sigma_tic <- solve(beta)
		I_mat <- J_mat <- matrix(0, 0.5*np*(np+1), 0.5*np*(np+1))
		II_mat <- JJ_mat <- matrix(0, np, np)
		Sigma_est <- var(x)
		count <- 0
		for(i in 1:np)
		{
			for(j in i:np)
			{
				for(k in 1:np)
				{
					for(l in k:np)
					{
						II_mat[l,k] <- ( Sigma_est[k, i]*Sigma_est[j, l] + Sigma_est[j, k]*Sigma_est[i, l] )
						JJ_mat[l,k] <- ( Sigma_tic[k, i]*Sigma_tic[j, l] + Sigma_tic[j, k]*Sigma_tic[i, l] )
					}
				}
				count <- count + 1
				I_mat[ ,count ] <- II_mat[upper.tri(II_mat)==FALSE]
				J_mat[ ,count ] <- JJ_mat[upper.tri(JJ_mat)==FALSE]
			}
		}
		if(sum(which( invSigma_tic[upper.tri(invSigma_tic)==FALSE] == 0)) != 0) {
			I_mat <- I_mat[-which( invSigma_tic[upper.tri(invSigma_tic)==FALSE] == 0), -which( invSigma_tic[upper.tri(invSigma_tic)==FALSE] == 0)]
			J_mat <- J_mat[-which( invSigma_tic[upper.tri(invSigma_tic)==FALSE] == 0), -which( invSigma_tic[upper.tri(invSigma_tic)==FALSE] == 0)]
		}
		aic <- -n*log(det(invSigma_tic)) + n*sum(diag(invSigma_tic%*%Sigma_est)) + 2*sum( diag(I_mat%*%solve(J_mat)) )
	}
    ans <- list( AIC=aic, call=match.call() )
    class(ans) <- "sAIC"
    return(ans)
}
