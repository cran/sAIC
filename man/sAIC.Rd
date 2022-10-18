\name{sAIC}
\alias{sAIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute the Akaike information criterion for the lasso in generalized linear models}
\description{This function computes the Akaike information criterion for generalized linear models estimated by the lasso. }
\usage{sAIC(x, y=NULL, beta, family=c("binomial","poisson","ggm"))} 
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A data matrix.
}
  \item{y}{
A response vector. If you select family="ggm", you should omit this argument.
}
  \item{beta}{
An estimated coefficient vector including the intercept.  If you select family="ggm", you should use an estimated precision matrix.
}
  \item{family}{
Response type (binomial, Poisson or Gaussian graphical model).
}
}
\value{
\item{AIC}{The value of AIC.}
}
\references{
%% ~put references to the literature/web site here ~
%Friedman, J. (2008).  Fast sparse regression and classification. \eqn{Technical report}, Standford University.\cr
Ninomiya, Y. and Kawano, S. (2016). 
 \emph{AIC for the Lasso in generalized linear models.}
 \emph{Electronic Journal of Statistics, 10, 2537--2560.}
 \doi{10.1214/16-EJS1179}
}
\author{
Shuichi Kawano\cr
\email{skawano@math.kyushu-u.ac.jp}
}
%\note{
%  ~~further notes~~
%}
% ~Make other sections like Warning with \section{Warning }{....} ~
%\seealso{\code{cv.spcr}}
\examples{
library(MASS)
library(glmnet)
library(glasso)

### logistic model
set.seed(3)
n <- 100; np <- 10; beta <- c(rep(0.5,3), rep(0,np-3))
Sigma <- diag( rep(1,np) )
for(i in 1:np) for(j in 1:np) Sigma[i,j] <- 0.5^(abs(i-j))
x <- mvrnorm(n, rep(0, np), Sigma)
y <- rbinom(n,1,1-1/(1+exp(x\%*\%beta)))
glmnet.object <- glmnet(x,y,family="binomial",alpha=1)
coef.glmnet <- coef(glmnet.object)
### coefficients
coef.glmnet[ ,10]
### AIC
sAIC(x=x, y=y, beta=coef.glmnet[ ,10], family="binomial")

### Poisson model
set.seed(1)
n <- 100; np <- 10; beta <- c(rep(0.5,3), rep(0,np-3))
Sigma <- diag( rep(1,np) )
for(i in 1:np) for(j in 1:np) Sigma[i,j] <- 0.5^(abs(i-j))
x <- mvrnorm(n, rep(0, np), Sigma)
y <- rpois(n,exp(x\%*\%beta))
glmnet.object <- glmnet(x,y,family="poisson",alpha=1)
coef.glmnet <- coef(glmnet.object)
### coefficients
coef.glmnet[ ,20]
### AIC
sAIC(x=x, y=y, beta=coef.glmnet[ ,20], family="poisson")

### Gaussian graphical model
set.seed(1)
n <- 100; np <- 10; lambda_list <- 1:100/50
invSigma <- diag( rep(0,np) )
for(i in 1:np) 
{
	for(j in 1:np)
	{
		if( i == j ) invSigma[i ,j] <- 1
		if( i == (j-1) || (i-1) == j ) invSigma[i ,j] <- 0.5
	}
}
Sigma <- solve(invSigma)
x <- scale(mvrnorm(n, rep(0, np), Sigma))
glasso.object <- glassopath(var(x), rholist=lambda_list, trace=0)
### AIC
sAIC(x=x, beta=glasso.object$wi[,,10], family="ggm")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
