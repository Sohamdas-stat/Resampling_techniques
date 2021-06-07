					##	Resampling Techniques 	##	Assignment
							##	Soham Das (MB2007)
##	Qn. 5.

library(tictoc)
library(MASS)
library(Rfast)

m = 6
q = 4
eta = 5
n = 100
Z = NULL
for(t in m:1){
	Z = rbind(Z,c(1, t, round(t^1.2,0), round(t^1.5)))
}
Z_inv = solve(t(Z)%*%Z)
theta = c(0,0,1,2)
sd = c(1,1.2,1.2,2)
Sigma = diag(sd) %*% matrix(c(1,0,0.25,-0.25,0,1,0,0,0.25,0,1,0.25,-0.25,0,0.25,1),q,q,byrow=T) %*% diag(sd)
sigma = 2
r_0 = NULL
for(t in 1:m){
	z = Z[t,]
	s = sqrt(as.numeric(t(z) %*% Sigma %*% z))
		x = (t(z) %*% theta - eta) / s
	r_0 = c(r_0, pnorm(x,0,1,lower.tail=T))
}

estimate <- function(n,Y){
	ybar = rowMeans(Y)
	theta_hat = Z_inv %*% t(Z) %*% ybar
	
	temp = 0
	for(i in 1:n){
		y = Y[,i]
		temp = temp + t(y) %*% y - t(y) %*% Z %*% Z_inv %*% t(Z) %*% y
	}
	sigma_hat = sqrt(temp/(n*(m-q)))
	
	temp = 0
	for(i in 1:n){
		y = Y[,i]
		vec = Z_inv %*% t(Z) %*% (y-ybar)
		temp = temp + vec %*% t(vec) - as.numeric(sigma_hat)^2 * Z_inv
	}
	Sigma_hat = temp/n
	
	r = NULL
	for(t in 1:m){
		z = Z[t,]
		s = sqrt(as.numeric(t(z) %*% Sigma_hat %*% z))
		x = (t(z) %*% theta_hat - eta) / s
		r = c(r, pnorm(x,0,1,lower.tail=T))
	}
	return(r)
}

simulate <- function(n){
	Theta = t(mvrnorm(n,theta,Sigma))
	epsilon = matrix(rnorm(n*m,0,sigma),m,n)
	Y = Z %*% Theta + epsilon
	
	## Jackknife
	
	R = NULL
	for(i in 1:n){
		jack = c(1:n)[-i]
		Y_jack = Y[,jack]
		R = rbind(R,estimate(n-1,Y_jack))
	}
	R_mean1 = colMeans(R)
	R_var1 = colVars(R)*(n-1)

	## Bootstrap
	
	R = NULL
	B = 1000
	for(b in 1:B){
		boot = sample(1:n,n,replace=T)
		Y_boot = Y[,boot]
		R = rbind(R,estimate(n,Y_boot))
	}
	R_mean2 = colMeans(R)
	R_var2 = colVars(R)
		
	return(cbind(R_mean1,R_var1,R_mean2,R_var2))
}

g <- function(X,t){
	temp = (X[1] - eta)/(sqrt(X[2]-X[1]^2-X[3]*(t(Z[t,])%*%Z_inv%*%Z[t,])))
	return(pnorm(temp,0,1,lower.tail=T))
}

grad_g <- function(X,t){
	w = (t(Z[t,])%*%Z_inv%*%Z[t,])
	s = sqrt(X[2] - X[1]^2 - X[3]*w)
	phi = dnorm((X[1]-eta)/s,0,1)
	g1 = (X[2] - X[1]*eta - X[3]*w)/s^3
	g2 = (eta - X[1])/(2*s^3)
	g3 = w*(X[1] - eta)/(2*s^3)
	g = c(g1,g2,g3)
	return(g%*%phi)
}

delta <- function(n){
	Theta = t(mvrnorm(n,theta,Sigma))
	epsilon = matrix(rnorm(n*m,0,sigma),m,n)
	Y = Z %*% Theta + epsilon
	R_hat = NULL
	var_R_hat = NULL
	
	for(t in 1:m){
		x1 = t(Z[t,]) %*% Z_inv %*% t(Z) %*% Y
		x2 = x1^2
		x3 = NULL
		for(j in 1:n){
			y = Y[,j]
			temp = t(y)%*%y - t(y) %*% Z %*% Z_inv %*% t(Z) %*% y
			x3 = cbind(x3,c(temp/(m-q)))
		}
		X = rbind(x1,x2,x3)
		X_bar = rowMeans(X)
		r_hat = g(X_bar,t)
		R_hat = c(R_hat,r_hat)
		
		temp = matrix(0,3,3)
		for(j in 1:n){
			x = X[,j]
			temp = temp + (x - X_bar) %*% t(x - X_bar)
		}
		Sigma_hat = temp/n
		
		del.g = grad_g(X_bar,t)
		var_hat = t(del.g) %*% Sigma_hat %*% del.g
		var_R_hat = c(var_R_hat,var_hat/n)
	}
	return(cbind(R_hat,var_R_hat))
}

main <- function(){
	v_boot = NULL
	v_jack = NULL
	S = 100
	
	se_L = NULL
	se_boot = NULL
	se_jack = NULL
	cat("Degradation models. Simulation output.\nPlease wait....\n0%")
	for(c in 1:S){
		cat("-")
	}
	cat("100%\n0%")
	for(c in 1:S){
		cat("=")
		se_L = cbind(se_L,sqrt(delta(n)[,2]))
		output = simulate(n)
		se_jack = cbind(se_jack,sqrt(output[,2]))
		se_boot = cbind(se_boot,sqrt(output[,4]))
	}
	cat("100%\n")
	f1 = rowMeans(se_L)
	f2 = sqrt(rowVars(se_L))
	f3 = rowMeans(se_jack)
	f4 = sqrt(rowVars(se_jack))
	f5 = rowMeans(se_boot)
	f6 = sqrt(rowVars(se_boot))
	cat("True value of R(t)\n")
	cat(r_0)
	cat("\n")
	R1 = delta(n)[,1]
	V1 = delta(n)[,2]
	output = simulate(n)
	R2 = output[,1]
	V2 = output[,2]
	R3 = output[,3]
	V3 = output[,4]
	
	tab1 = data.frame(R1,V1,R2,V2,R3,V3)
	tab2 = data.frame(f1,f2,f3,f4,f5,f6)
	names(tab1) = c("mean(R)_L","var(R)_L","mean(R)_jack","var(R)_jack","mean(R)_boot","var(R)_boot")
	names(tab2) = c("mean(se(R))_L","se(se(R))_L", "mean(se(R))_jack","se(se(R))_jack", "mean(se(R))_boot","se(se(R))_boot")
	print(tab1)
	print(tab2)
}

tic()
main()
toc()