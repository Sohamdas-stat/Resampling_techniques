					##	Resampling Techniques 	##	Assignment
							##	Soham Das (MB2007)
##	Qn. 4.

library(tictoc)
library(LaplacesDemon)
library(rootSolve)

low.r <- function(x) return(max(2,x))
low.s <- function(x) return(max(3,x))
low.u <- function(x) return(max(4,x))

income.gen <- function(urb.cat){
	inc.fac = 40
	sigma = c(1,0.7,0.5)
	x = NULL
	for(i in urb.cat){
		y = rlnormp(1, 0.6, sigma[i]) * inc.fac / 2
		if(i==1) y = sapply(y,low.r)
		if(i==2) y = sapply(y,low.s)
		if(i==3) y = sapply(y,low.u)
		r = (2 - pnorm(10*y/inc.fac,5,2,lower.tail=T)) / 2
		z = y * r
		x = cbind(x,c(y,z))
	}
	return(x)
}

family.gen <- function(urb.cat){
	lambda = c(5,4.5,4)
	x = NULL
	for(i in urb.cat){
		y = rpois(1,lambda[i])
		while(y == 0 || y > 8){
			y = rpois(1,lambda[i])
		}
		x = c(x,y)
	}
	return(x)
}

estimate <- function(Data,fam,urb){
	LM = lm(log(expend) ~ 1 + log(income) + family.size + as.factor(urb.cat), data=Data)
	gamma = mean(Data$expend)/mean(Data$income)
	gamma = c(gamma,LM$coef[1:3])
	beta = c(0,LM$coef[4:5])
	fun <- function (theta) (log(gamma[1]*theta + 0.2) - (gamma[2] + gamma[3] * log(theta) + gamma[4]*fam + beta[urb]))
	x = seq(1,100,by=10)
	uni <- as.vector(uniroot.all(fun, c(1,100)))	####### attention
	if(length(uni)>1) print(lengh(uni))
	return(uni)
}

simulation <- function(n,req.fam.size,req.urb.cat,B){
	urb.cat = sample(c(1:3),n,replace=T)
	income = income.gen(urb.cat)
	expend = income[2,]
	income = income[1,]
	family.size = family.gen(urb.cat)
	Data = data.frame(urb.cat, income, expend, family.size)
	theta_hat = estimate(Data,req.fam.size,req.urb.cat)
	
	##	Bootstrap
	
	theta_boot = NULL
	for(b in 1:B){
		boot = sample(1:n,n,replace=T)
		Dat = Data[boot,]
		theta_boot = c(theta_boot,estimate(Dat,req.fam.size,req.urb.cat))
	}
	m_boot = mean(theta_boot)
	v_boot = var(theta_boot)
	
	##	Jackknife
	
	theta_jack = NULL
	for(i in 1:n){
		jack = c(1:n)[-i]
		Dat = Data[jack,]
		theta_jack = c(theta_jack,estimate(Dat,req.fam.size,req.urb.cat))
	}
	m_jack = mean(theta_jack)
	v_jack = var(theta_jack)*(n-1)
	return(c(sqrt(v_boot),sqrt(v_jack)))
}

main <- function(){
	n = 500
	B = 1000
	se_boot = NULL
	se_jack = NULL
	S = 100
	req.fam.size = 6
	req.urb.cat = 1
	cat("Estimation of poverty line. Simulation output.\nPlease wait....\n")
	for(i in 1:S){
		cat("=")
		output = simulation(n,req.fam.size,req.urb.cat,B)
		se_boot = c(se_boot,output[1])
		se_jack = c(se_jack,output[2])
	}
	cat("\nmean(se(theta))_boot = ",mean(se_boot),"\n")
	cat("se(se(theta))_boot = ",sqrt(var(se_boot)),"\n")
	cat("mean(se(theta))_jack = ",mean(se_jack),"\n")
	cat("se(se(theta))_jack = ",sqrt(var(se_jack)),"\n")
	
	return(data.frame(se_boot,se_jack))
}

tic()
output = main()
toc()