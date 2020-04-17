# Particle marginal Metropolis-Hastings algorithm to locate the x-coordinate 
# of the "third" sensor in an example tracking problem
# Note: Code requires PFilter.R to have been run first.

set.seed(200)
# Prior Normal distribution
log.prior<-function(theta,prior.mean,prior.sd){
  l.prior<-dnorm(theta,prior.mean,prior.sd,log = TRUE)
  l.prior
}

prior.theta <- function(theta) {
  return(dnorm(theta,mean=1,sd=0.1))  
}

my.q <- function(theta.old,theta.new) {
  q.sd <- 0.05;
  return(dnorm(theta.new,mean=theta.old,sd=q.sd))
}

# Data
log.likelihood<-function(theta,X){ 
  l.lik<-dnorm(X,theta,log=TRUE)
  l.lik
}

log.posterior<-function(theta,X,prior.mean,prior.sd){
  l.prior<-log.prior(theta,prior.mean,prior.sd)
  l.lik<-log.likelihood(theta,X)
  l.post<-l.prior+l.lik # this is Bayes theorm in log form
  l.post
}

PMH<-function(theta0=1,sigma.prop,N,prior.mean,prior.sd,num_particles){
  start <- Sys.time ()
  # Theta <- array(0,c(N.pmcmc,1))
  # Paths <- list()
  ILikes <- array(0,c(N,1))
  ar <- array(0,N)
  theta<-array(0,N)
  phi <- array(0,c(120,N))
  
  # This is essentially the first sample for theta
  theta[1]<-theta0  
  
  # Run the particle filter
  G <- pf(num_particles, Y, s1, s2, s3=c(theta0,0), f, X, Sigma,h, plot1 = F)
  phi[,1] <- as.vector(G$path)
  ilike <- 1
  for( i in 1:n ) {
    ilike <- ilike * sum(G$weights[[i]])/num_particles
  }
  cat("Marginal likelihood: ",ilike,"\n")
  ILikes[1] <- ilike
  
  alpha<-rep(0,N) 
  alpha[1] <- 1
  for (i in 2:N){ # MH loop
    if(i%%1000 == 0 ){
        cat(i,"\n");
        }
    
    current_theta<-theta[i-1]
    
    # Sample theta
    proposed_theta<-current_theta + rnorm(1,sd=sigma.prop)   
    # Run the particle filter
    G <- pf(num_particles, Y, s1, s2, s3=c(proposed_theta,0), f, X, Sigma,h)
    
    # Here we do the marginal likelihood and sampling steps just as before
    ilike <- 1
    for( j in 1:n ) {
      ilike <- ilike * sum(G$weights[[j]])/num_particles
    }
    
    # Acceptance test
    accept.prob <- # eqn (12) 
      min(1, 
          (ilike*prior.theta(proposed_theta)) / (ILikes[i-1]*prior.theta(current_theta)) * 
            my.q(current_theta,proposed_theta) / my.q(proposed_theta,current_theta))
    if (is.na(accept.prob)) {
      print("Something is wrong")
    }
    
    alpha[i] <- accept.prob
    u <- runif(1);
    if (u < accept.prob) {
      ILikes[i] <- ilike;
      theta[i] <- proposed_theta;
      phi[,i] <- as.vector(G$path)
      ar[i] <- 1
    } else {
      ILikes[i] <- ILikes[i-1];
      theta[i] <- current_theta;
      phi[,i] <- phi[,i-1]
    }
  }
  cat("Time per loop: ",difftime(Sys.time(), start, units = "secs")[[1]]/N, "secs. \n")
  cat("Total Time: ",difftime(Sys.time(), start, units = "mins")[[1]], "Mins. \n")
  list(theta=theta,alpha=alpha,ar=ar,phi=phi)
}

theta0 = 0
N.MCMC = 1000 #number of simulations
prior.mean = 1
prior.sd = .1
sigma.prop <- 0.1
num_particles=2000
out<-PMH(theta0,sigma.prop,N.MCMC,prior.mean,prior.sd,num_particles)
mean(out$ar)

time<-rep(0, 8)
num_particles <- c(100,200,400,800,1600,3200,6400,12800)
i=0
for(d in num_particles){
    i=i+1
    start<-Sys.time()
    out<-PMH(theta0,sigma.prop,100,prior.mean,prior.sd,d)
    time[i] <- difftime(Sys.time(), start, units = "secs")[[1]]
}
plot(num_particles,time,type="l")
lines(num_particles,num_particles*(time[7]/num_particles[7]),col="red")

# Find prop.sd
accept <- matrix(0,nrow = 5,ncol = 6)
i=0

for(num_particles in c(100,200,400,800,1600,2000)){
    i=i+1
    j=0
    for(sigma.prop in c(0.01,0.05,0.10,0.15,0.20)){
        j=j+1
        print(sigma.prop)
        print(num_particles)
        accept[j,i] <- mean(PMH(theta0,sigma.prop,N.MCMC,prior.mean,prior.sd,num_particles)$ar)
    }
}
accept

par(mfrow=c(1,1),mai = c(0.8, 0.7, 0.2, 0.2), cex.axis=1)
layout(1)
# Check burnin
plot(out$theta[1:1000],type = "l",
     ylab = expression(theta),
     xlab = "Timestep")

burn = 500
theta_burn <- out$theta[burn:length(out$theta)]

#ACF and ess for theta
par(mfrow = c(1,2),mai = c(0.8, 0.7, 0.2, 0.2), mgp=c(2.4,1,0))
max_lag = max(which(acf(theta_burn,plot=F)[[1]]>2/sqrt(length(theta_burn))))
acf(theta_burn,lag.max = max_lag,
    main = "",xlab=expression(tau))
ACL=2*sum(acf(theta_burn,plot=F)[[1]])+1
(ess = length(theta_burn)/ACL)
theta_thin <- theta_burn[seq(1,length(theta_burn),ACL)]
acf(theta_thin,main="",ylab="",lag.max=max_lag,xlab=expression(tau))
ks.test(theta_thin[1:floor(length(theta_thin)/2)],
        theta_thin[(1+ceiling(length(theta_thin)/2)):length(theta_burn)])

#Histograms of samples to visually check ks test
par(mfrow=c(1,1),mai = c(0.8, 0.9, 0.5, 0.4))
hist(theta_thin[1:floor(length(theta_thin)/2)],col=rgb(1,0,0,0.5),breaks = 30,
     main = "",
     xlab=expression(theta),freq = F )
hist(theta_thin[(1+ceiling(length(theta_thin)/2)):length(theta_thin)],col=rgb(0,0,1,0.5),breaks=30,freq = F,add=T)
legend("topright",legend = c("1st Half","2nd Half"),col=c("red","blue"),pch=20)


# Output for sensor x-coordinate
par(mfrow = c(2,1),mai = c(0.7, 0.6, 0.2, 0.2), mgp=c(2,1,0))
plot(theta_thin,type = "l",
     ylab = expression(theta),
     xlab = "",
     xaxt="n")

xtick<-seq(burn, N.MCMC, by=1000)
axis(side=1, at=(xtick-burn), labels = FALSE)
text(x=(xtick-burn),  par("usr")[3], 
     labels =xtick, pos = 1, xpd = TRUE,cex=1)

plot(cumsum(theta_thin)/1:length(theta_thin),type="l",main="",xlab="Timestep",
     ylab=expression(theta),xaxt="n")
abline(h=cumsum(theta_thin)[length(theta_thin)]/length(theta_thin),col="red")
axis(side=1, at=(xtick-burn), labels = FALSE)
text(x=(xtick-burn),  par("usr")[3], 
     labels =xtick, pos = 1, xpd = TRUE,cex=1)
par(mfrow = c(1,1), mai = c(.8,1,0.3,.8), cex.axis=1)
hist(theta_thin,breaks = 20, freq = F,
     main = "",
     xlab=expression(theta))
abline(v=1,col='red')
# curve(dnorm(x, mean=mean(theta_thin), sd=sqrt(var(theta_thin))), col="darkblue",
#             lwd=2, add=TRUE, yaxt="n")
abline(v=mean(theta_thin)+2*sqrt(var(theta_thin)),lty=2)
abline(v=mean(theta_thin)-2*sqrt(var(theta_thin)),lty=2)
legend("topright", legend = c("True location", "2 s.d."), col=c("red","black"),lty=c(1,2))

# qqplot(theta_thin,rnorm(length(theta_thin),mean(theta_thin),sqrt(var(theta_thin))))
# xy=seq(0.9,1.15,0.01)
# lines(xy,xy,col="red")

layout(1)

phi_burn <- out$phi[57,burn:length(out$theta)]

#ACF and ess for object
par(mfrow = c(1,2),mai = c(0.8, 0.7, 0.2, 0.2), mgp=c(2.4,1,0),cex.axis=1)
max_lag = max(which(acf(phi_burn,plot=F)[[1]]>2/sqrt(length(phi_burn))))
acf(phi_burn,lag.max = max_lag,
    main = "",xlab=expression(tau))
ACL=2*sum(acf(phi_burn,plot=F)[[1]])+1
(ess = length(phi_burn)/ACL)
phi_thin <- phi_burn[seq(1,length(phi_burn),ACL)]
acf(phi_thin,main="",ylab="",lag.max=max_lag,xlab=expression(tau))
ks.test(phi_thin[1:floor(length(phi_thin)/2)],
        phi_thin[(1+ceiling(length(phi_thin)/2)):length(phi_burn)])

#Histograms of samples to visually check ks test
par(mfrow=c(1,1),mai = c(0.8, 0.7, 0.5, 0.2))
hist(phi_thin[1:floor(length(phi_thin)/2)],col=rgb(1,0,0,0.5),breaks = 30,
     main = "",
     xlab=expression(x[15]),freq = F )
hist(phi_thin[(1+ceiling(length(phi_thin)/2)):length(phi_thin)],col=rgb(0,0,1,0.5),breaks=30,freq = F,add=T)
legend("topright",legend = c("1st Half","2nd Half"),col=c("red","blue"),pch=20)


# Output for object x-coordinate
par(mfrow = c(2,1),mai = c(0.7, 0.6, 0.2, 0.2), mgp=c(2,1,0))
plot(phi_thin,type = "l",
     ylab = expression(x[15]),
     xlab = "",
     xaxt="n")


xtick<-seq(burn, 4500, by=1000)
axis(side=1, at=(xtick-burn), labels = FALSE)
text(x=(xtick-burn),  par("usr")[3], 
     labels =xtick, pos = 1, xpd = TRUE,cex=1)

plot(cumsum(phi_thin)/1:length(phi_thin),type="l",main="",xlab="Timestep",
     ylab=expression(x[15]),xaxt="n")
abline(h=cumsum(phi_thin)[length(phi_thin)]/length(phi_thin),col="red")
axis(side=1, at=(xtick-burn), labels = FALSE)
text(x=(xtick-burn),  par("usr")[3], 
     labels =xtick, pos = 1, xpd = TRUE,cex=1)
par(mfrow = c(1,1), mai = c(.8,1,0.3,.8), cex.axis=1)
hist(phi_thin,breaks = 20, freq = F,
     main = "",
     xlab=expression(x[15]))
# curve(dnorm(x, mean=mean(phi_thin), sd=sqrt(var(phi_thin))), col="darkblue",
#       lwd=2, add=TRUE, yaxt="n")
abline(v=X[1,15],col='red')
abline(v=mean(phi_thin)+2*sqrt(var(phi_thin)),lty=2)
abline(v=mean(phi_thin)-2*sqrt(var(phi_thin)),lty=2)
legend("topright", legend = c("True location", "2 s.d."), col=c("red","black"),lty=c(1,2))

# qqplot(phi_thin,rnorm(length(phi_thin),mean(phi_thin),sqrt(var(phi_thin))))
# xy=seq(0.55,0.7,0.01)
# lines(xy,xy,col="red")
# ,
#      yaxt="n",
#      ylab="")

#check the acf of the plot
traj = apply(out$phi[,burn:length(out$theta)],MARGIN = 1,FUN = mean)
traj=matrix(traj,nrow=4,ncol=30)

traj_sd = apply(out$phi[,burn:length(out$theta)],MARGIN = 1,FUN = sd)
traj_sd=matrix(traj_sd,nrow=4,ncol=30)

par(mfrow = c(1,1), mai = c(0.9,0.9,0.2,0.1), cex.axis=1)
plot(X[1,],X[2,], col = 'black', xlim=c(0,1.2),
     ylim=c(0,1.2), main = "",
     xlab="x-coordinate",ylab="y-coordinate")
points(s1[1],s1[2], pch=4)
points(s2[1],s2[2], pch=4)
points(s3[1],s3[2], pch=4)
points(mean(theta_thin),0,col="red",pch=4)
points(mean(theta_thin)+2*sd(theta_thin),0,col="red",pch=124)
points(mean(theta_thin)-2*sd(theta_thin),0,col="red",pch=124)
points(traj[1,],traj[2,], col="red",pch=20)
lines(traj[1,]+2*traj_sd[1,],traj[2,]+2*traj_sd[2,],col="red",lty=2)
lines(traj[1,]-2*traj_sd[1,],traj[2,]-2*traj_sd[2,],col="red",lty=2)
legend(max(X[1,]),max(X[2,])+.7,
       legend = c("Sensor","Sensor Estimate    ","Object","Object Estimate  ","2 s.d."),
       col=c("black","red","black","red","red"),pch=c(4,4,1,20,45))

plot(X[1,],ylim=c(0.45,0.9),ylab="x-coordinate",xlab="Timestep")
points(traj[1,],col="red",pch=20)
lines(traj[1,]+2*traj_sd[1,],col="red",lty=2)
lines(traj[1,]-2*traj_sd[1,],col="red",lty=2)
legend(1,0.9,
       legend = c("Object","Object Estimate","2 s.d."),
       col=c("black","red","red"),pch=c(1,20,45))

# Final Plot zoomed
plot(X[1,],X[2,], col = 'black', xlim=c(min(X[1,])-.1,max(X[1,])+.1),
     ylim=c(min(X[2,])-.2,max(X[2,])+.2), main = "",
     xlab="x-coordinate",ylab="y-coordinate")
points(traj[1,],traj[2,], col="red",pch=20)
lines(traj[1,]+2*traj_sd[1,],traj[2,]+2*traj_sd[2,],col="red",lty=2)
lines(traj[1,]-2*traj_sd[1,],traj[2,]-2*traj_sd[2,],col="red",lty=2)
legend(min(X[1,])+0.2,max(X[2,])+.2,
       legend = c("Object","Object Estimate     ","2 s.d."),
       col=c("black","red","red"),pch=c(1,20,45))

