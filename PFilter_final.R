# Bootstrap particle filter to estimate the trajectory of an object and provide
# exact approximation for the PMMH algorithm

library('abind')
library('MASS')
library('tictoc')
library('grid')
set.seed(200)
# Subfunctions to calculate distance ####
h <- function(s,x){ # Find the euclidean distance between sensor and object (for vectors)
    sqrt((x[1,]-s[1])^2+(x[2,]-s[2])^2)
}
h1 <- function(s,x){ # Find the euclidean distance between sensor and object (for scalars)
    sqrt((x[1]-s[1])^2+(x[2]-s[2])^2)
}

# Particle Filter ####
pf <- function(num_particles, Y, s1, s2, s3, f, X, Sigma,h, plot1 = F, plot2 = F, plot3 = F){
    if(plot2 == T){
        par(mfrow = c(3,1),mai = c(0.6, 0.6, 0.1, 0.5),cex.axis=1)
    }#
    # Create Stores ####
    a<-.1           # variance in dnorm (for the weighting of particles)
    n <- dim(Y)[2]  # Dimension of store based on number of timesteps
    x <- t(mvrnorm(num_particles,c(0.5,0.5,0,0),Sigma))
    Genealogies <-list()
    Genealogies[[1]] <- x
    weights <- list()
    #####
    
    for(i in 1:n){
        # Weights and normalise ####
        w1 <- dnorm(h(s1,x[1:2,]),array(Y[1,i],num_particles),a,log = F)
        w2 <- dnorm(h(s2,x[1:2,]),array(Y[2,i],num_particles),a,log = F)
        w3 <- dnorm(h(s3,x[1:2,]),array(Y[3,i],num_particles),a,log = F)
        w  <- w1 * w2 * w3
        
        weights[[i]] <- w 
        w <-  w/(sum(w)) #normalise weights
        
        # Resample ####
        I <- sample(x = num_particles,
                    size = num_particles,
                    replace = T,
                    prob = w)
        x.res <- x[,I]
        
        # Mutate ####
        x <- f%*%x.res + t(mvrnorm(num_particles,c(0,0,0,0),Sigma))
        
        # Store Genealogy ####
        if(i == 1){
            Genealogies[[i+1]] <- abind(Genealogies[[i]][,I],x,along=3 )
        }else{
            Genealogies[[i+1]] <- abind(Genealogies[[i]][,I,],x, along=3)
        }
        
        # Plot ####
        # Plot the path of the object and the particle filter estimation
        if(plot1 == T){
            par(mfrow = c(1,1),cex.axis=1)
            plot(X[1,1:i],X[2,1:i], col = 'black', xlim=c(min(X[1,])-.5,max(X[1,])+.5),
                 ylim=c(min(X[2,])-.5,max(X[2,])+.5), main = "True Object Path against Particle Filter",
                 xlab="",ylab="")
            
            points(s1[1],s1[2], pch=4)
            points(s2[1],s2[2], pch=4)
            points(s3[1],s3[2], pch=4)
            if(i==1){
                points(Genealogies[[i]][1,],Genealogies[[i]][2,],pch=20)
                points(mean(Genealogies[[i]][1,]),mean(Genealogies[[i]][2,]), col="red",pch=20)
                }
            else if(i==n){
                for(j in 2:i){
                    points(mean(Genealogies[[i]][1,,j]),mean(Genealogies[[i]][2,,j]), col="red",pch=20)
                }
            }
            else{
                points(Genealogies[[i]][1,,i],Genealogies[[i]][2,,i],pch=20)
                for(j in 2:i){
                    points(mean(Genealogies[[i]][1,,j]),mean(Genealogies[[i]][2,,j]), col="red",pch=20)
                }
            }
            layout(mat=1)
            Sys.sleep(0.2)
        }
        #plot the genealogies of the x and y coordinates
        if(plot2 == T){
            if(i == 5 ){
                # Plot the X coordinates
                plot(X[1,], ylim=c(min(X[1,])-.5,max(X[1,])+.5),
                     ylab = "X position", xlab = "",col="black",xaxt="n",
                     main = "")#
                for(j in 1:(num_particles/10)){
                    lines(Genealogies[[i]][1,j,])
                    points(array(i,length(Genealogies[[i]][1,j,i])),Genealogies[[i]][1,j,i],
                           cex = w[j]*100,pch=19,col='red')
                }
                
            }
            if(i == 15){
                # Plot the X coordinates
                plot(X[1,], ylim=c(min(X[1,])-.5,max(X[1,])+.5),xaxt="n",
                     ylab = "X position", xlab = "",col="black")
                for(j in 1:(num_particles/10)){
                    lines(Genealogies[[i]][1,j,])
                    points(array(i,length(Genealogies[[i]][1,j,i])),Genealogies[[i]][1,j,i],
                          cex = w[j]*100,pch=19,col='red')
                }
            }
            if(i == 30){
                # Plot the X coordinates
                plot(X[1,], ylim=c(min(X[1,])-.5,max(X[1,])+.5),
                     ylab = "X position", xlab = "Timestep",col="black")
                for(j in 1:(num_particles/10)){
                    lines(Genealogies[[i]][1,j,])
                    points(array(i,length(Genealogies[[i]][1,j,i])),Genealogies[[i]][1,j,i],
                           cex = w[j]*100,pch=19,col='red')
                }
            }
            Sys.sleep(.1)
        }
    }
    par(mfrow = c(1,1))

    I.path <- sample(num_particles, size = 1, replace = T, prob = weights[[n]])
    path <- Genealogies[[n]][,I.path,]
    
    list(path=path,weights=weights) 
}

# Signal and Data ####

dt <- 0.05   
n  <- 30    # number of time steps
d  <- .01   # variance in velocity

# How the object moves between timesteps
f  <- matrix(c(1,0,0,0,0,1,0,0,dt,0,1,0,0,dt,0,1),
             nrow = 4,
             ncol = 4)
# Variance and mean of signal
Sigma <- matrix(c(0,0,0,0,0,0,0,0,0,0,d,0,0,0,0,d),4,4)
mu <- c(.5,.5,0,0)

# Sensor Locations
s1 <- c(0,0)
s2 <- c(0,1)
s3 <- c(1,0)

# Generate the coordiantes and velocity of the object (Hidden State)
X <- matrix(0,4,n)
X[,1] <- mvrnorm(1,mu,Sigma)
for(i in 2:n){
    X[,i] <- f%*%X[,i-1] + mvrnorm(1,c(0,0,0,0),Sigma)
}

# Generate observed data (Distances from the sensor)
Y <- matrix(0,3,n)
for(i in 1:n){
    Y[1,i] <- h1(s1,X[,i])
    Y[2,i] <- h1(s2,X[,i])
    Y[3,i] <- h1(s3,X[,i])
}
# Add noise to observed data
for(i in 1:n){
    Y[,i] <- Y[,i] + mvrnorm(1,c(0,0,0),diag(c(d,d,d)))
}

# Run Particle Filter ####
num_particles = 200
G <- pf(num_particles, Y, s1, s2, s3, f, X, Sigma,h, plot1 = F, plot2 = F, plot3 = F)

# Plot of object against BPF path
par(pty="s")
plot(X[1,],X[2,], col = 'black', xlim=c(0,1.1),
     ylim=c(0,1.1), main = "",
     xlab="x-coordinate",ylab="y-coordinate")
points(s1[1],s1[2], pch=4)
points(s2[1],s2[2], pch=4)
points(s3[1],s3[2], pch=4)
points(G$path[1,],G$path[2,], col="red",pch=20)
legend(max(X[1,]),max(X[2,])+.6,
       legend = c("Object","BPF Estimate   ","Sensor"),
       col=c("black","red","black"),pch=c(1,20,4))

mse <- matrix(0,nrow=6,ncol=100)
bias <- matrix(0,nrow=6,ncol=100)

for(j in 1:100){
    paths <- list()
    i=0
    for(num_particles in c(100,200,400,800,1600,2000)){
        i=i+1
        paths[[i]] <- pf(num_particles, Y, s1, s2, s3, f, X, Sigma,h, plot1 = F, plot2 = F, plot3 = F)$path
        mse[i,j]<-mean((paths[[i]][1,]-X[1,])^2+(paths[[i]][2,]-X[2,])^2)
        bias[i,j]<-mean(sqrt((paths[[i]][1,]-X[1,])^2+(paths[[i]][2,]-X[2,])^2))
    }
}
mse <- apply(mse,MARGIN = 1,mean)
bias <- apply(bias,MARGIN = 1,mean)
mse/mse[1]
bias/bias[1]

G <- pf(num_particles, Y, s1, s2, s3, f, X, 5*Sigma,h, plot1 = F, plot2 = T, plot3 = F)

