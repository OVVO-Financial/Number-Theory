Newton_attractor<- function(constant, tol=1E-12,x0=0,N=500) {

### X^2 + C = 0
f<- function(x) {x^2 + constant}

### Coefficients of (i) for plot reference only!
solution <- sqrt(abs(constant))

    h <- .01
    i <- 1; x1 <- x0
    p <- numeric(N)
    while (i<=N) {
      df.dx <- (f(x0+h)-f(x0))/h
      x1 <- (x0 - (f(x0)/df.dx))
      p[i] <- x1
      i <- i + 1
      if (abs(x1-x0) < tol) break
      x0 <- x1
    }

### Newton-Raphson estimates 
    Est=p[1:(i-1)]
    print(Est)
    
    L<- length(Est)-1

### Plot estimates
    for(i in 1:L){
      par(mfrow=c(1,2))
      
### Plot total estimates
      plot(Est[1:i],Est[2:(i+1)],type = "l",col="blue", main = "Total",xlim=c(min(Est),max(Est)),ylim=c(min(Est),max(Est)),
           xlab="X1",ylab = "X0")
      abline(v = (solution), col="red")
      abline(v = -(solution), col="red")
      
      
### Zoom in on solution area
      plot(Est[1:i],Est[2:(i+1)],type = "l",col="blue", main = "Zoom", 
           xlab="X1",ylab = "X0",
           xlim=c(-(2*solution), 2*solution)
           ,ylim=c(-(2*solution), 2*solution))
      abline(v = -(solution), col="red")
      abline(v = (solution), col="red")

### Watch iterations      
      Sys.sleep(.05)
      
      
    }
    
}