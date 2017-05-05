Roots.of.unity <- function(nth.root = N, radius=1, simultaneous.reals=TRUE,plot.roots=TRUE){
  
  Simultaneous.reals = data.frame()
  
  Complex.points = numeric()
  Reals = numeric()
  Imaginaries = numeric()
  
  r = radius^(1/nth.root)
  
  
  for(j in 0:(nth.root-1)){
    theta= ((360*j)/nth.root)*pi/180

    Simultaneous.reals[(j+1),1] = r*cos(theta) - r*sin(theta)
    Simultaneous.reals[(j+1),2] = r*cos(theta) + r*sin(theta)
    
    Complex.points[j+1] = complex(real = r*cos(theta), imaginary = r*sin(theta))
    Complex.points = c(Complex.points,Complex.points[j+1])

    Reals[j+1] = r*cos(theta)
    Imaginaries[j+1] = r*sin(theta)
      
  }
  
  if(plot.roots==TRUE){
  circleFun <- function(center = c(0,0),diameter = 2*radius, npoints = 100){
    r = diameter / 2
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
  }
  dat <- circleFun(c(0,0),2*radius^(1/nth.root),npoints = 100)

  plot(dat,type = 'l',main = substitute(paste("roots of " , "Z"^nth.root ,"= ", Z),list(nth.root=nth.root,Z=radius)),
       xlab = "Real", ylab = "Imaginary",
       xlim = c(min(Simultaneous.reals[,1]),max(Simultaneous.reals[,1])),
       ylim = c(min(Simultaneous.reals[,1]),max(Simultaneous.reals[,1]))
       )
  if(simultaneous.reals==TRUE){
      legend('topleft',c("Complex Points","Simultaneous Reals"),pch = c(19,19), col=c('red','blue'),bty='n')
  
      points(Simultaneous.reals[,1],rep(0,length(Simultaneous.reals[,1])),
         pch=19,
         col='blue')
      points(Reals,Imaginaries,col='red',pch=19)
     
       
  } else {
  
           legend('topleft',"Complex Points",pch = 19, col='red',bty='n')
           
  points(Reals,Imaginaries,col='red',pch=19)}
  }
  colnames(Simultaneous.reals) = c("Real 1","Real 2")
  
  list(c("Reals"=Reals))
  
  return(cbind(Simultaneous.reals,"Complex Points" = Complex.points[1:nth.root]))
 
  
}


#Factor check
Roots.check <- function(x,y){
  a<- Roots.of.unity(x,plot.roots = FALSE)[-1,1]
  b<- Roots.of.unity(y,plot.roots = FALSE)[-1,1]
  
plot(a,rep(-0.01,x-1),pch="|",bg=ifelse(a%in%b,'red','white'), main = "Simultaneous Reals for Roots of Unity",ylim = c(-.1,.1),cex=3,ylab='',xlab = "Real Numbers",
     yaxt='n')
points(b,rep(0.01,y-1),pch="|",col='blue',cex=3)
points(a,rep(0,length(a%in%b)),pch=19,col=ifelse(a%in%b,'red','transparent'))
legend('top',legend=list(y,x), pch = "|",cex=2,col = c('blue','black'),bty='n')
  }



DFT <- function(x){
  l = length(x)
  f.hat = numeric()
  Simultaneous.reals = data.frame()
  magnitude = numeric()  
  
  for(i in 1:(l)){
    
    f.hat[i]=sum(((Roots.of.unity(l,plot.roots = FALSE)$"Complex Points")^(i-1))*x)
    Simultaneous.reals[i,1] = Re(f.hat[i]) - Im(f.hat[i])
    Simultaneous.reals[i,2] = Re(f.hat[i]) + Im(f.hat[i])
    if(i <= l/2){
    magnitude[i] = sqrt(Re(2*f.hat[i])^2 + Im(2*f.hat[i])^2)} else {
      magnitude[i] = 0 
    }
  }
  
  Amplitude = magnitude/l
  par(mfrow=c(2,1))
  plot((1:(l/2))-1,Amplitude[1:(l/2)],col='steelblue',pch=19,xlab="Frequency",ylab = "Amplitude")
  
  plot(Simultaneous.reals[,1],rep(0,l),pch="|",cex=2,col='blue',ylab="",xlab="Simultaneous Reals")
  
  colnames(Simultaneous.reals) = c("Real 1","Real 2")
  return(cbind(Simultaneous.reals,"f hat"=f.hat, "Magnitude"=format(magnitude,digits = 4),"Amplitude"=magnitude/l))
  
}  
