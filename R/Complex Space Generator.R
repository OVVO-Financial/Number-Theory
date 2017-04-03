Complex.Space.Generator <- function(N){
  
  top = ceiling(sqrt(N))
  bottom = floor(sqrt(N))
  middle = mean(c(top,bottom))
  coordinates = data.frame(ncol=2)
  
  segtop_x = ceiling(N/6)+1
  segtop_y = floor(N/6)-1

  par(mfrow=c(1,1))
  
  
  #These are the complex coordinates, for example (48 + 15i) for N=1219
  r=3:as.integer(N/3)
  im=0:as.integer(N/3)
  plane=expand.grid(r,im)
  
  ### Reduce plane
  plane=plane[(plane[,1]-plane[,2])>=3,]
  plane=plane[(plane[,1]+plane[,2])<=(N/3),]
  plane=plane[(plane[,1]-plane[,2])<=top,]
  
  complex.plane=complex(real=plane[,1],imaginary = plane[,2])
  

  plot(plane,pch = 17,
          col = ifelse((plane[,1]-plane[,2])*(plane[,1]+plane[,2]) < N,'blue',
          ifelse((plane[,1]-plane[,2])*(plane[,1]+plane[,2]) == N,'green', 
          ifelse(plane[,1]==top && plane[,2]==(top-3),'orange' , 'red') )),
          xlab = "Real",ylab="Imaginary",main=paste("N= ",N,sep = "")
          )
  
  

  
}  
