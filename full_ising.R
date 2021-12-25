#Metropolis-Hastings
J=1
l=10
iters=10000
temp=1
calc_Energy<-function(config){
  sum1=0
  
  for(i in 1:l){
    for(j in 1:l){
      if(i==1){
        if(j==1){
          sum1=config[i,j+1]*config[i,j]
          sum1=config[i+1,j]*config[i,j]
        }
        else if(j==l){
          sum1=config[i,j-1]*config[i,j]
          sum1=config[i+1,j]*config[i,j]
        }
        else{
          sum1=config[i,j-1]*config[i,j]
          sum1=config[i,j+1]*config[i,j]
          sum1=config[i+1,j]*config[i,j]
        }
        
        
      }
      else if(i==l){
        if(j==1){
          sum1=config[i,j+1]*config[i,j]
          sum1=config[i-1,j]*config[i,j]
        }
        else if(j==l){
          sum1=config[i,j-1]*config[i,j]
          sum1=config[i-1,j]*config[i,j]
        }
        else{
          sum1=config[i,j-1]*config[i,j]
          sum1=config[i,j+1]*config[i,j]
          sum1=config[i-1,j]*config[i,j]
        }
        
        
      }
      else if(j==1){
        sum1=config[i,j+1]*config[i,j]
        sum1=config[i-1,j]*config[i,j]
        sum1=config[i+1,j]*config[i,j]
      }
      else if(j==l){
        sum1=config[i,j-1]*config[i,j]
        sum1=config[i-1,j]*config[i,j]
        sum1=config[i+1,j]*config[i,j]
      }
      else{
        sum1=config[i,j-1]*config[i,j]
        sum1=config[i,j+1]*config[i,j]
      }
    }}
  -J*(sum1)
}

lattice<-function(l){
  MArray = array(c(0), dim=c(l,l))
  
  for(x in 1:10){
    
    MArray[x,1:l]=runif(10)
    MArray[x,1:l]=(MArray[x,1:l]<0.5)*2-1
    
    
    
    
  }
  
  
  MArray
}
alg<-function(lattice,iters,temp){
  e_array=c(calc_Energy(lattice))
  mag_array=c(sum(lattice)/(l**2))
  for(n in 2:iters){
    result<-mh(e_array[n-1],lattice,temp)
    e_array[n]=result$energy
    lattice=result$config
    print("result")
    print(result$config)
    
    mag_array[n]=sum(lattice)/(l**2)
  }
  mag_array
}
mh<-function(init_e,state,temp){                 
  x<-round(runif(1, 1, l))
  y<-round(runif(1,1,l))
  print("x")
  print(x)
  print("y")
  print(y)
  print(state)
  print("x,y")
  print(state[x,y])
  state[x,y]=-1*state[x,y]
  print("x,y")
  print(state[x,y])
  
  energy<-calc_Energy(state)-init_e
  if (energy<0){
    return(list(config=state,energy=calc_Energy(state)))
  }
  if (runif(1,0,1)<exp(-1*energy/temp)){
    return(list(config=state,energy=calc_Energy(state)))
  }
  state[x,y]=-1*state[x,y]
  return(list(config=state,energy=init_e))
  
}
state<-lattice(l) #for 1st chain
state2<-lattice(l)#for second chain
energy<-calc_Energy(state)#for first chain
energy2<-calc_Energy(state2)#for second chain
results<-alg(state,iters,temp)#for first chain
results2<-alg(state2,iters,temp)#for second chain
plot(results,main='Metropolis Trace Plot chain 1',xlab='Iteration #',ylab='Magnetization')
plot(results2,main='Metropolis Trace Plot chain 2',xlab='Iteration #',ylab='Magnetization')


acf(results, main = 'Metropolis, chain 1')
acf(results2, main = 'Metropolis, chain 2')

d<-500
l<-10000-500
mean<-mean(results[500:length(results)])#for first chain
mean2<-mean(results2[500:length(results2)])#for second chain
hist(results[501:length(results)],breaks=100,main='Histogram (Metropolis),Chain 1',xlab='Magnetization Estimates')#for first chain
hist(results2[501:length(results2)],breaks=100,main='Histogram (Metropolis), Chain 2',xlab='Magnetization Estimates')#for second chain


chainmean=(mean+mean2)/2
btw_chain_variance=((length(results)-500))*((mean-chainmean)**2+(mean2-chainmean)**2)
var1=var(results[500:length(results)])
var2=var(results[500:length(results2)])
w=(1/2)*(var1+var2)
r=(((l-1)/l)*w+(1/l)*btw_chain_variance)/(w)

#Gibbs Sampling

J=1
l=10
iters=10000
temp=1
calc_Energy<-function(config){
  sum1=0
  
  for(i in 1:l){
    for(j in 1:l){
      if(i==1){
        if(j==1){
          sum1=config[i,j+1]*config[i,j]
          sum1=config[i+1,j]*config[i,j]
        }
        else if(j==l){
          sum1=config[i,j-1]*config[i,j]
          sum1=config[i+1,j]*config[i,j]
        }
        else{
          sum1=config[i,j-1]*config[i,j]
          sum1=config[i,j+1]*config[i,j]
          sum1=config[i+1,j]*config[i,j]
        }
        
        
      }
      else if(i==l){
        if(j==1){
          sum1=config[i,j+1]*config[i,j]
          sum1=config[i-1,j]*config[i,j]
        }
        else if(j==l){
          sum1=config[i,j-1]*config[i,j]
          sum1=config[i-1,j]*config[i,j]
        }
        else{
          sum1=config[i,j-1]*config[i,j]
          sum1=config[i,j+1]*config[i,j]
          sum1=config[i-1,j]*config[i,j]
        }
        
        
      }
      else if(j==1){
        sum1=config[i,j+1]*config[i,j]
        sum1=config[i-1,j]*config[i,j]
        sum1=config[i+1,j]*config[i,j]
      }
      else if(j==l){
        sum1=config[i,j-1]*config[i,j]
        sum1=config[i-1,j]*config[i,j]
        sum1=config[i+1,j]*config[i,j]
      }
      else{
        sum1=config[i,j-1]*config[i,j]
        sum1=config[i,j+1]*config[i,j]
      }
    }}
  -J*(sum1)
}

lattice<-function(l){
  MArray = array(c(0), dim=c(l,l))
  
  for(x in 1:10){
    
    MArray[x,1:l]=runif(10)
    MArray[x,1:l]=(MArray[x,1:l]<0.5)*2-1
    
    
    
    
  }
  
  
  MArray
}
alg<-function(lattice,iters,temp){
  e_array=c(calc_Energy(lattice))
  mag_array=c(sum(lattice)/(l**2))
  for(n in 2:iters){
    result<-mh(e_array[n-1],lattice,temp)
    e_array[n]=result$energy
    lattice=result$config
    print("result")
    print(result$config)
    
    mag_array[n]=sum(lattice)/(l**2)
  }
  mag_array
}
mh<-function(init_e,state,temp){                 
  x<-round(runif(1, 1, l))
  y<-round(runif(1,1,l))
  print("x")
  print(x)
  print("y")
  print(y)
  print(state)
  print("x,y")
  print(state[x,y])
  state[x,y]=-1*state[x,y]
  print("x,y")
  print(state[x,y])
  
  energy<-calc_Energy(state)-init_e
  if (energy<0){
    return(list(config=state,energy=calc_Energy(state)))
  }
  if (runif(1,0,1)<(exp(-1*energy/temp)/(exp(-1*energy/temp)+exp(1*energy/temp)))){
    return(list(config=state,energy=calc_Energy(state)))
  }
  state[x,y]=-1*state[x,y]
  return(list(config=state,energy=init_e))
  
}
state<-lattice(l)#for first chain
state2<lattice(l)#for second chain
energy<-calc_Energy(state)#for first chain
energy2<-calc_Energy(state)#for second chain

res<-alg(state,iters,temp)#for first chain
res2<-alg(state,iters,temp)#for second chain
print(res)
plot(res,main='Gibbs Trace Plot chain 1',xlab='Iteration #',ylab='Magnetization')
plot(res2,main='Gibbs Trace Plot chain 2',xlab='Iteration #',ylab='Magnetization')

res
acf(res,main='Gibbs Sampling, chain 1')
acf(res2,main='Gibbs Sampling, chain 2')
mean<-mean(res[500:length(results)])#for first chain
mean2<-mean(res[500:length(results)])#for second chain


hist(res[501:(length(results))],breaks=100,main='Histogram (Gibbs Sampling), Chain 1',xlab='Magnetization Estimates')
hist(res2[501:(length(results))],breaks=100,main='Histogram (Gibbs Sampling), Chain 2',xlab='Magnetization Estimates')
l<-length(res)-500

chainmean=(mean+mean2)/2
btw_chain_variance=((length(res)-500))*((mean-chainmean)**2+(mean2-chainmean)**2)
var1=var(results[500:length(res)])
var2=var(results[500:length(res2)])
w=(1/2)*(var1+var2)
r=(((l-1)/l)*w+(1/l)*btw_chain_variance)/(w)



