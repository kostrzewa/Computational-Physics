ptm <- proc.time()

#Function: H_get. Calculate the energy of one spin. Here we applied the boundary condirion for the first and last collumn and the first and last row, otherwise we compare it to their neighbor. 
H_get<-function(y,x,s)
  {
  h=0
  if (x==1)
  {
    h=h+s[y,1]*s[y,ncol(s)]+s[y,1]*s[y,2]
  }
  else if(x==ncol(s))
  {
    h=h+s[y,ncol(s)]*s[y,ncol(s)-1]+s[y,ncol(s)]*s[y,1]
  }
  else
  {
    h=h+s[y,x]*s[y,x-1]+s[y,x]*s[y,x+1]
  }
  if(y==1)
  {
    h=h+s[1,x]*s[nrow(s),x]+s[1,x]*s[2,x]
  }
  else if(y==ncol(s))
  {
    h=h+s[ncol(s),x]*s[ncol(s)-1,x]+s[y,x]*s[1,x]
  }
  else
  {
    h=h+s[y,x]*s[y-1,x]+s[y,x]*s[y+1,x]
  }
  return(-h)
}

  
# Metropolis algorithm. We apply two conditions to flip the spins.
MA<-function(s,temp)
{
  y<-sample(1:nrow(s),1)
  x<-sample(1:ncol(s),1)
  del_H=-2*H_get(y,x,s)
  if(exp(-del_H)<1)
  {
    p_accept=exp(-del_H/temp)
  }
  else
  {
    p_accept=1  
  }
  r<-runif(1,0,1)
  if(r<p_accept)
  {
    s[y,x]=-s[y,x]
  }
  return (s)
}


#Generates 1000000 random 35*35 configruations and get image of every matrix. Without generating the images, it will spend about 17.8s for 1*10^6.  
e<-sample(c(-1,1),4*4,replace=T)
s<-matrix(e,4,4)
nc<-1000000*4*4
for (n in 1 : nc)
  {
    s=MA(s,0.0001)
    #image(s)
  }

proc.time() - ptm  
  
#Total H of systems.
get_TH<-function (S)
{
 H<-0 
 for(i in 1:nrow(S))
 {
  for (j in 1:ncol(S))
  {
   for(x in 1:ncol(S))
   {
     if (j-x==1)
     {
       H=H+S[i,j]*S[i,x]
     }
   }
 }
} 
for(j in 1:ncol(S))
{
 for (i in 1:nrow(S))
 {
   for(x in 1:nrow(S))
   {
     if (i-x==1)
     {
       H=H+S[i,j]*S[x,j]
     }
   }
 }
} 
for(n in 1:nrow(S))
{
 H=H+S[n,1]*S[n,ncol(S)]
}
for(m in 1 : ncol(S))
{
 H=H+S[1,m]*S[nrow(S),m]
}

return(-H)
}

  
#calculation of |m|
m_ab<-function(s)
{
  H<-0
  for(i in 1:nrow(s))
  {
    for(j in 1:ncol(s))
    {
      H=H+s[i,j]
    }
  }
  F=abs(H)/nrow(s)/ncol(s)
  return(F)
}

  
#calculation of m_bar, with R=J/T
m_bar<-function(N,T,a,b,R)
{
e<-sample(c(-1,1),a*b,replace=T)
s<-matrix(e,a,b)
M<-0
Mu<-0
Mb<-0
for (n in 1 : N)
{
  for(i in 1:a)
    {
    for(j in 1:b)
      {
      s=MA(s,T)
      }
    }
  Mu=Mu+exp(-R*get_TH(s))*m_ab(s)
  Mb=Mb+exp(-R*get_TH(s))
}
M=Mu/Mb
return(M)
}



#average energy per spin
E_bar_spin<-function(s)
{
  E=get_TH(s)/nrow(s)/ncol(s)
  return(E)
}
