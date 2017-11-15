Function of ising model (a3)

#get Hamiltonian coefficient  Ham= H*J
get_H<-function (S){           
H<-0 #get initial value 0
for(i in 1:nrow(S)){for (j in 1:ncol(S)){for(x in 1:ncol(S)){if (j-x==1){H=H+S[i,j]*S[i,x]}}}} #sum all the sisj in row 
for(j in 1:ncol(S)){for (i in 1:nrow(S)){for(x in 1:nrow(S)){if (i-x==1){H=H+S[i,j]*S[x,j]}}}} #sum all the sisj in col
return(-H)}# H=-sisj



                                                                                     #return H=-sisj.






#function: get_random:get a random matrix  if not same with before, add the system to the matrix. 
#in the matrix every row repreasents a system.

get_random<-function(x){
l<-0                                             # get initial l 0
r<-0                                             # get initial r 0
h<-sample(c(-1,1),ncol(x),replace=T)             # get a random system-t
t<-matrix(h,nrow=1,ncol= ncol(x))                
for(j in 1:nrow(x))                              # compared t to every existed system
{
for(k in 1:ncol(x)){if (x[j,k]!=t[1,k]){l=l+1}}      ## campared every element between t and existed system [j,]. l represents the number of different elements.
if (l>0){r=r+1}                                      ##if l>1, which means system j is different with t ,r+1  r is  the number of different system.
l<-0                                                 ##l come to zero and get next comparing
}
if(r==nrow(x)){x<-rbind(x,t)}                    # if the number of different systems equal to the number of existed systems, add t into the matrix
return(x)}                                       # retrun the operated matrix.









#function: Generate N random a*b system and form a matrix. every row of the matrix repreasents a system.


get_ns<-function(a,b,n){
d<-sample(c(-1,1),a*b,replace=T)                     #get a random a*b system x
x<-matrix(d,1,a*b)
z<-0                                                 # generate a system q which is different with x
                                                     ## set z as a initial value 0 to active the loop random generator, z is the number of a different system elements
while(z==0){                                         ## if the q is same with x 
z<-0                                                 ## set z=0 
q<-matrix(sample(c(-1,1),a*b,replace=T),1,a*b)       ## get a random system q 
for(m in 1:a*b){if (x[1,m]!=q[1,m]){z=z+1            ## compare q to x, z is the number of the diiffrent elements
}}
}
x<-rbind(x,q)                                        # add q into x.
i<-2                                                 # i: the number of system is 2
while(i<n){
if(nrow(x)==i){x<-get_random(x)                      # for row i of the matrix, get a random system 
i=i+1}                                               # if the random system is same with existed system, redo row i or do i+1
}
return(x)}                                           # return x







#selects a system of x and make it into a a*b matrix.
make_m<-function(a,b,x){
h<-matrix(x,a,b)
return(h)}





#get |m|(s)=1/(LR*LC)|SUM Sx| of a system n.

m_s<-function(n){
h<-0
x<-0
for(i in 1:nrow(n)){            #get Sum Sx.
for(j in 1:ncol(n)){
h=h+n[i,j]}
}
x=1/(nrow(n)*ncol(n))*abs(h)    #get |m|(s)=1/(LR*LC)|SUM Sx|  .  
return(x)
}





#get m bar for random system.
#n is the number of the random system.a,b is the dimansion of the matrix.

m_bar<-function(a,b,n){
x<-get_ns(a,b,n)                 #get n a*b system and add it in rows of a matrix
h<-0
s<-0
s2<-0
M<-0
for (m in 1:n){                  # operate the m system. 
h<-make_m(a,b,x[m,])             # pick the m system of n systems.And make it as a matrix.
s=s+(-0.5)*get_H(h)*m_s(h)   # Sum the  numerator
s2<-s2+exp((-0.5)*get_H(h))  # Sum the denominator
}           
M=s/s2                           # m bar = numerator/denominator
return(M)                        # return M
}  



#m bar of a  stable system with n RANDOM SYSTEM.
x<-get_ns(a,b,n)               #get a random system group
m_bar_n<-function(a,b,x){
h<-0
s<-0
s2<-0
M<-0
for (m in 1:nrow(x)){          #for every system, m is the number of system
h<-make_m(a,b,x[m,])           #pick m system and make it in a matrix form.
s=s+(-0.5)*get_H(h)*m_s(h)     #calculate the sum of numerator
s2<-s2+exp((-0.5)*get_H(h))    #calculate the sum of diamenator 
}        
M=s/s2                         # m bar = numerator/denominator
return(M)                      #return M
}



#Metropolis algotithm
M_a<-function(x){
H_i<-get_H(x)                             # same initial Hamiltonian in H_i
i<-sample(1:nrow(x),1)                    # get a random number i and j
j<-sample(1:ncol(x),1)
x_i<-x                                    # save initial system in x_i
P_a<-0
s_1<-0
x[i,j]<-(-x[i,j])                         # flip a random element x[i,j] and form system x bar
H_f<-get_H(x)                             # get hamiltonian of system x bar          
d_H<-H_f-H_i                              # get hamiltonian difference d_H
if(exp(-d_H)<1){P_a<-exp(-d_H)}           # make P_a= min (i ,-exp(-d_H))
else{P<-1}
r<-runif(1,0,1)                           # get a random number r from the uniform distribution on a unit intervel.
if(r<P_a){s_1<-x}                         # if pa>r get s bar or get s
else{s_1<-x_i}                               
return(s_1)                               #turn to the transform system
}






# simple m bar calculation
m_average_simple<-function(a,b,N){
x<-get_ns(a,b,N)             # get N system
S<-0
for(m in 1:N){               # select m system
h<-make_m(a,b,x[m,])         # pick m system and make it in a matrix
S<-S+m_s(h)                  # Sum the |m|(s) of all the system
}
M<-S/N                       # m bar = Sum the |m|(s)/N
return(M)}



d<-sample(c(-1,1),10*10,replace=T)
x<-matrix(d,10,10)
p<-0
q<-0
for (i in 1:10000){
d<-x
image(x)
x<-M_a(x)
for(m in 1:10){for (n in 1:10){if(x[m,n]==d[m,n]){q=q+1}}}
if(q!=10*10){
p=p+1
}
q<-0
}

> p
[1] 66







