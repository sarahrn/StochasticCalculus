euler<-function(T,t,N,r,X0,Y0,a,xi,rho,b){
    #initialize variables
    deltat<-(T-t)/N
    
    #construct empty vectors X, Y.
    Y<-rep(0,(N+1))
    X<-rep(0,(N+1))
    eta<-rnorm(N+1)
    nu<-rnorm(N+1)
    
    #initialize first entry in the vector to X0
    Y[1]<-Y0
    X[1]<-X0
    
    #use Euler discretization to calculate the next 20 values of X,Y
    for(k in 2:(N+1)){
        Y[k]<-Y[k-1]+a*(b-Y[k-1])*deltat+xi*sqrt(deltat)*(rho*eta[k-1]+sqrt(1-rho^2)*nu[k-1])
        X[k]<-X[k-1]*(1+r*deltat+exp(Y[k-1])*sqrt(deltat)*eta[k-1])
    }
    
    return(data.frame(cbind(X,Y)))
}

getvols<-function(T,t,N,r,X0,Y0,a,xi,rho,b){
    meanprice<-c()
    price<-c()
    vols<-c()
    for(K in seq(90,110,by=1)){
        print(K)
        for(n in 1:1000){
            XY<-euler(T,t,N,r,X0,Y0,a,xi,rho,b)
            price<-append(price,exp(-r*(T-t))*max((XY$X[N+1]-K),0))
        }
        vols<-append(vols,impvol(T,t,X0,K,r,mean(price),.01))
        meanprice<-append(meanprice,mean(price))
    }
    return(data.frame(K=seq(90,110,by=1),vols,meanprice))
}

impvol<-function(T,t,S,K,r,C,alpha){
    Ch<-0;sigma<-0
    while(abs(Ch-C)>alpha){
        sigma<-sigma+alpha/100
        d1<-(log(S/K)+(r+sigma^2/2)*(T-t))/(sigma*sqrt(T-t))
        d2<-d1-sigma*sqrt(T-t)
        Ch<-pnorm(d1)*S-pnorm(d2)*K*exp(-r*(T-t))
        #print(paste(Ch," ",C))
    }
    return(sigma)
}

vasicek<-function(r0,a,b,sigma,lambda,T){
    #declare variables
    rbar=0; A=0; G=0; B=0
    
    rbar<-b+(sigma*lambda)/a-sigma^2/(2*a^2)
    G<-(1-exp(-a*T))/a
    A<-exp(rbar*(G-T)-(G^2*sigma^2)/(4*a))
    B<-A*exp(-G*r0)
    return(B)
}

cir<-function(r0,a,b,sigma,lambda,T){
    #declare variables
    rbar<-0; A<-0; G<-0; B<-0; gamma<-0;sigma<-sigma*sqrt(r0)
    
    B<-A*exp(-G*r0)
    gamma<-sqrt(a^2+2*sigma^2)
    A<-((2*gamma*exp(((a+gamma)*T)/2))/((gamma+a)*(exp(gamma*T)-1)+2*gamma))^((2*a*b)/(sigma^2))
    G<-(2*(exp(gamma*T)-1))/((gamma+a)*(exp(gamma*T)-1)+2*gamma)
    B<-A*exp(-G*r0)
    return(B)
}

run<-function(){
    binom(.15,1,1000,.03,100,100)
    explicit(exp(10),exp(10),1,5000,100,.2,9,11,0)
    implicit(exp(10),exp(10),1,100,2000,.2,9,11,0)
    geoblackscholes(380,376,.04,1,0,.4)
    arithmetic(380,376,.04,1,0,.4,1000,1000)
    geometric(380,376,.04,1,0,.4,1000,1000)
    bstar(380,376,.04,1,0,.4,100,100,1000)
    controlvariate(380,376,.04,1,0,.4,100,100,1000)
    findk(380,376,.04,1,0,.4,1000,1000,1000)
    plot(380,376,.04,1,0,.4,10,10)
}

blackscholes<-function(K,S,r,T,t,sigma){
    d1<-(log(S/K)+(r+sigma^2/2)*(T-t))/(sigma*sqrt(T-t))
    d2<-d1-sigma*sqrt(T-t)
    return(pnorm(d1)*S-pnorm(d2)*K*exp(-r*(T-t)))
}
binom<-function(Sigma, t, n, R, S0, K){
    #calculate u, d, p, q, discount factor
    u = exp(Sigma*sqrt(t/n))
    d = exp(-Sigma*sqrt(t/n))
    p = (exp((R*t)/n)-d)/(u-d)
    q = 1-p
    disc = exp(-R*(t/n))
    
    #declare matrices for stock price and option price over time
    S<-diag(0,n)
    C<-diag(0,n)
    
    #initialize first row of S
    for(j in 0:(n-1)){
        S[1,j+1]=S0*u^j
    }
    
    #get remaining rows of S
    for (i in 2:n){
        b = i
        for (j in seq(b,n,by=1)){
            S[i,j]=S[i-1,j-1]*d;
        }
    }
    
    #initialize final column of C
    for (i in 0:(n-1)){
        C[i+1,n]=max((K-S[i+1,n]),0.0);
        #    cout << S[i][n-1] << endl;
    }
    
    #get remaining columns of C
    for (j in seq(n-1,1,by=-1)){
        b = j
        for (i in 1:b){
            C[i,j]=max((disc*(p*C[i,j+1]+q*C[i+1,j+1])),(K-S[i,j]))
        }
    }
    
    #return the put price
    return(C[1,1])
}
h<-function(S,T,t,r){
    return(S*exp((r+(1/2)*sigma^2)*(T-t)+sigma*sqrt(T-t)*(rnorm(1))))
}
explicit<-function(S,K,t,n,m,sigma,xl,xu,r){
    #define variables
    deltat=t/n
    deltax=(xu-xl)/m
    alpha=deltat/deltax^2
    
    #define matrices
    A<-diag(0,m-1)
    U<-matrix(nrow=(n+1),ncol=(m+1),0)
    u<-c()
    
    #set A
    for (i in 1:(m-1)){
        A[i,i]=(1-alpha*sigma^2)
    }
    for (j in 1:(m-2)){
        A[j+1,j]=(alpha*sigma^2)/2
        A[j,j+1]=(alpha*sigma^2)/2
    }
    
    #boundary conditions for U
    for (q in 1:(n+1)){
        U[q, 1]=0
        U[q, m+1]=max((exp(xu)-K),0.0)
    }
    for (v in 1:m){
        U[1,v+1]=max(((exp(xl+v*deltax)-K)),0.0)
    }
    
    #set u, b
    u<-U[1,2:m]
    b<-.5*alpha*sigma^2*U[1,(m+1)]
    
    #calculate the rest of U
    for(i in 2:(n+1)){
        U[i,2:m]<-A %*% u
        U[i,m]<-U[i,m]+b
        u<-U[i,2:m]
    }
    return(exp(-r*t)*U[n+1,round(log(S)-xl)/deltax])
    #return(U[n+1, ][which(abs((U[n+1, ]-blackscholes(K,S,r,t,0,sigma)))<.5)])
}
implicit<-function(S,K,t,n,m,sigma,xl,xu,r){
    #define variables
    deltat=t/n
    deltax=(xu-xl)/m
    alpha=deltat/deltax^2
    
    #define matrices
    A<-diag(0,m-1)
    U<-matrix(nrow=(n+1),ncol=(m+1),0)
    u<-c()
    
    #set A
    for (i in 1:(m-1)){
        A[i,i]=(1+alpha*sigma^2)
    }
    for (j in 1:(m-2)){
        A[j+1,j]=-(alpha*sigma^2)/2
        A[j,j+1]=-(alpha*sigma^2)/2
    }
    
    #invert A
    A<-solve(A)
    
    #boundary conditions for U
    for (q in 1:(n+1)){
        U[q, 1]=0
        U[q, m+1]=max(K-(exp(xu)),0.0)
    }
    for (v in 1:m){
        U[1,v+1]=max((K-(exp(xl+v*deltax))),0.0)
    }
    
    #set u, b
    u<-U[1,2:m]
    b<-.5*alpha*sigma^2*U[1,(m+1)]
    u[m-1]<-u[m-1]+b
    
    #calculate the rest of U
    for(i in 2:(n+1)){
        U[i,2:m]<-A %*% u
        u<-U[i,2:m]
        u[m-2]<-u[m-2]+b
    }
    return(exp(-r*t)*U[n+1,round((log(S)+(r-.5*sigma^2)*t-xl)/deltax)])
    #ceiling(1000*(deltax+deltat/10))
}

geoblackscholes<-function(K,S,r,T,t,sigma){
    delta=.5*(r-1/6*sigma^2)
    sigma=sigma/sqrt(3)
    d1<-(log(S/K)+(r-delta+.5*sigma^2)*(T-t))/(sigma*sqrt(T-t))
    d2<-d1-sigma*sqrt(T-t)
    C<-pnorm(d1)*S-pnorm(d2)*K*exp(-r*(T-t))
    return(C)
}

arithmetic<-function(K,S0,r,T,t,sigma,n,N){
    #initialize variables
    h=(T-t)/n; n_av<-c(); N_av<-0;
    
    #use discretization of prompt function to get
    #subsequent values of S
    #Put the average of every S_t1,j...S_tn,j, less the call,
    #and maximized versus zero, into n_av
    for(i in 1:N){
        S<-rep(0,(n+1))
        S[1]<-S0
        for(k in 2:(n+1)){
            S[k]<-S[k-1]*exp((r-(1/2)*sigma^2)*h+sigma*sqrt(h)*(rnorm(1)))
        }
        #print(mean(S[1:(n+1)]))
        n_av<-append(n_av,max((mean(S[2:(n+1)])-K),0))
        #print(n_av)
    }
    #print(n_av)
    
    #average n_av and store as N_av
    N_av<-mean(n_av)
    
    return(exp(-r*T)*N_av)
}

geometric<-function(K,S0,r,T,t,sigma,n,N){
    #initialize variables
    h=(T-t)/n; n_av<-c(); N_av<-0;
    
    #use discretization of prompt function to get
    #subsequent values of S
    #Put the average of every S_t1,j...S_tn,j, less the call,
    #and maximized versus zero, into n_av
    for(i in 1:N){
        S<-rep(0,(n+1))
        S[1]<-S0
        for(k in 2:(n+1)){
            S[k]<-S[k-1]*exp((r-(1/2)*sigma^2)*h+sigma*sqrt(h)*(rnorm(1)))
        }
        #print(mean(S[1:(n+1)]))
        n_av<-append(n_av,max((exp(mean(log(S[2:(n+1)])))-K),0))
        #print(n_av)
    }
    #print(n_av)
    
    #average n_av and store as N_av
    N_av<-mean(n_av)
    
    return(exp(-r*T)*N_av)
}

bstar<-function(K,S0,r,T,t,sigma,n,N,M){
    #calculate hte bstar used in control variate
    #by approximating the covariance btw the arithmetic
    #and geometric prices
    n_ar<-c()
    n_geo<-c()
    for(i in 1:M){
        n_ar<-append(n_ar,arithmetic(K,S0,r,T,t,sigma,n,N))
        n_geo<-append(n_geo,geometric(K,S0,r,T,t,sigma,n,N))
    }
    return(cov(n_ar,n_geo)/var(n_geo))
}

controlvariate<-function(K,S0,r,T,t,sigma,n,N,M){
    #define arrays for arithmetic and geometric functions
    n_ar<-c()
    n_geo<-c()
    #simulate various arithmetic and geometric functions
    #and store in arrays
    for(i in 1:M){
        n_ar<-append(n_ar,arithmetic(K,S0,r,T,t,sigma,n,N))
        n_geo<-append(n_geo,geometric(K,S0,r,T,t,sigma,n,N))
    }
    return(mean(n_ar)+bstar(K,S0,r,T,t,sigma,n,N,M)*(mean(n_geo)-geoblackscholes(K,S0,r,T,t,sigma)))
}

findk<-function(K,S0,r,T,t,sigma,n,N,M){
    #define an approximate k that is yet unknown
    kapprox<-c()
    #get control variate price
    price<-controlvariate(K,S0,r,T,t,sigma,n,N,M)
    #iterate over several k values and mark the right ones to k
    for(k in seq(0,(3*S0),by=1)){
        if(abs(blackscholes(k,S0,r,T,t,sigma)-price)<.1){
            kapprox<-append(kapprox,k)
        }
    }
    #return mean of all the k's close to K
    return(mean(kapprox))
}

