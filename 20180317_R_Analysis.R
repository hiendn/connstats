library(mixdist)

RATES <- c(0.01, 0.05, 0.1, 0.2)

SPEC <- SENS <- FDR <- array(0,c(3,3,4))
SPEC_B <- SENS_B <- FDR_B <- array(0,c(3,3,4))

for (rr in 1:4) {
  for (ss in 1:100) {
    for (tt in 1:3) {
      for (pp in 1:3) {
        
        Rate <- RATES[rr]
        
        pvalues <- as.numeric(read.csv(paste('pvalues',tt,pp,ss,'.csv',sep='_'),header=F)[[1]])
        mask <- as.numeric(read.csv(paste('mask',tt,pp,ss,'.csv',sep='_'),header=F)[[1]])
        
        ZVAL_M <- qnorm(1-pvalues)
        # hist(ZVAL_M,xlab='z-score',probability=T,main='')
        
        HISTDAT <- mixgroup(ZVAL_M,breaks=hist(ZVAL_M,breaks='Sturges')$breaks)
        HISTDAT[1,2] <- length(which(ZVAL_M<=HISTDAT[1,1]))
        HISTDAT[dim(HISTDAT)[1],2] <- HISTDAT[dim(HISTDAT)[1],2]+sum(ZVAL_M==Inf)
        
        if(ss>1) {plot(MMIX)}
        
        # %% Initialize Mixing Probabilities
        # PI_A = 1-RATE
        # PI_B = RATE;
        # %% Initialize Means
        # MEAN_A = 0;
        # MEAN_B = mean(Z_VALUE)/(1-PI_A);
        # %% Initialize Standard Deviations
        # SD_A = 1;
        # SD_B = sqrt((std(Z_VALUE)^2-PI_A-PI_A*(1-PI_A)*MEAN_A^2)/(1-PI_A));
        PI_0 <- 1-Rate
        PI_1 <- Rate
        Mean_0 <- 0
        Mean_1 <- 3*sd(ZVAL_M[which(ZVAL_M<Inf)])
        SD_0 <- 1
        SD_1 <- 1
        PARA <- mixparam(c(Mean_0,Mean_1),
                         c(SD_0,SD_1),
                         c(PI_0,PI_1))
        MMIX <- mix(HISTDAT,PARA,
                    emsteps = 100,iterlim=10,steptol=1e-3)
        
        N <- length(ZVAL_M)
        OUTPUT <- unlist(c(MMIX$parameters[1,c(1:3)],MMIX$parameters[2,2:3]))
        TAU <- OUTPUT[1]*dnorm(ZVAL_M,OUTPUT[2],OUTPUT[3])/
          (OUTPUT[1]*dnorm(ZVAL_M,OUTPUT[2],OUTPUT[3])+(1-OUTPUT[1])*dnorm(ZVAL_M,OUTPUT[4],OUTPUT[5]))
        TAU[which(is.na(TAU)&(ZVAL_M==-Inf))] <- 1
        TAU[which(is.na(TAU)&(ZVAL_M==Inf))] <- 0
        CUMTAU <- cumsum(sort(TAU))/(1:N)
        REJECT <- rep(0,N)
        REJECT[which(TAU<=CUMTAU[which(CUMTAU>=Rate)[1]-1])] <- 1
        
        
        ADJ <- p.adjust(pvalues,method = 'fdr')
        ADJ <- ADJ<=Rate
        
        SPEC[tt,pp,rr] <- SPEC[tt,pp,rr] + sum(REJECT*mask)/sum(mask)
        SENS[tt,pp,rr] <- SENS[tt,pp,rr] + 1 - sum(REJECT*(1-mask))/sum(1-mask)
        FDR[tt,pp,rr] <- FDR[tt,pp,rr] + sum(REJECT*(1-mask))/max(sum(REJECT),1)
        
        SPEC_B[tt,pp,rr] <- SPEC_B[tt,pp,rr] + sum(ADJ*mask)/sum(mask)
        SENS_B[tt,pp,rr] <- SENS_B[tt,pp,rr] + 1 - sum(ADJ*(1-mask))/sum(1-mask)
        FDR_B[tt,pp,rr] <- FDR_B[tt,pp,rr] + sum(ADJ*(1-mask))/max(sum(ADJ),1)
        
        print(c(ss,tt,pp,rr))
        print(cbind(SPEC[,,rr],SPEC_B[,,rr])/ss)
        print(cbind(SENS[,,rr],SENS_B[,,rr])/ss)
        print(cbind(FDR[,,rr],FDR_B[,,rr])/ss)
        
      }
    }
  } 
}


