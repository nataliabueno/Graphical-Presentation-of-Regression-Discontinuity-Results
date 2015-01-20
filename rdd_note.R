############################################################################
# Bueno & Tuñón, 2015
# Replication code 
# Note to reader: We intend for this material to be published as an online Appendix or Supplemental Material. 
# It will also be publicly available at github by the time of publication.
############################################################################
###Estimators functions
#Difference of means 
#Polynomial Regressions
#Local Linear (with Kernel)

#Optimal bandwidth function

###Plotting function
#RD_plot

###Main text (Using Caughey and Sekhon, 2011, replication data)
#(a)Main Estimates (with number of observations and no additional estimator)
#(b)Main Estimates (with number of observations and an additional estimator)
#(c)Balance (without number of observations and without additional estimator, 
#			for 2 covariates)

#Online Appendix
#(A.1)Main Estimates (with number of observations and an additional estimator, 
#					 polynomial regression)
#(A.2)Main Estimates (with number of observations and no additional estimator, 
#					 but standardized outcome variable.)
#(A.3)Balance (without number of observations and without additional estimator, 
#			for 12 covariates)

#Additional figures (types of plot not shown in text nor appendix)
#(d)Main Estimates (without number of observations and without additional estimator)
#(e)Main Estimates (with number of observations and with an additional estimator,
#					loca linear regression as the main estimator and 
#					difference of means as an additional)
############################################################################

rm(list = ls())

#Preambule
options(scipen=999) # supressing scientific notation

#Libraries
library(foreign)
library(sandwich)

#Data directory
#User should adjust these accordingly
data.dir <- setwd("~/Dropbox/Graphical-Presentation-of-Regression-Discontinuity-Results/Caughey_Sekhon_2011")

#Estimator Functions

#Difference of mean (using OLS with robust standard errors)
dom <- function(rescaled, treat, outcome){
  model <- lm(outcome~treat)
  est <- NA
  est[1] <- model$coefficients[2]
  est[2] <- sqrt(diag(vcovHC(model,type="HC3"))[2])
  return(est)
}

#Regression with fourth-degree polynomials
polynomial <- function(rescaled, treat, outcome, deg=4){
  model <- lm(outcome ~ treat + poly(rescaled, degree=deg,raw=T) +
                poly(rescaled*treat, degree=deg,raw=T))
  est <- NA
  est[1] <- model$coefficients[2]
  est[2] <- sqrt(diag(vcovHC(model,type="HC3"))[2])
  return(est)
}

### MODIFIED CODE FROM:
### Replication code for "Testing the Accuracy of Regression Discontinuity Analysis 
### Using Experimental Benchmarks" by Donald P. Green, Terence Y. Leong, Holger L. Kern,
### Alan S. Gerber, and Christopher W. Larimer

# Local linear regression with triangular (edge) kernel

loc_lin <- function(outcome, rescaled, cutoff, running, h=windows, treat)   {
  
  require(sandwich)
  
  temp <- (running - cutoff) / h
  kern <- as.numeric(abs(temp) <= 1) * (1 - abs(temp))
  
  runningCen <- running - cutoff
  RunningCenTreat <- running
  sel <- kern > 0
  
  wls.pool <- lm(outcome ~ treat + runningCen + RunningCenTreat, weights = kern, subset = sel)
  summary(wls.pool)
  
  est <- NA
  est[1] <- summary(wls.pool)$coeff[2,1]
  est[2] <- sqrt(diag(vcovHC(wls.pool, type = "HC3")))[2]
  return(est)
}

### CODE FROM:
### Replication code for "Testing the Accuracy of Regression Discontinuity Analysis 
### Using Experimental Benchmarks" by Donald P. Green, Terence Y. Leong, Holger L. Kern,
### Alan S. Gerber, and Christopher W. Larimer

#Optimal Bandwidth
optimal.bw <- function(Y,X,c,reg)   {
  
  # STEP 1
  above <- X >= c
  N <- length(X)
  S2x <- sum((X-mean(X))^2) / (N-1)
  h1 <- 1.84 * sqrt(S2x) * N^(-1/5)
  
  Nh1.neg <- length(X[X < c & X >= (c-h1)])
  Nh1.pos <- length(X[X >= c & X <= (h1+c)])
  
  Yh1.neg <- sum(Y[X < c & X >= (c-h1)]) / Nh1.neg
  Yh1.pos <- sum(Y[X >= c & X <= (h1+c)]) / Nh1.pos
  
  fhatx <- (Nh1.neg + Nh1.pos) / (N*h1*2)
  
  fhatx.temp1 <- sum((Y[X < c & X >= (c-h1)] - Yh1.neg)^2)
  fhatx.temp2 <- sum((Y[X >= c & X < (h1+c)] - Yh1.pos)^2)
  
  sigmahat2 <- (fhatx.temp1 + fhatx.temp2) / (Nh1.neg + Nh1.pos)
  
  
  # STEP 2
  
  medianX.pos <- median(X[X >= c])
  medianX.neg <- median(X[X < c])
  
  keep <- !(X < medianX.neg | X > medianX.pos)
  
  X.cen <- X - c
  X.cen2 <- X.cen^2
  X.cen3 <- X.cen^3
  
  lm.out <- lm(Y ~ above + X.cen + X.cen2 + X.cen3, subset = keep)
  summary(lm.out)
  
  mhat3 <- summary(lm.out)$coef[5,1] * 6
  
  
  N.pos <- sum(above)
  N.neg <- length(above) - sum(above)
  
  h2.pos <- (sigmahat2 / (fhatx * max(mhat3^2, .01)))^(1/7) * 3.56 * N.pos^(-1/7)
  h2.neg <- (sigmahat2 / (fhatx * max(mhat3^2, .01)))^(1/7) * 3.56 * N.neg^(-1/7)
  
  
  Yhat.pos <- Y[X >= c & X <= (h2.pos+c)]
  Xhat.pos <- X[X >= c & X <= (h2.pos+c)]
  
  Yhat.neg <- Y[X < c & X >= (c-h2.neg)]
  Xhat.neg <- X[X < c & X >= (c-h2.neg)]
  
  N2.pos <- length(Xhat.pos)
  N2.neg <- length(Xhat.neg)
  
  Xhat.pos.cen <- Xhat.pos - c
  Xhat.pos.cen2 <- Xhat.pos.cen^2
  
  Xhat.neg.cen <- Xhat.neg - c
  Xhat.neg.cen2 <- Xhat.neg.cen^2
  
  
  lm.out <- lm(Yhat.pos ~ Xhat.pos.cen + Xhat.pos.cen2)
  summary(lm.out)
  mhat2.pos <- summary(lm.out)$coef[3,1] * 2
  
  lm.out <- lm(Yhat.neg ~ Xhat.neg.cen + Xhat.neg.cen2)
  summary(lm.out)
  mhat2.neg <- summary(lm.out)$coef[3,1] * 2
  
  
  # STEP 3
  
  rhat.pos <- (720 * sigmahat2) / (N2.pos * h2.pos^4) # NOT identical to IK paper
  rhat.neg <- (720 * sigmahat2) / (N2.neg * h2.neg^4) # NOT identical to IK paper
  
  if(reg == F)    { rhat.pos <- 0; rhat.neg <- 0 }    # turn off regularization
  
  
  hhat.opt <- ((2 * sigmahat2) /
                 (fhatx * ((mhat2.pos - mhat2.neg)^2 + (rhat.pos + rhat.neg))))^(1/5) *
    3.4375 * N^(-1/5)
  
  
  return(hhat.opt)
  
}




############################################################################
# Plotting function
RD_plot <- function(
  
  # Data
  running=running,
  treat=treat,
  outcome=outcome,
  cutoff=0,
  min_running=F,
  max_running=F,
  
  # Intervals and Opt bw
  opt_bw=F, # opt bw
  nr_windows=50,
  
  # Estimator/s
  # Main Estimator
  main_est="dom",
  # estimator
  #confidence intervals
  ci="95%", #option: ci="90%", ci="F"
  add_est=F, # no additional estimator
  
  #Plot aesthetics
  nr_obs=T,
  nr_obs_lab,
  label_x="Absolute distance from the cutoff",
  label_y="Estimated effect",
  plot_label="RD plot",
  label_size=1.2,
  main_size=2,
  legend=F # Keep false if desired output is grid of plots
  ){
  
  # WARNINGS
  if (cutoff==0) print("Assuming that cutoff is at the value of 0 for the running variable. \n 
                       To change this, input new value for the cutoff argument.")
  if (cutoff!=0) print("Rescaling the running variable to place cutoff at zero.")
  library("sandwich")
  print("This function requires the sandwich package")

  
  
  # DATA
  # 0. rescaling running var and organizing data 
  rescaled <- running - cutoff
  abs_running <- abs(rescaled)
  data <- cbind(rescaled, treat, outcome, abs_running, running)
  
  # 1. get windows from intervals and add optimal bandwidth
  min_running <- ifelse(min_running==F, 0, min_running-cutoff)
  max_running <- ifelse(max_running==F, max(abs_running), max_running-cutoff)
  windows <- seq(min_running, max_running,by=(max_running/nr_windows))[-1]
  if (opt_bw!=F) windows <- sort(c(windows, opt_bw))
  if (opt_bw!=F) print("Assuming optimal bandwidth measured as distance from the cutoff.")
  
  # 2. calculate estimate(s) and CI
  # 2.a main estimate
  if (main_est=="dom") main <- dom
  if (main_est=="poly") main <- polynomial
  if (main_est=="loc_lin") main <- loc_lin

  # 2.b additional
  if (add_est=="dom") second <- dom
  if (add_est=="poly") second <- polynomial
  if (add_est=="loc_lin") second <- loc_lin

  #2.c loop
  ests <- matrix(NA,length(windows),3)
  for (i in 1:length(windows)){
    # select data
    temp <- as.data.frame(data[abs_running<=windows[i],])
    # get estimates
    if (main_est=="loc_lin"){ 
        ests[i,1:2] <- with(temp, main(rescaled=rescaled, treat=treat,
                       outcome=outcome, running=running,
                       cutoff=0, h=windows[i]))}
    else{
        ests[i,1:2] <- with(temp, main(rescaled=rescaled, treat=treat, outcome=outcome))}
    
    if (add_est==F) next 
    if (add_est=="loc_lin"){
        ests[i,3] <- with(temp, second(rescaled=rescaled, treat=treat,
                       outcome=outcome, running=running,
                       cutoff=0, h=windows[i]))[1]}
    else{    
    ests[i,3] <- with(temp, second(rescaled=rescaled, treat=treat, outcome=outcome))[1]}
  }
  #2.d confidence intervals for main estimate
  if (ci=="95%") CI <- cbind(ests[,1]+1.96*ests[,2],ests[,1]-1.96*ests[,2])
  if (ci=="90%") CI <- cbind(ests[,1]+1.64*ests[,2],ests[,1]-1.64*ests[,2])
  # 2.b Define max and min for y axis
  if (add_est!=F){
    lim_y <- NA
    lim_y[1] <- ifelse(max(CI)>=max(ests[,3]), max(CI), max(ests[,3]))
    lim_y[2] <- ifelse(min(CI)<=min(ests[,3]), min(CI), min(ests[,3]))}
  if (add_est==F){ lim_y <- c(max(CI), min(CI))}
  lim_y[3] <- lim_y[2] - (lim_y[1]-lim_y[2])/7
  lim_y[4] <- lim_y[1] - (lim_y[1]-lim_y[2])/11
  
  # 2.c Define max for x axis
  if (abs(opt_bw) > abs(max_running)) max_axis <- opt_bw
  if (abs(opt_bw) < abs(max_running)) max_axis <- max_running
  
  # 3. ordering to know nr of observations by value of running variable
  data <- as.data.frame(data[order(abs_running),])
  if (nr_obs==T) {
    if (missing(nr_obs_lab))
      stop("Missing specification for number of observations to be included in the plot. \n
           If no axis for the number of observations is needed, use nr_obs=F.")
    nr_obs_lab <- cbind(nr_obs_lab, data$abs_running[nr_obs_lab])}
  
  
  # PLOT
  # Main line (add types, mins and maxs)
  if(legend==F){
    par(mar=c(4.1,6.1,6.1,2.1))
    plot(windows, ests[,1], pch=16, axes=T, 
         #limits
         ylim=c(lim_y[3],lim_y[1]),
         xlim=c(0, max_axis),
         #labels
         xlab=label_x, ylab=label_y, main=plot_label,
         bty='n', cex.main=main_size,
         cex.lab=label_size, 
         cex=1.2)
    #Zero effect line
    abline(h=0, lty=4, lw=3)
    #opt bw vertical line
    if(opt_bw!=F){
      abline(v=opt_bw, lty=5, col="gray80", lw=3)}
    
    # Add estimates line
    lines(windows, ests[,3], lwd=3,col="gray30")
    # CIs 
    if(ci!=F){
      lines(windows, CI[,1], lty=2,col="gray30")
      lines(windows, CI[,2], lty=2,col="gray30")}
    # Number of observations
    if (nr_obs==T) {
      axis(3, at=c(nr_obs_lab[,2]), labels=c(nr_obs_lab[,1]), cex=.6, col="grey50", 
           lwd = 0.5, padj=1, line=1, cex.axis=.7, col.axis="grey50")
      mtext("Number of observations", side=3, col="grey50", cex=.7, adj=0)}
  }
  
  if(legend==T){
  layout(rbind(1,2), heights=c(7,1))
  par(mar=c(4.1,6.1,6.1,2.1))
  plot(windows, ests[,1], pch=16, axes=T, 
       #limits
       ylim=c(lim_y[3],lim_y[1]),
       xlim=c(0, max_axis),
       #labels
       xlab=label_x, ylab=label_y, bty='n', 
       cex.lab=label_size, 
       cex=1.2)
  # title
  title(plot_label, cex.main=main_size, line=4)
  #Zero effect line
  abline(h=0, lty=4, lw=3)
  #opt bw vertical line
  if(opt_bw!=F){
  abline(v=opt_bw, lty=5, col="gray80", lw=3)}
  
  # Add estimates line
  lines(windows, ests[,3], lwd=3,col="gray30")
  # CIs 
  if(ci!=F){
    lines(windows, CI[,1], lty=2,col="gray30")
    lines(windows, CI[,2], lty=2,col="gray30")}
  # Number of observations
  if (nr_obs==T) {
    axis(3, at=c(nr_obs_lab[,2]), labels=c(nr_obs_lab[,1]), cex=.6, col="grey50", 
         lwd = 0.5, padj=1, line=1, cex.axis=.7, col.axis="grey50")
    mtext("Number of observations", side=3, col="grey50", cex=.7, adj=0)}
  
  # Legend
  par(mar=c(0.5,0,0,0))
  plot.new()
  if (main_est=="dom") main_leg <- "Difference of Means"
  if (main_est=="poly") main_leg <- "Polynomial Regression"
  if (main_est=="loc_lin") main_leg <- "Local Linear Regression"
  
  if (add_est!=F){
    if (add_est=="dom") second_leg <- "Difference of Means"
    if (add_est=="poly") second_leg <- "Polynomial Regression"
    if (add_est=="loc_lin") second_leg <- "Local Linear Regression"}
  
  if (ci=="95%") ci_leg <- paste("95% Confidence Intervals ", "(", main_leg, ")", sep="")
  if (ci=="90%") ci_leg <- paste("90% Confidence Intervals ", "(", main_leg, ")", sep="")
  
  if (add_est!=F){
    legend("center", "groups", c(main_leg, second_leg, ci_leg), pch=c(16, NA, NA), 
           lty=c(NA, 1, 2), lwd=c(2, 2, 2),
           col=c("black", "gray", "gray"), bty="n", cex=0.8)}
  
  if (add_est==F){
    legend("center", "groups", c(main_leg, ci_leg), pch=c(16, NA), lty=c(NA, 2), lwd=c(2, 2),
           col=c("black", "gray"), bty="n", cex=0.8)}
  }
   
}

############################################################# Analysis

#Getting data
rd <- read.dta("~/Dropbox/Graphical-Presentation-of-Regression-Discontinuity-Results/Caughey_Sekhon_2011/RDReplication.dta")
data <- rd[rd$Use == 1, ]

#Getting optimal bandwidth following imbens and Kalynaranam (2009). Using function from Green et al 2009.
data <- data[!is.na(data$DPctNxt),]
#using only non-missing cases.  
opt_bw <- optimal.bw(Y=data$DPctNxt, X=data$DifDPct, c=0, reg=T) #-1590 observation from DPctNxt

#Recoding pre-treatment covariates for balance plots
data$DWinPrv_r <- (data$DWinPrv - mean(data$DWinPrv, na.rm=T))/sd(data$DWinPrv, na.rm=T)
data$DPctPrv_r <- (data$DPctPrv - mean(data$DPctPrv, na.rm=T))/sd(data$DPctPrv, na.rm=T)
data$DemInc_r <- (data$DemInc - mean(data$DemInc, na.rm=T))/sd(data$DemInc, na.rm=T)
data$NonDInc_r <- (data$NonDInc - mean(data$NonDInc, na.rm=T))/sd(data$NonDInc, na.rm=T)
data$PrvTrmsD_r <- (data$PrvTrmsD - mean(data$PrvTrmsD, na.rm=T))/sd(data$PrvTrmsD, na.rm=T)
data$PrvTrmsO_r <- (data$PrvTrmsO - mean(data$PrvTrmsO, na.rm=T))/sd(data$PrvTrmsO, na.rm=T)
data$DExpAdv_r <- (data$DExpAdv - mean(data$DExpAdv, na.rm=T))/sd(data$DExpAdv, na.rm=T) 
data$ElcSwing_r <- (data$ElcSwing - mean(data$ElcSwing, na.rm=T))/sd(data$ElcSwing, na.rm=T)
data$DSpndPct_r <- (data$DSpndPct - mean(data$DSpndPct, na.rm=T))/sd(data$DSpndPct, na.rm=T) 
data$DDonaPct_r <- (data$DDonaPct - mean(data$DDonaPct, na.rm=T))/sd(data$DDonaPct, na.rm=T) 
data$OpenSeat_r <- (data$OpenSeat - mean(data$OpenSeat, na.rm=T))/sd(data$OpenSeat, na.rm=T) 
data$VtTotPct_r <- (data$VtTotPct - mean(data$VtTotPct, na.rm=T))/sd(data$VtTotPct, na.rm=T)

#Recoding outcome variable
data$DPctNxt_r <- (data$DPctNxt - mean(data$DPctNxt, na.rm=T))/sd(data$DPctNxt, na.rm=T)

#Main text - Figures

#(a)Main Estimates (with number of observations and without additional estimator)
pdf("~/Dropbox/Graphical-Presentation-of-Regression-Discontinuity-Results/paper_latex/no_add_est_plot.pdf", width=8, height=8)
RD_plot(running=data$DifDPct, treat=data$DemWin, 
        outcome=data$DPctNxt, 
        cutoff=0, min_running=0.25, max_running=10,
        opt_bw=opt_bw, nr_windows=50, main_est="dom", ci="95%",
        add_est=F, nr_obs=T, 
        nr_obs_lab=c(30, 250, 500, 1000, 1250, 1500),
        label_x="Absolute Distance from the Cutoff (Vote Margin %)",
        label_y="Average Difference in Vote Share \n Between Bare Losers and Bare Winners",
        plot_label=" ",
        legend=T)
dev.off()

#(b)Main Estimates (with number of observations and additional estimator - local linear)
pdf("~/Dropbox/Graphical-Presentation-of-Regression-Discontinuity-Results/paper_latex/with_add_est_plot.pdf", width=8, height=8)
RD_plot(running=data$DifDPct, treat=data$DemWin, 
        outcome=data$DPctNxt, 
        cutoff=0, min_running=0.25, max_running=10,
        opt_bw=opt_bw, nr_windows=50, main_est="dom", ci="95%",
        add_est="loc_lin", nr_obs=T, 
        nr_obs_lab=c(30, 250, 500, 1000, 1250, 1500),
        label_x="Absolute Distance from the Cutoff",
        label_y="Estimated Effect",
        plot_label=" ",
        legend=T)
dev.off()


#(c)Balance (without number of observations and without additional estimator, 2 covariates)

#(c.1)Balance - Dem Win t-1
pdf("~/Dropbox/Graphical-Presentation-of-Regression-Discontinuity-Results/paper_latex/balance_plot_1.pdf", height=5, width=6.5)
par(mfrow=c(1,1))
RD_plot(running=data$DifDPct, treat=data$DemWin, 
        outcome=data$DWinPrv_r, 
        cutoff=0, min_running=0.25, max_running=10,
        opt_bw=opt_bw, nr_windows=50, main_est="dom", ci="95%",
        add_est=F, nr_obs=F, 
        nr_obs_lab=NA,
        label_x="Absolute Distance from the Cutoff (Vote margin %)",
        label_y="Difference between Treatment \n and Control",
        label_size=1.5,
        plot_label=" ",
        legend=F) 
dev.off()

#(c.2)Balance - Voter turnout %
pdf("~/Dropbox/Graphical-Presentation-of-Regression-Discontinuity-Results/paper_latex/balance_plot_2.pdf", height=5, width=6.5)
par(mfrow=c(1,1))
# Voter turnout %
RD_plot(running=data$DifDPct, treat=data$DemWin, 
        outcome=data$VtTotPct_r, 
        cutoff=0, min_running=0.25, max_running=10,
        opt_bw=opt_bw, nr_windows=50, main_est="dom", ci="95%",
        add_est=F, nr_obs=F, 
        nr_obs_lab=NA,
        label_x="Absolute Distance from the Cutoff (Vote margin %)",
        label_y="Difference between Treatment \n and Control",
        label_size=1.5,
        plot_label=" ",
        legend=F) 
dev.off()


#Online Appendix - Figures

#(A.1)Main Estimates (with number of observations and additional estimator - polynomial)
pdf("~/Dropbox/Graphical-Presentation-of-Regression-Discontinuity-Results/paper_latex/with_add_est_plot_poly.pdf", width=8, height=8)
RD_plot(running=data$DifDPct, treat=data$DemWin, 
        outcome=data$DPctNxt, 
        cutoff=0, min_running=0.25, max_running=max(data$DifDPct),
        opt_bw=opt_bw, nr_windows=50, main_est="dom", ci="95%",
        add_est="poly", nr_obs=T, 
        nr_obs_lab=c(30, 1000, 25000, 5000, 7000, 8500),
        label_x="Absolute Distance from the Cutoff (Vote Margin %)",
        label_y="Estimated Effect",
        plot_label=" ",
        legend=T)
dev.off()

#(A.2)Main Estimates (with number of observations and without additional estimator): standardized outcome variable 

pdf("~/Dropbox/Graphical-Presentation-of-Regression-Discontinuity-Results/paper_latex/no_add_est_plot_standardized.pdf", width=8, height=8)
RD_plot(running=data$DifDPct, treat=data$DemWin, 
        outcome=data$DPctNxt_r, 
        cutoff=0, min_running=0.25, max_running=10,
        opt_bw=opt_bw, nr_windows=50, main_est="dom", ci="95%",
        add_est=F, nr_obs=T, 
        nr_obs_lab=c(30, 250, 500, 1000, 1250, 1500),
        label_x="Absolute Distance from the Cutoff (Vote Margin %)",
        label_y="Average Difference in Vote Share \n Between Bare Losers and Bare Winners",
        plot_label=" ",
        legend=T)
dev.off()

#(A.3)Balance (without number of observations and without additional estimator, 12 covariates)

outcomes   <-  list(data$DWinPrv_r, data$DPctPrv_r,  
                   data$DemInc_r, data$NonDInc_r, data$PrvTrmsD_r, 
                   data$PrvTrmsO_r, data$DExpAdv_r,
                   data$ElcSwing_r, data$DSpndPct_r, data$DDonaPct_r, 
                   data$OpenSeat_r, data$VtTotPct_r)


names <- c('Dem Win t - 1', 'Dem % t - 1',
            'Dem % Margin t - 1',
            'Dem Inc in Race', 'Rep Inc in Race',
            'Dem\'s # Prev Terms', 'Rep\'s # Prev Terms',  'Dem Experience Adv',
            'Partisan Swing', 'Dem Donation %',
            'Open Seat', 'Voter Turnout %')


pdf("~/Dropbox/Graphical-Presentation-of-Regression-Discontinuity-Results/paper_latex/balance_plots_12.pdf", height=12, width=22)
par(mfrow=c(4,3))
for (i in 1:length(outcomes)){
RD_plot(running=data$DifDPct, treat=data$DemWin, 
          outcome=outcomes[[i]], 
          cutoff=0, min_running=0.25, max_running=10,
          opt_bw=opt_bw, nr_windows=50, main_est="dom", ci="95%",
          add_est=F, nr_obs=F, 
          nr_obs_lab=NA,
          label_x="Abs. Distance from the Cutoff",
          label_y="Difference between Treatment \n and Control",
          plot_label=names[i],
          label_size=1.8,
          main_size=2.3,
          legend=F) 
}
dev.off()

#Additional Figures

#(d)Main Estimates (without number of observations and without additional estimator)
pdf("~/Dropbox/Graphical-Presentation-of-Regression-Discontinuity-Results/paper_latex/no_add_est_no_obs_plot.pdf", width=8, height=8)
RD_plot(running=data$DifDPct, treat=data$DemWin, 
        outcome=data$DPctNxt, 
        cutoff=0, min_running=0.25, max_running=10,
        opt_bw=opt_bw, nr_windows=50, main_est="dom", ci="95%",
        add_est=F, nr_obs=F, 
        nr_obs_lab=NA,
        label_x="Absolute Distance from the Cutoff (Vote margin %)",
        label_y="Average Difference in Vote Share \n Between Bare Losers and Bare Winners",
        plot_label="Effect of Incumbency on Vote Share",
        legend=T)
dev.off()

#(e)Main Estimates (with number of observations and with an additional estimator,
#					loca linear regression as the main estimator and 
#					difference of means as an additional)

pdf("~/Dropbox/Graphical-Presentation-of-Regression-Discontinuity-Results/paper_latex/no_add_est_no_obs_plot.pdf", width=8, height=8)
RD_plot(running=data$DifDPct, treat=data$DemWin, 
        outcome=data$DPctNxt, 
        cutoff=0, min_running=0.25, max_running=10,
        opt_bw=opt_bw, nr_windows=50, main_est="loc_lin", ci="95%",
        add_est="dom", nr_obs=F, 
        nr_obs_lab=NA,
        label_x="Absolute Distance from the Cutoff (Vote margin %)",
        label_y="Estimated Effect",
        plot_label="Effect of Incumbency on Vote Share",
        legend=T)
dev.off()
