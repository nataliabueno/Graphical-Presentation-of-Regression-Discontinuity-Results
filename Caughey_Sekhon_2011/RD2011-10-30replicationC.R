rm(list = ls())
#####################
#### Directories ####
#####################
## User should adjust these accordingly
plot.dir <- "."
data.dir <- "."

###################
#### Libraries ####
###################
library(MASS)
library(foreign)
library(Matching)
library(xtable)
library(snow)
library(stats)
#library(arm)
library(car)
library(np)
library(KernSmooth)
library(boot)
library(sandwich)
library(coda)
library(rbounds)
library(gdata)
library(gtools)
library(gmodels)
library(gplots)
library(coin)
library(gsubfn)
library(mice)
library(fields)
library(pscl)
library(MCMCpack)
library(Synth)  
library(beanplot)
library(Kendall)
library(ggplot2)

### Function for naming plot files
mypdf <- function(name, ...) {
  pdf(paste(name, Sys.Date(), ".pdf", sep = ""), ...)
}

### Load data from data directory
setwd(data.dir)
rd <- read.dta("./RDReplication.dta")
## Includes correction to DifDPct and political experience variables
## in LA-3-2004

#########################################################
#### Fig. 1: Incumbent histogram broken out by party ####
#########################################################
### Define subsets of data
use <- rd$Use == 1
close.25 <- abs(rd$DifDPct) < .25 & !is.na(rd$DifDPct)
close.5 <- abs(rd$DifDPct) < .5 & !is.na(rd$DifDPct)
close1 <- abs(rd$DifDPct) < 1 & !is.na(rd$DifDPct)
close2 <- abs(rd$DifDPct) < 2 & !is.na(rd$DifDPct)
close5 <- abs(rd$DifDPct) < 5 & !is.na(rd$DifDPct)
close10 <- abs(rd$DifDPct) < 10 & !is.na(rd$DifDPct)
dinc <- rd$DWinPrv == 1 & !is.na(rd$DWinPrv)
rinc <- rd$DWinPrv == 0 & !is.na(rd$DWinPrv)
use.5 <- use & close.5

sum(use & rinc & abs(rd$DifDPct) < .5)
sum(use & dinc & abs(rd$DifDPct) < .5)
summary(rd$DWinPrv[use])
summary(rd$DWinNxt[use])
sum(use & abs(rd$DifDPct) < .5, na.rm = T)
with(rd[close.5 & use, ], table(x = DemWin, y = DWinPrv))
with(rd[close.5, ], table(x = DemWin, y = DWinPrv))

### Plot histograms

## Incumbent margin
## DMC 2011-11-30: The sample of elections used in the following
## figure differs slightly from the sample used in the figure
## `IncHistByPartyUse', resulting in small discrepancies between the
## figures. To create a version of the figure that is exactly
## compatible with the by-party figure, uncomment the following three
## lines of code:
#inc.margin <- NA
#inc.margin[dinc] <- rd$DifDPct[dinc]
#inc.margin[rinc] <- -rd$DifDPct[rinc]
## and comment the following line:
inc.margin <- ifelse(dinc, rd$DifDPct, -rd$DifDPct)

brks <- seq(-10, 10, 0.5)
setwd(plot.dir)
pdf(file = paste("IncHistUse", Sys.Date(), ".pdf", sep = ""),
    height = 3.5, width = 6.5)
par(mar = c(4.5, 4.5, 1.5, 0), mfrow = c(1, 1), cex = .6)
hist(x = inc.margin[close10 & use],
     breaks = brks,
     freq = TRUE,
     labels = TRUE,
     right = FALSE,
     xlim = c(-10, 10),
     col = "grey",
     main = "",
     xlab = "Incumbent Party Margin (%)",
     ylab = 'Frequency Count in 0.5% Bins',
     las = 1,
     cex.lab = 1.2,
     cex.main = 1.3)
abline(v = 0, lwd = 2) 
dev.off()

brks <- seq(-5, 5, 0.5)
setwd(plot.dir)
pdf(file = paste("IncHistByParty", Sys.Date(), ".pdf", sep = ""),
    height = 7, width = 7)
par(mar = c(2.5, 4.5, 1.5, 0), mfrow = c(2, 1), cex = .6)
hist(x = rd$DifDPct[close5 & dinc],
     breaks = brks,
     freq = TRUE,
     labels = TRUE,
     right = FALSE,
     xlim = c(-10, 10),
     col = "grey",
     main = "Democrat-Held Seats",
     xlab = '',
     ylab = 'Frequency Count in 0.5% Bins',
     las = 1,
     cex.lab = 1.2,
     cex.main = 1.3)
abline(v = 0, lwd = 2)
par(mar = c(4.5, 4.5, 1.5, 0))
hist(x = rd$DifDPct[close5 & rinc],
     breaks = brks,
     freq = TRUE,
     labels = TRUE,
     right = FALSE, ## excludes DifDPct = 0.5 from the closest bin,
     xlim = c(-10, 10),
     col = "grey",
     main = "Republican-Held Seats",
     xlab = 'Democratic Margin (%)',
     ylab = 'Frequency Count in 0.5% Bins',
     las = 1,
     cex.lab = 1.2,
     cex.main = 1.3)
abline(v = 0, lwd = 2)
dev.off()   

brks <- seq(-5, 5, 0.5)
pdf(file = paste("IncHistByPartyUse", Sys.Date(), ".pdf", sep = ""),
    height = 4, width = 7)
par(mar = c(4.5, 4.5, 2, 0), mfrow = c(1, 2), cex = .6)
hist(x = rd$DifDPct[close5 & use & dinc],
     breaks = brks,
     freq = TRUE,
     labels = TRUE,
     right = FALSE,
     xlim = c(-5, 5),
     ylim = c(0, 45),
     col = "grey",
     main = "Democrat-Held Seats",
     xlab = "Democratic Margin (%)",
     ylab = 'Frequency Count in 0.5% Bins',
     las = 1,
     cex.lab = 1.2,
     cex.main = 1.3)
abline(v = 0, lwd = 2)
par(mar = c(4.5, 3.5, 2, 1))
hist(x = rd$DifDPct[close5 & use & rinc],
     breaks = brks,
     freq = TRUE,
     labels = TRUE,
     right = FALSE, ## excludes DifDPct = 0.5 from the closest bin
     xlim = c(-5, 5),
     ylim = c(0, 45),
     col = "grey",
     main = "Republican-Held Seats",
     xlab = 'Democratic Margin (%)',
     ylab = '',
     las = 1,
     cex.lab = 1.2,
     cex.main = 1.3)
abline(v = 0, lwd = 2)
dev.off()

####################################################################
#### Table 2: Cross-tabulation of current and lagged Democratic ####
#### victory, 1944–2006                                         ####
####################################################################
with(rd[use & close2, ], table(x = DemWin, y = DWinPrv))
with(rd[use & close1, ], table(x = DemWin, y = DWinPrv))
with(rd[use & close.5, ], table(x = DemWin, y = DWinPrv))
sum(use & close.5)
with(rd[use & close.5, ],
     sum(DWinPrv == 1 & DemWin == 1) / sum(DemWin == 1))
with(rd[use & close.25, ], table(x = DemWin, y = DWinPrv))

######################################
#### CQ Rating in close elections ####
######################################
with(rd[use & close.5, ], table(x = DemWin, y = CQRating3))
toss <- rd$CQRating3 == 0
with(rd[use & close.5 & toss, ], table(x = DemWin, y = DWinPrv))
t.test(DWinPrv ~ DemWin, rd[use & close.5 & toss, ])
t.test(DWinNxt ~ DemWin, rd[use & close.5 & toss, ])

######################################################################
#### Table 3: Cross-tabulation of CQ Rating and Democratic        ####
#### Victory in elections decided by less than 0.5%.              ####
######################################################################
with(rd[use & close.5, ], cbind(by(DWinPrv, CQRating3, length),
                                by(DWinPrv, CQRating3, mean),
                                by(DemWin, CQRating3, mean),
                                by(DWinNxt, CQRating3, mean)))

######################################################################
#### Figure 2: Balance Plot                                       ####
######################################################################
## Outcome variables (post-treatment)
dvs <- matrix(c('DWinNxt', 'Dem Win t + 1',
                'DPctNxt', 'Dem % t + 1',
                'DifDPNxt', 'Dem % Margin t + 1'),
              ncol = 2, byrow = TRUE)
## Pre-treatment covariates
covs <- matrix(c('DWinPrv', 'Dem Win t - 1',
                 'DPctPrv', 'Dem % t - 1',
                 'DifDPPrv', 'Dem % Margin t - 1',
                 'IncDWNOM1', 'Inc\'s D1 NOMINATE',
                 'DemInc', 'Dem Inc in Race',
                 'NonDInc', 'Rep Inc in Race',
                 'PrvTrmsD', 'Dem\'s # Prev Terms',
                 'PrvTrmsO', 'Rep\'s # Prev Terms', 
                 'RExpAdv', 'Rep Experience Adv',
                 'DExpAdv', 'Dem Experience Adv',
                 'ElcSwing', 'Partisan Swing',
                 'CQRating3', 'CQ Rating {-1, 0, 1}',
                 'DSpndPct', 'Dem Spending %',
                 'DDonaPct', 'Dem Donation %',
                 'SoSDem', 'Dem Sec of State',
                 'GovDem', 'Dem Governor',
                 'DifPVDec', 'Dem Pres % Margin', ## average over decade
                 'DemOpen', 'Dem-held Open Seat',
                 'NonDOpen', 'Rep-held Open Seat',
                 'OpenSeat', 'Open Seat',
                 'VtTotPct', 'Voter Turnout %',
                 'GovWkPct', 'Pct Gov\'t Worker',
                 'UrbanPct', 'Pct Urban',
                 'BlackPct', 'Pct Black',
                 'ForgnPct', 'Pct Foreign Born'),
               ncol = 2, byrow = TRUE)
## Parameters for plotting
r <- 10000
varline <- 6
nline <- 4
tline <- 2
cline <- 0

## Balance plot
setwd(plot.dir)
pdf(paste("BalancePlot", Sys.Date(), ".pdf", sep = ""),
    width = 7, height = 8)
par(mar = c(2, 12, 1, 0))
plot(x = NULL, y = NULL, xlim = c(0, 1), ylim = c(1, nrow(covs) + 3),
     ylab = '', xlab = '', xaxt = "n", yaxt = "n", bty = 'n')
mtext(text = c('Variable\nName', 'Valid\nCases',
        'Treated\nMean', 'Control\nMean'),
      side = 2, font = 2,
      line = c(varline + 3, nline + .7, tline + .7, cline + 0.3),
      adj = .5, las = 2, at = 29, cex = .7)
## For each covariate...
for(i in 1:nrow(covs)) {
  print(covs[i, 2])
  print(sum(is.na(rd[use.5, covs[i, 1]])))
  aty <- nrow(covs) - i + 1
  print(aty)
  mtext(text = covs[i, 2], side = 2, line = varline, adj = 1, las = 2,
        at = aty, cex = .7)
  ## Number of valid cases
  nonna <- sum(!is.na(rd[use.5, covs[i, 1]]))
  ## Mean of treated
  meanT <- signif(mean(rd[(use.5 & rd$DemWin == 1), covs[i, 1]],
                      na.rm = TRUE),
                 digits = 2)
  ## Adding/subtracting digits for presentation purposes
  if(abs(meanT) < 0.1) {
    meanT <- signif(meanT, digits = 1)
  }
  if(meanT %% 1 == 0 & abs(meanT) < 10) {
    meanT <- paste(meanT, '.0', sep = '')
  }
  if(abs(as.numeric(meanT)) >= .1 & abs(as.numeric(meanT)) < 1 &
     nchar(meanT) == 3) {
    meanT <- paste(meanT, '0', sep = '')
  }
  ## Mean of control
  meanC <- signif(mean(rd[(use.5 & rd$DemWin == 0), covs[i, 1]],
                      na.rm = TRUE),
                 digits = 2)
  ## Presentation adjustments
  if(abs(meanC) < 0.1) {
    meanC <- signif(meanC, digits = 1)
  }
  if(meanC %% 1 == 0 & abs(meanC) < 10) {
    meanC <- paste(meanC, '.0', sep = '')
  }
  if(as.numeric(meanC) %% .1 == 0 & abs(as.numeric(meanC)) < 1) {
    meanC <- paste(meanC, '0', sep = '')
  }
  mtext(text = c(nonna, meanT, meanC), side = 2,
        line = c(nline, tline, cline),
        adj = 1, las = 2, at = aty, cex = .7) 
  ## Gray bands
  if(aty %% 2 == 1) {
    polygon(x = c(0, 0, 1, 1),
            y = c(aty - .5, aty + .5, aty + .5, aty - .5),
            border = FALSE,
            col = 'lightgray')
  }
  ## If the variable has three or more levels
  if(length(levels(factor(rd[, covs[i, 1]]))) >= 3) {
    p1 <- pvalue(wilcox_test(rd[, covs[i, 1]][use.5] ~
                             factor(rd$DemWin[use.5]),
                             distribution = 'exact'))
    print(p1)
    sym <- 18
  } else
  ## If the variable is dichotomous
  {
    p1 <- fisher.test(x = factor(rd[, covs[i, 1]][use.5]),
                      y = factor(rd$DemWin[use.5])
                      )$p.value
    sym <- 20
  }
  ## Plot p-value
  points(pch = sym, x = p1, y = aty)
}
## For each outcome variable...
for(i in 1:nrow(dvs)) {
  print(dvs[i, 2])
  aty <- nrow(covs) + nrow(dvs) - i + 1
  print(aty)
  mtext(text = dvs[i, 2], side = 2, line = varline, adj = 1, las = 2,
        at = aty, cex = .7, font = 3)
  nonna <- sum(!is.na(rd[use.5, dvs[i, 1]]))
  meanT <- signif(mean(rd[(use.5 & rd$DemWin == 1), dvs[i, 1]],
                      na.rm = TRUE),
                 digits = 2)
  if(abs(meanT) < 0.1) {
    meanT <- signif(meanT, digits = 1)
  }
  if(meanT %% 1 == 0 & abs(meanT) < 10) {
    meanT <- paste(meanT, '.0', sep = '')
  }
  if(as.numeric(meanT) %% .1 == 0 & abs(as.numeric(meanT)) < 1) {
    meanT <- paste(meanT, '0', sep = '')
  }
  meanC <- signif(mean(rd[(use.5 & rd$DemWin == 0), dvs[i, 1]],
                      na.rm = TRUE),
                 digits = 2)
  if(abs(meanC) < 0.1) {
    meanC <- signif(meanC, digits = 1)
  }
  if(meanC %% 1 == 0 & abs(meanC) < 10) {
    meanC <- paste(meanC, '.0', sep = '')
  }
  if(as.numeric(meanC) %% .1 == 0 & abs(as.numeric(meanC)) < 1) {
    meanC <- paste(meanC, '0', sep = '')
  }
  mtext(text = c(nonna, meanT, meanC), side = 2,
        line = c(nline, tline, cline),
        adj = 1, las = 2, at = aty, cex = .7, font = 3) 
  if(aty %% 2 == 1) {
    polygon(x = c(0, 0, 1, 1),
            y = c(aty - .5, aty + .5, aty + .5, aty - .5),
            border = FALSE,
            col = 'lightgray')
  }
  if(length(levels(factor(rd[, covs[i, 1]]))) >= 3) {
    p1 <- pvalue(wilcox_test(rd[, dvs[i, 1]][use.5] ~
                             factor(rd$DemWin[use.5]),
                             distribution = 'exact'))
    print(p1)
    sym <- 18
  } else {
    p1 <- p2 <- NA
    p1 <- fisher.test(x = factor(rd[, dvs[i, 1]][use.5]),
                      y = factor(rd$DemWin[use.5])
                      )$p.value
    sym <- 20
  }
  points(pch = sym[1], x = p1, y = aty)
  points(pch = sym[2], x = p2, y = aty)
}
segments(x0 = 0, x1 = 0, y0 = .5, y1 = 28.5)
segments(x0 = 0, x1 = 1, y0 = .49, y1 = .49)
segments(x0 = c(.05, .1), x1 = c(.05, .1),
         y0 = .5, y1 = 28.5,
         lty = 'dotted')
segments(x0 = 0, x1 = 1, y0 = 25.5, y1 = 25.5, lty = 'dashed')
mtext(side = 1, at = c(0, .05, .1, 1), text = c('0', '.05', '.1', '1'),
      cex = .7, line = -.75)
mtext(side = 1, at = .5, text = 'p-value')
dev.off()

##################
#### Spending ####
##################
spend.pct.winner <-
  ifelse(rd$DemWin == 1, rd$DSpndPct, 100 - rd$DSpndPct)
## CQ tossups
summary(spend.pct.winner[toss])
mean(spend.pct.winner[toss] > 50, na.rm = T)
t.test(spend.pct.winner[toss] - 50)
wilcox.test(spend.pct.winner[toss] - 50, exact = TRUE)
## 0.5% window
summary(spend.pct.winner[use.5])
mean(spend.pct.winner[use.5] > 50, na.rm = T)
t.test(spend.pct.winner[use.5] - 50)
wilcox.test(spend.pct.winner[use.5] - 50, exact = TRUE)
## CQ tossups in 0.5% window
summary(spend.pct.winner[use.5 & toss])
mean(spend.pct.winner[use.5 & toss] > 50, na.rm = T)
t.test(spend.pct.winner[use.5 & toss] - 50)
wilcox.test(spend.pct.winner[use.5 & toss] - 50, exact = TRUE)

###################
#### Donations ####
###################
donation.pct.winner <-
  ifelse(rd$DemWin == 1, rd$DDonaPct, 100 - rd$DDonaPct)
## CQ tossups
summary(donation.pct.winner[toss])
mean(donation.pct.winner[toss] > 50, na.rm = T)
t.test(donation.pct.winner[toss] - 50)
wilcox.test(donation.pct.winner[toss] - 50, exact = TRUE)
## 0.5% window
summary(donation.pct.winner[use.5])
mean(donation.pct.winner[use.5] > 50, na.rm = T)
t.test(donation.pct.winner[use.5] - 50)
wilcox.test(donation.pct.winner[use.5] - 50, exact = TRUE)
## CQ tossups in 0.5% window
summary(donation.pct.winner[use.5 & toss])
mean(donation.pct.winner[use.5 & toss] > 50, na.rm = T)
t.test(donation.pct.winner[use.5 & toss] - 50)
wilcox.test(donation.pct.winner[use.5 & toss] - 50, exact = TRUE)

####################
#### Experience ####
####################
with(rd[use.5, ], mean(DExpAdv, na.rm = T) + mean(RExpAdv, na.rm = T))
with(rd[use.5 & (rd$DExpAdv == 1 | rd$RExpAdv == 1), ],
     mean((DExpAdv == 1 & DemWin == 1) |
          (RExpAdv == 1 & DemWin == 0),
          na.rm = T))

#############################
#### Best balance window ####
#############################
## Windows
window <- seq(10, 0.25, -0.01)
### DWinPrv
dwp.prop.diff <- rep(NA, length(window)) 
dwn.prop.diff <- rep(NA, length(window)) 
dwp.min.pv <- rep(NA, length(window))
for (i in 1:length(window)) {
  print(i)
  win <- (abs(rd$DifDPct) < window[i] + 0.25 &
          abs(rd$DifDPct) > window[i] - 0.25 &
          use)
  tt <- t.test(DWinPrv ~ DemWin, data = rd[win, ])
  dwp.prop.diff[i] <- tt$estimate[2] - tt$estimate[1]
  dwp.min.pv[i] <- tt$p.value
  tt2 <- t.test(DWinNxt ~ DemWin, data = rd[win, ])
  dwn.prop.diff[i] <- tt2$estimate[2] - tt2$estimate[1]
}
bal.df <- data.frame(dwp.diff = dwp.prop.diff,
                     dwn.diff = dwn.prop.diff,
                     window = window)

mypdf("OptBalWinDWP", width = 9, height = 6)
ggplot(aes(window, dwp.diff), data = bal.df) +
  geom_hline(yintercept = 0, color = "black", size = .5) +
  geom_vline(xintercept = 0, color = "black", size = .5) +
  geom_point(colour = alpha("black", 1/2)) +
  geom_smooth(se = FALSE, color = "black", size = 1.5,
              method = "loess", span = .3) +
  ylab("Treated-Control Difference") +
  xlab("Absolute Value of Democratic Margin (%)") +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  scale_y_continuous(breaks = seq(0, 0.6, 0.1)) +
  opts(title = "Previous Democratic Victory")
dev.off()

## In color        
setwd(plot.dir)
mypdf("OptBalWinDWPvsDWN", width = 9, height = 6)
ggplot(aes(window, dwn.diff), data = bal.df) +
  geom_hline(yintercept = 0, colour = "black", size = .5) +
  geom_vline(xintercept = 0, colour = "black", size = .5) +
  geom_point(alpha = I(1/3), shape = 19, colour = "red") +
  geom_smooth(se = FALSE, colour = "red",
              size = 1,
              method = "loess", span = .15) +
  labs(x = "Absolute Value of Democratic Margin (%), Election t",
       y = "Difference in Proportion") +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  #scale_y_continuous(breaks = seq(0, 0.6, 0.1)) +
  opts(title = "Treated-Control Difference in Proportion of Democratic Victories, t+1 vs. t-1") +
  geom_point(mapping = aes(window, dwp.diff),
             shape = 19,
             colour = I("blue"),
             alpha = I(1/3),
             data = bal.df) +
  geom_smooth(mapping = aes(window, dwp.diff),
              colour = I("blue"),
              linetype = "solid",
              data = bal.df,
              se = FALSE, size = 1,
              method = "loess", span = .15) +
  geom_text(data = data.frame(x = c(9.75, 9.75), y = c(0.79, 0.47)),
            mapping = aes(x, y,  label = c("t+1", "t-1")))
dev.off()

### More covariates
covs <- c("DWinPrv", "DPctPrv", "DifDPPrv", "IncDWNOM1", "ElcSwing",
           "SoSDem", "GovDem", "YearElec", "DifPVDec", "DemInc", "NonDInc")
summary(rd[use, covs])

min.pv <- rep(NA, length(window))
for (i in 1:length(window)) {
  ttp <- vector(mode = "numeric", length = length(covs))
  #ksp <- vector(mode = "numeric", length = length(covs))
  print(i)
  win <- (abs(rd$DifDPct) < window[i] + 0.25 &
          abs(rd$DifDPct) > window[i] - 0.25 &
          use)
  for (j in 1:length(covs)) {
    ttp[j] <- t.test(rd[win, covs[j]] ~ rd$DemWin[win])$p.value
    #ksp[j] <- ks.boot(rd[win, covs[j]], rd$DemWin[win],
    #                 nboots = 1000)$ks.boot.value
  min.pv[i] <- min(ttp)
  }
}
cov.df <- data.frame(pvalue = min.pv, window = window)
                              
## In color
setwd(plot.dir)
mypdf("OptBalWinCovs", width = 9, height = 6)
ggplot(aes(window, pvalue), data = cov.df) +
  geom_hline(yintercept = 0, color = "black", size = .5) +
  geom_vline(xintercept = 0, color = "black", size = .5) +
  geom_point(alpha = (1/2)) +
  geom_smooth(se = FALSE, size = 1,
              method = "loess", span = .15) +
  ylab("Minimum p-Value") +
  xlab("Absolute Value of Democratic Margin (%), Election t") +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  opts(title = "Balance Tests of Ten Covariates in Disjoint 0.5% Intervals")
dev.off()

### Differences in campaign spending in "optimal" window, among all
### elections (n = 53) for which data are available
opt.win <- with(rd, abs(DifDPct) < 1.75 & abs(DifDPct) > 1.25 &
                !is.na(DifDPct))
t.test(DSpndPct ~ DemWin, data = rd[opt.win, ], alternative = "less")


##############################
#### Sensitivity analysis ####
##############################
sumTestSensUnmatch <- function(T, q, n, m, Gamma) { 
  ## Created by Devin Caughey on 26 February 2010
  ## Last modified on 2 December 2010
  ##
  ## The function 'sumTestSensUnmatch' is an implementation of the method of
  ## sensitivity analysis for comparing two unmatched groups described in
  ## Section 4.6 of Paul Rosenbaum "Observational Studies" (2nd Ed., 2002).
  ## It is designed to be used for sum statistics, such as Wilcoxon's rank sum
  ## statistic.
  ##
  ## It takes five arguments:
  ##   T, the observed value of the test statistic (e.g., the sum of the ranks
  ## of the responses of the treated group);
  ##   q, a vector of functions of the responses (e.g., their ranks; note 
  ## that a higher rank corresponds to a higher response), sorted in
  ## decreasing order (don't forget to do this).
  ##   n, the total number of observations
  ##   m, the number of treated observations
  ##   Gamma, an upper limit on the ratio of the a priori odds of treatment
  ## assignment between the treated and control groups.
  ##
  ## This function prints the upper bound of the normal approximation
  ## p-value for the test at the given value of Gamma. It also
  ## invisibly returns a list of important intermediate statistics.
  ##
  ## An example of how this function may be used, taken from Section
  ## 4.6 of Rosenbaum (2002), is provided at the end of this code.
  
  K <- 0:n
  u <- matrix(nrow = n + 1, ncol = n)
  gamma <- log(Gamma)
  rho.i <- matrix(data = NA, nrow = n + 1, ncol = n)
  rho.ij <- array(data = NA, dim = c(n + 1, n, n))
  mu.T <- rep(NA, n + 1)
  var.T <- rep(NA, n + 1)
  sd.T <- rep(NA, n + 1)
  deviate <- rep(NA, n + 1)                   
  
  for(k in K) {
    ## Fix u[k+1,] so that it consists of k 0's followed by n-k 1's.
    u[k + 1,] <- c(rep(1, k), rep(0, n - k))
    ## Define function 'zeta'.
    zeta <- function(n, m, k, Gamma) {
      zeta.out <- 0
      if(m >= 0 & k >= 0) { 
        for(a in max(0, m + k - n):min(m, k)) { 
          zeta.out <- (zeta.out +
                       (choose(k, a) *
                        choose((n - k), (m - a)) *
                        (Gamma ^ a)))
        }
      }
      return(zeta.out)
    }
    ## Calculate unit i's probability of treatment
    for(i in 1:n) {
      rho.i[k + 1, i] <-
        exp(gamma * u[k + 1, i]) *
          zeta(n - 1, m - 1, k - u[k + 1, i], Gamma) /
            zeta(n, m, k, Gamma)
    }
    ## Calculate i and j's joint probability of treatment.
    for(i in 1:n) {
      for(j in 1:n) {
        if(i == j) {
          rho.ij[k + 1, i, j] <- rho.i[k + 1, i]
        } else {
          rho.ij[k + 1, i, j] <-
            ((exp(gamma * (u[k + 1, i] + u[k + 1, j])) *
              zeta(n - 2,
                   m - 2,
                   k - u[k + 1, i] - u[k + 1, j],
                   Gamma)) / zeta(n, m, k, Gamma))
        }
      } 
    }
    ## Calculate mean of T under the null
    mu.T[k + 1] <- q %*% rho.i[k + 1, ]
    ## Calculate standard deviation of T under the null
    var.T[k + 1] <- 0
    for(i in 1:n) {
      for(j in 1:n) {
        var.T[k + 1] <-
          (var.T[k + 1] +
           (q[i] * q[j] *
           (rho.ij[k + 1, i, j] - rho.i[k + 1, i] * rho.i[k + 1, j])))
      }
    }
    sd.T[k + 1] <- sqrt(var.T[k + 1])
    ## Calculate deviate
    deviate[k + 1] <- (T - mu.T[k + 1]) / sd.T[k + 1]
  }
  ## Main result
  minDeviate <- min(deviate)
  pValueUB <- pnorm(q = minDeviate,
                    mean = 0,
                    sd = 1,
                    lower.tail = FALSE)
  pValueUB.print <- ifelse(pValueUB < 0.0001,
                           "< 0.0001",
                           as.character(round(pValueUB, 4)))
  ## Collect output
  kMin <- K[which(deviate == min(deviate))]
  output <- list(pValueUB,
                 minDeviate,
                 deviate,
                 kMin,
                 T,
                 mu.T,
                 var.T,
                 sd.T,
                 rho.i,
                 rho.ij)
  names(output) <- c("pValueUB",
                     "minDeviate",
                     "deviate",
                     "kMin",
                     "T",
                     "mu.T",
                     "var.T",
                     "sd.T",
                     "rho.i",
                     "rho.ij")
  cat("For Gamma = ",
      Gamma, 
      ", the upper-bound on the p-value of the sum test is: ", 
      pValueUB.print,
      ".\n",
      sep = "")                 
  invisible(output)
}

valid.obs <- close.5 & use & !is.na(rd$DPctNxt) & !is.na(rd$DPctPrv)

wilcox_test(DifDPNxt ~ factor(DemWin), data = rd[valid.obs, ],
            alternative = "less", distribution = "exact")
wilcox_test(DifDPPrv ~ factor(DemWin), data = rd[valid.obs, ],
            alternative = "less", distribution = "exact")
wilcox_test(DPctNxt ~ factor(DemWin), data = rd[valid.obs, ],
            alternative = "less", distribution = "exact")
wilcox_test(DPctPrv ~ factor(DemWin), data = rd[valid.obs, ],
            alternative = "less", distribution = "exact")

rd.q <- sort(rank(rd$DPctNxt[valid.obs]), decreasing = TRUE)
rdT <- rank(rd$DPctNxt[valid.obs]) %*% rd$DemWin[valid.obs]
rdN <- sum(valid.obs)
rdM <- sum(rd$DemWin[valid.obs])

# system.time(rdOut <- sumTestSensUnmatch(T = rdT,
#                                         q = rd.q,
#                                         n = rdN,
#                                         m = rdM,
#                                         Gamma = 1))
## Takes over four minutes per value of Gamma.

system.time(for (g in seq(3.6, 3.7, .1)) {
  print(g)
  sumTestSensUnmatch(T = rdT, q = rd.q, n = rdN, m = rdM, Gamma = g)
})
## > For Gamma = 3.7, the upper-bound on the p-value of the sum test is: 0.0541.
## > For Gamma = 3.6, the upper-bound on the p-value of the sum test is: 0.0489.

### Sensitivity for DWinPrv and DWinNxt
with(rd[valid.obs, ], table(DemWin, DWinNxt))
with(rd[valid.obs, ], table(DemWin, DWinPrv))
with(rd[valid.obs, ],
     fisher.test(DemWin, DWinNxt, alternative = "greater"))
with(rd[valid.obs, ],
     fisher.test(DemWin, DWinPrv, alternative = "greater"))

FisherSens <- function(totalN,
                       numberTreated,
                       numberSuccesses,
                       treatedSuccesses,
                       Gammas) { 
  ## 'FisherSens' 
  ## Created by Devin Caughey on 1 March 2010
  ## Last modified: 4 December 2010
  ##
  ## This function performs a sensitivity analysis for Fisher's Exact Test
  ## for count (binary) data. It is derived from Section 4.4 of Paul
  ## Rosenbaum’s ``Observational Studies" (2nd Ed., 2002). The test is
  ## one-sided; that is, the alternative hypothesis to the null is
  ## that the number of "successes" (however defined) greater in among
  ## the "treated" (however defined).
  ##
  ## It takes five arguments:
  ##   'totalN': total number of observations
  ##   'numberTreated': number of observations that received treatment
  ##   'numberSuccesses': total number of "successes"--i.e., 1's
  ##   'treatedSuccesses': number of successes among the treated
  ##   'Gammas': a vector Gammas (bounds on the differential odds of 
  ##  treatment) at which to test the significance of the results.
  ##
  ## It returns a matrix of Gammas and upper and lower bounds on the exact 
  ## p-value for Fisher's test.
  
  n <- totalN
  m <- numberTreated
  c_plus <- numberSuccesses
  a <- treatedSuccesses
  
  Upsilon <- function(n, m, c_plus, a, Gamma) {
    ## Prob A >= a 
    numer <- 0
    for(k in max(a, (m + c_plus - n)):min(m, c_plus)) {
      numer <- numer +
        choose(c_plus, k) * choose((n-c_plus), m - k) * (Gamma ^ k)
    }
    denom <- 0
    for(k in max(0, (m + c_plus - n)):min(m, c_plus)) {
      denom <- denom +
        choose(c_plus, k) * choose((n - c_plus), (m - k)) * Gamma ^ k
    }
    return(numer / denom)
  }
  p_plus <- rep(NA, length(Gammas))
  p_minus <- rep(NA, length(Gammas))
  
  for(g in 1:length(Gammas)) {
    p_plus[g] <- Upsilon(n = n, m = m, c_plus = c_plus,
                         a = a, Gamma = Gammas[g])
    p_minus[g] <- Upsilon(n = n, m = m, c_plus = c_plus,
                          a = a, Gamma = (1 / Gammas[g]))
  }
  
  output <- cbind(Gammas, p_minus, p_plus)
  dimnames(output)[[2]] <- c("Gamma", "P-Value LB", "P-Value UB")
  return(output)
}

trSuc <- sum(rd$DWinNxt[valid.obs] == 1 & rd$DemWin[valid.obs] == 1)

FisherSens(totalN = rdN,
           numberTreated = rdM,
           numberSuccesses = sum(rd$DWinNxt[valid.obs]),
           treatedSuccesses = trSuc,
           Gammas = seq(1, 5, .1))
## > Crosses one-side p-value of 0.05 at Gamma = 2.4

####################################################
#### Bootstrapped local linear regression plots ####
####################################################

### Define functions
EdgeKernel <- function(X, mid, bw, norm = FALSE) {
  ## Created by: Devin Caughey
  ## Last modified: 12 April 2011
  
  ## Invisibly returns a vector of edge (triangular) kernel weights.
  ## If an observation's distance from the midpoint is greater than
  ## bw, it receives zero weight. If it is within the bandwidth, it
  ## receives a weight equal to 1 minus the distance over the
  ## bandwidth.
  ## 'norm' regulates whether the weights should be normalized to a
  ## mean of 1 (and thus to a sum equal to the number of non-missing
  ## observations)
  
  wts.out <- ifelse(abs(X - mid) > bw, 0, 1 - (abs(X - mid) / bw))
  if(norm) {
    mean.wts <- mean(wts.out, na.rm = TRUE)
    wts.out <- wts.out /  mean.wts
  }
  
  invisible(wts.out)
}

EdgeLocalEst <- function(Y, X, bw, point,
                         norm = FALSE, print = FALSE,
                         verbose = FALSE) {
  ## Returns the estimate of the local regression function at
  ## 'point', with weights determined by an edge kernel with a
  ## bandwidth of bw.
  Xcentered <- X - point
  wts.c <- EdgeKernel(X = Xcentered, mid = 0, bw = bw, norm = TRUE)
  llr.c <- lm(Y ~ Xcentered, weights = wts.c)
  if(verbose) {
    print(summary(llr.c))
  }
  PointEst <- coef(llr.c)[[1]]
  SEofEst <- sqrt(diag(vcov(llr.c)))[1]
  names(SEofEst) <- NULL
  out <- list(PointEst, SEofEst)
  names(out) <- c("Point Estimate", "Standard Error of Estimate")
  if(print) {
    return(out)
  }
  if(!print) {
    invisible(out)
  }
}

EdgeLocalDiffEst <- function(Y, X, bw, point,
                             norm = FALSE, print = FALSE,
                             verbose = FALSE) {
  ## Returns the difference in the value of a function at a given
  ## point estimated from above and from below, using local linear
  ## regression with an edge kernel.
  Xcentered <- X - point
  wts.c <- EdgeKernel(X = Xcentered, mid = 0, bw = bw, norm = TRUE)
  trmt <- as.numeric(Xcentered >= 0)
  llr.c <- lm(Y ~ trmt + Xcentered + I(trmt * Xcentered),
              weights = wts.c)
  if(verbose) {
    print(summary(llr.c))
  }
  DiffEst <- coef(llr.c)[[2]]
  SEofEst <- sqrt(diag(vcov(llr.c)))[2]
  names(SEofEst) <- NULL 
  out <- list(DiffEst, SEofEst)
  names(out) <- c("Estimated Difference", "Standard Error of Estimate")
  if(print) {
    return(out)
  }
  if(!print) {
    invisible(out)
  }
}

EdgeSmoothBoot <- function(dat, inds, y.name, x.name, points, bw) {
  y <- dat[inds, y.name]
  x <- dat[inds, x.name]
                                        # print(summary(y))
                                        # print(summary(x))
  fit <- vector(mode = 'numeric', length = length(points))
  for(i in 1:length(points)) {
    p <- points[i]
                                        # print(p)
    fit[i] <-
      unlist(EdgeLocalEst(Y = y, X = x, point = p, bw = bw)[1])
  }
                                        # fit.df <- as.data.frame(fit)
                                        # rownames(fit.df) <- points
  invisible(fit)
}

## Test weights (normalized to mean = 1)
wts.test <- EdgeKernel(X = rd$DifDPct, mid = 0, bw = 10, norm = TRUE)
summary(wts.test)
plot(x = rd$DifDPct, y = wts.test)
## Test difference estimator
EdgeLocalDiffEst(Y = rd$DPctPrv, X = rd$DifDPct,
                 point = 0, bw = 6, print = TRUE, verbose = TRUE)

## Definitions
breaks.5 <- seq(-100, 100, .5)
mids.5 <- hist(rd$DifDPct, freq = TRUE, breaks = breaks.5)$mids
cuts.5 <- cut(rd$DifDPct, breaks = breaks.5)
covs.llr <- matrix(c('DWinPrv', 'Dem Win t - 1',
                 'DPctPrv', 'Dem % t - 1',
                 'DifDPPrv', 'Dem % Margin t - 1',
                 'IncDWNOM1', 'Inc\'s 1st Dim NOMINATE',
                 'DemInc', 'Dem Inc in Race',
                 'NonDInc', 'Rep Inc in Race',
                 'RExpAdv', 'Rep Experience Adv',
                 'DExpAdv', 'Dem Experience Adv',
                 'ElcSwing', 'Partisan Swing',
                 'DSpndPct', 'Dem Spending %',
                 'DDonaPct', 'Dem Donation %',
                 'SoSDem', 'Dem Sec of State',
                 'GovDem', 'Dem Governor',
                 'DifPVDec', 'Dem Pres % Margin', ## Average over decade
                 'DemOpen', 'Dem-held Open Seat',
                 'NonDOpen', 'Rep-held Open Seat',
                 'OpenSeat', 'Open Seat',
                 'VtTotPct', 'Voter Turnout %',
                 'GovWkPct', 'Pct Gov\'t Worker',
                 'UrbanPct', 'Pct Urban',
                 'BlackPct', 'Pct Black',
                 'ForgnPct', 'Pct Foreign Born'
                 ),
               ncol = 2, byrow = TRUE
               )
bws <- c(1, 2, 5, 10, 15, 20)
left <- rd$DifDPct < 0 & !is.na(rd$DifDPct)
right <- rd$DifDPct > 0 & !is.na(rd$DifDPct)
left.points <- seq(-25, 0, 0.5)
right.points <- seq(0, 25, 0.5)
bw.names <- paste('bw', bws, sep = '')
nboots <- 1000
boot.out <- vector(mode = 'list', length = nrow(covs.llr))
names(boot.out) <- covs.llr[, 1]
rd.sub <- rd[, c('DifDPct', covs.llr[, 1])]
dim(rd.sub)


## Note: Each plot may take a half an hour or more to calculate.
## Loop over covariates---uncomment below to do all covariates.
#for (i in 1:nrow(covs.llr)) {
for (i in 1) {
  cat(covs.llr[i, 1], ' (', as.character(Sys.time()), ')', '\n', sep = '')
  boot.out[[i]] <- vector(mode = 'list', length = length(bws))
  names(boot.out[[i]]) <- bw.names
  mean.5 <- tapply(rd.sub[, covs.llr[i, 1]], cuts.5, function(x) {
    mean(x, na.rm = TRUE)
  })
  plot.data <- data.frame(Midpoints = mids.5[abs(mids.5) < 25],
                          BinMeans = mean.5[abs(mids.5) < 25])     
  pdf.name <- paste("EdgeBoot_", covs.llr[i, 1], "_", sep = "")
  mypdf(pdf.name)
  ## Loop over bandwidths---uncomment below to do all bandwidths
  #for (h in 1:length(bws)) {
  for (h in 1) {
    cat(bws[h], '\n')
    boot.out[[i]][[h]] <- vector(mode = 'list', length = 2)
    names(boot.out[[i]][[h]]) <- c('Left', 'Right')
    boot.out[[i]][[h]]$Left <- vector(mode = 'list')
    boot.out[[i]][[h]]$Right <- vector(mode = 'list')
    left.txt <-
      paste("boot(data = rd.sub[left, ], statistic = EdgeSmoothBoot, R = nboots, sim = 'ordinary', y.name = '",
            covs.llr[i, 1],
            "', x.name = 'DifDPct', points = left.points, bw = ",
            bws[h], ")",
            sep = '')
    right.txt <-
      paste("boot(data = rd.sub[right, ], statistic = EdgeSmoothBoot, R = nboots, sim = 'ordinary', y.name = '",
            covs.llr[i, 1],
            "', x.name = 'DifDPct', points = right.points, bw = ",
            bws[h], ")",
            sep = '')
    boot.out[[i]][[h]]$Left <- eval(parse(text = left.txt))
    boot.out[[i]][[h]]$Right <- eval(parse(text = right.txt))
    file.name <- paste('EdgeKernel', covs.llr[i, 1], bws[h], 'boot',
                       Sys.Date(), '.RData',
                       sep = '')
    ## Assign boot output to stand-alone object
    eval(parse(text = paste(covs.llr[i, 1],
                 '.boot <- boot.out[[i]][[h]]',
                 sep = '')))
    ## Save boot object to disk
    eval(parse(text = paste('save(',
                 covs.llr[i, 1],
                 '.boot, file = file.name)',
                 sep = '')))
    left.ci <- matrix(NA, nrow = 2,
                      ncol = length(left.points))
    for(p in 1:length(left.points)) {
      left.ci[, p] <- quantile(boot.out[[i]][[h]]$Left$t[, p],
                               c(.025, .975))
    }
    right.ci <- matrix(NA, nrow = 2,
                       ncol = length(right.points))
    for(p in 1:length(right.points)) {
      right.ci[, p] <- quantile(boot.out[[i]][[h]]$Right$t[, p],
                                c(.025, .975))
    }
    ldata <-
      data.frame(Points = left.points,
                 FittedEst = boot.out[[i]][[h]]$Left$t0,
                 LowerCI = left.ci[1, ],
                 UpperCI = left.ci[2, ])
    rdata <-
      data.frame(Points = right.points,
                 FittedEst = boot.out[[i]][[h]]$Right$t0,
                 LowerCI = right.ci[1, ],
                 UpperCI = right.ci[2, ])
    binplot <- qplot(x = Midpoints, y = BinMeans, data = plot.data,
                     ylab = 'Local Average',
                     xlab = 'Democratic Vote Percentage Margin',
                     main = paste(covs.llr[i, 2],' (h = ', bws[h], 
                       ')', sep = ''))
    ci.aes <- aes(x = c(Points, rev(Points)),
                  y = c(LowerCI, rev(UpperCI)))
    fit.aes <- aes(x = Points, y = FittedEst)
    print(binplot +
          geom_polygon(mapping = ci.aes, data = ldata,
                       fill = 'grey') +
          geom_line(mapping = fit.aes, data = ldata, size = .5) +
          geom_polygon(mapping = ci.aes, data = rdata,
                       fill = 'grey') +
          geom_line(mapping = fit.aes, data = rdata, size = .5) +
          geom_point(mapping = aes(x = Midpoints, y = BinMeans),
                     data = plot.data) + 
          geom_vline(xintercept = 0))    
    ## Drop object so doesn't take up memory
    eval(parse(text = paste('rm(', covs.llr[i, 1], '.boot)', sep = '')))
    boot.out[[i]][[h]]$Left <- boot.out[[i]][[h]]$Left$t
    boot.out[[i]][[h]]$Right <- boot.out[[i]][[h]]$Right$t
  } ## End bandwidth loop
  dev.off()
} ## End covariate loop

#############################
#### Imbalance over time ####
#############################
inc.win <- with(rd, as.numeric(DemWin == DWinPrv))
winner.pct.prv <- with(rd, ifelse(DemWin == 1, DifDPPrv, -DifDPPrv))
imbal.df.5 <- data.frame(year = rd$YearElec[close.5],
                         inc.win = inc.win[close.5],
                         winner.pct.prv = winner.pct.prv[close.5])

pdf(file = paste("ImbalIncWin.5Time", Sys.Date(), ".pdf", sep = ""),
    height = 4, width = 7)
jit <- position_jitter(width = 0, height = .01)
ggplot(data = imbal.df.5, aes(x = year, y = inc.win)) +
  geom_jitter(position = jit, colour = alpha("black", 1/2)) +
  geom_smooth() +
  ylab("Incumbent Party Win") +
  xlab("Year") +
  opts(title = "Imbalance over Time in 0.5% Window") +
  geom_hline(yintercept = .5, color = "black")
dev.off()
## >> Two peaks (1940 & mid-1980s) and two troughs (mid-1960s & 2006),
## >> but almost always statistically different from 0.   

setwd(plot.dir)
pdf(file = paste("ImbalIncWin05TimeNoSmooth", Sys.Date(), ".pdf", sep = ""),
    height = 4, width = 7) 
set.seed(4)
jit <- position_jitter(width = 0, height = .2)
ggplot(data = imbal.df.5, aes(x = year, y = inc.win)) +
  geom_jitter(position = jit, colour = alpha("black", 1/2)) +  
  ylab("Incumbent Party Win (Jittered)") +
  xlab("Year") +
  opts(title = "Incumbent Victories and Losses in 0.5% Window, Over Time") +
  geom_hline(yintercept = .5, color = "black") +
  scale_y_continuous(breaks = seq(0, 1, 1), )
dev.off()

pdf(file = paste("ImbalPrevMarg.5TimeNoSmooth", Sys.Date(), ".pdf", sep = ""),
    height = 4, width = 7)
ggplot(data = imbal.df.5, aes(x = year, y = winner.pct.prv)) +
  geom_jitter(position = jit, colour = alpha("black", 1/2)) +
  ylab("Winning Party's Previous Margin") +
  xlab("Year") +
  opts(title = "Imbalance over Time in 0.5% Window") +
  geom_hline(yintercept = 0, color = "black")
dev.off()
## >> Peak in mid-1970s and trough in late 1980s - 2000, during which
## >> not statistically different from 0.

imbal.df1 <- data.frame(year = rd$YearElec[close1],
                         inc.win = inc.win[close1],
                         winner.pct.prv = winner.pct.prv[close1])

pdf(file = paste("ImbalIncWin1Time", Sys.Date(), ".pdf", sep = ""),
    height = 4, width = 7)
jit <- position_jitter(width = 0.5, height = 0)
ggplot(data = imbal.df1, aes(x = year, y = inc.win)) +
  geom_jitter(position = jit, colour = alpha("black", 1/3)) +
  geom_smooth() +
  ylab("Incumbent Party Win") +
  xlab("Year") +
  opts(title = "Imbalance over Time in 1% Window") +
  geom_hline(yintercept = .5, color = "black")
dev.off()
## >> Two peaks (1940 & mid-1980s) and two troughs (mid-1960s & 2006),
## >> like above, except that the drop from the mid-1980s is much more
## >> dramatic (fitted line crosses 0).   

pdf(file = paste("ImbalPrevMarg1Time", Sys.Date(), ".pdf", sep = ""),
    height = 4, width = 7)
jit <- position_jitter(width = 0, height = 0)
ggplot(data = imbal.df1, aes(x = year, y = winner.pct.prv)) +
  geom_jitter(position = jit, colour = alpha("black", 1/2)) +
  geom_smooth() +
  ylab("Winning Party's Previous Margin") +
  xlab("Year") +
  opts(title = "Imbalance over Time in 1% Window") +
  geom_hline(yintercept = 0, color = "black")
dev.off()
## >> Peak in early 1970s, then gradual decline (not significant from
## >> 0 by 1990)

####################
#### Open seats ####
#################### 
### Previous Democratic Victory
open <- rd$OpenSeat == 1
table(rd$DWinPrv[open & close.25 & use],
      rd$DemWin[open & close.25 & use])
table(rd$DWinPrv[open & close.5 & use],
      rd$DemWin[open & close.5 & use])
table(rd$DWinPrv[open & close1 & use],
      rd$DemWin[open & close1 & use])
t.test(DWinPrv ~ DemWin, data = rd[close.5 & use, ])
t.test(DWinPrv ~ DemWin, data = rd[open & close.25 & use, ])
t.test(DWinPrv ~ DemWin, data = rd[open & close.5 & use, ])
t.test(DWinPrv ~ DemWin, data = rd[open & close1 & use, ])
t.test(DWinPrv ~ DemWin, data = rd[open & close.5, ])
t.test(DWinPrv ~ DemWin, data = rd[open & close1, ])
binom.test(sum(rd$DWinPrv[open & close.5 & use] ==
               rd$DemWin[open & close.5 & use]),
           sum(open & close.5 & use))
binom.test(sum(rd$DWinPrv[close.5 & use] ==
               rd$DemWin[close.5 & use]),
           sum(close.5 & use))

### Spending 
## All open seats 
t.test(DSpndPct ~ DemWin, data = rd[open, ]) 
wilcox_test(DSpndPct ~ factor(DemWin), data = rd[open, ])
## Open seats in 1% window 
t.test(DSpndPct ~ DemWin, data = rd[open & close1, ])   
wilcox_test(DSpndPct ~ factor(DemWin),
            data = rd[open & close1, ],
            distribution = exact())
## Open seats in 1% window and in `use' sample
t.test(DSpndPct ~ DemWin, data = rd[open & close1 & use, ])   
wilcox_test(DSpndPct ~ factor(DemWin),
            data = rd[open & close1 & use, ],
            distribution = exact())
## Open seats in 0.5% window    
t.test(DSpndPct ~ DemWin, data = rd[open & close.5, ])   
wilcox_test(DSpndPct ~ factor(DemWin),
            data = rd[open & close.5, ],
            distribution = exact())     
## Open seats in 0.5% window and in `use' sample
t.test(DSpndPct ~ DemWin, data = rd[open & close.5 & use, ])   
wilcox_test(DSpndPct ~ factor(DemWin),
            data = rd[open & close.5 & use, ],
            distribution = exact())                 

#################################################### 
#### EFFECT ESTIMATES CONDITIONAL ON COVARIATES ####
####################################################
### 0.5% Window, all races, DPctNxt 
summary(lm(DPctNxt ~ DemWin, data = rd[close.5 & use, ]))  
summary(lm(DPctNxt ~ DemWin + DifDPct + I(DemWin * DifDPct), 
  data = rd[close.5 & use, ]))
summary(lm(DPctNxt ~ DemWin + DWinPrv + DifDPct + 
	I(DemWin * DifDPct) + DPctPrv, 
	data = rd[close.5 & use, ])) 
summary(lm(DPctNxt ~ DemWin + DWinPrv + DifDPct + 
	I(DemWin * DifDPct) + DPctPrv + factor(CQRating3),
	data = rd[close.5 & use, ]))
summary(lm(DPctNxt ~ DemWin + DWinPrv + DifDPct + 
	I(DemWin * DifDPct) + DPctPrv + factor(CQRating3) + DSpndPct, 
	data = rd[close.5 & use, ]))        
### 0.5% Window, all races, DWinNxt
summary(glm(DWinNxt ~ DemWin, 
	data = rd[close.5 & use, ],
	family = binomial(link = "logit")))  
summary(glm(DWinNxt ~ DemWin + DifDPct + I(DemWin * DifDPct), 
	data = rd[close.5 & use, ],
	family = binomial(link = "logit")))
summary(glm(DWinNxt ~ DemWin + DifDPct + I(DemWin * DifDPct) + DWinPrv + 
  DPctPrv, 
	data = rd[close.5 & use, ],
	family = binomial(link = "logit"))) 
## >> Marginal significance
summary(glm(DWinNxt ~ DemWin + DifDPct + I(DemWin * DifDPct) + DWinPrv + 
  DPctPrv + factor(CQRating3), 
	data = rd[close.5 & use, ],
	family = binomial(link = "logit")))           
## >> Marginal significance
summary(glm(DWinNxt ~ DemWin + DifDPct + I(DemWin * DifDPct) + DWinPrv + 
  DPctPrv + factor(CQRating3) + DSpndPct, 
	data = rd[close.5 & use, ],
	family = binomial(link = "logit")))      
### 0.5% Window, open seats, DPctNxt 
summary(lm(DPctNxt ~ DemWin, data = rd[close.5 & use & open, ]))  
summary(lm(DPctNxt ~ DemWin + DifDPct + I(DemWin * DifDPct), 
	data = rd[close.5 & use & open, ]))
## >> Marginal significance
summary(lm(DPctNxt ~ DemWin + DWinPrv + DifDPct + I(DemWin * DifDPct) + 
  DPctPrv, 
	data = rd[close.5 & use & open, ])) 
## >> Marginal significance
summary(lm(DPctNxt ~ DemWin + DWinPrv + DifDPct + I(DemWin * DifDPct) + 
  DPctPrv + factor(CQRating3),
	data = rd[close.5 & use & open, ])) 
## >> Marginal significance
summary(lm(DPctNxt ~ DemWin + DWinPrv + DifDPct + I(DemWin * DifDPct) + 
  DPctPrv + factor(CQRating3) + DSpndPct, 
	data = rd[close.5 & use & open, ]))        
### 0.5% Window, open seats, DWinNxt
summary(glm(DWinNxt ~ DemWin, 
	data = rd[close.5 & use & open, ],
	family = binomial(link = "logit")))  
summary(glm(DWinNxt ~ DemWin + DifDPct + I(DemWin * DifDPct), 
	data = rd[close.5 & use & open, ],
	family = binomial(link = "logit")))
## Linear probability model  
summary(lm(DWinNxt ~ DemWin + DifDPct + I(DemWin * DifDPct), 
	data = rd[close.5 & use & open, ]))
summary(lm(DWinNxt ~ DemWin + DifDPct, 
	data = rd[close.5 & use & open, ]))
## >> Negative estimate
summary(glm(DWinNxt ~ DemWin + DifDPct + I(DemWin * DifDPct) + DWinPrv + DPctPrv, 
	data = rd[close.5 & use & open, ],
	family = binomial(link = "logit"))) 
## >> Not significant
summary(glm(DWinNxt ~ DemWin + DifDPct + I(DemWin * DifDPct) + DWinPrv + DPctPrv + 
	factor(CQRating3), 
	data = rd[close.5 & use & open, ],
	family = binomial(link = "logit"))) 
## Can't estimate model below
summary(glm(DWinNxt ~ DemWin + DifDPct + I(DemWin * DifDPct) + DWinPrv + DPctPrv + 
	factor(CQRating3) + DSpndPct, 
	data = rd[close.5 & use & open, ],
	family = binomial(link = "logit"))) 
     
