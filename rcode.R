load(file = "assign3data.RData")
library(MASS)
# The data starts in January 1988, with 12 observations per year (every 28 days)
z <- ts(as.vector(t(snatch)), start = c(1988, 1), frequency = 12)


# Plot the raw time series
plot(z, type = "b", pch = 16,
     main = "Reported Purse Snatchings in Hyde Park (1988–1993)",
     xlab = "Time", ylab = "Number of Incidents")

par(mfrow = c(1, 3))  # Set up a 2x2 plot layout
plot(z, main = "Raw Data", ylab = "Snatchings")
lines(lowess(as.vector(time(z)),z),type='l')
plot(log(z), main = "Log Transformation", ylab = "log(Snatchings)")
lines(lowess(as.vector(time(z)),log(z)),type='l')
plot(sqrt(z), main = "Sqaure Root Transformation", ylab = "sqrt(Snatchings)")
lines(lowess(as.vector(time(z)),sqrt(z)),type='l')

dev.off()

# Plot the raw snatchings time series
plot(z, type = "b", pch = 16,
     main = "Reported Purse Snatchings in Hyde Park (1988–1993)",
     xlab = "Time", ylab = "Number of Incidents")
lines(lowess(time(z), z), col = "red", lwd = 2)
legend("topright", legend = "lowess trend", col = "red", lty = 1, lwd = 2,bty = "n")

#筛选用什么polynomial
# Linear trend: use time(z) directly
fit_lin <- rlm(z ~ time(z))

# Quadratic trend: still no need to center
fit_quad <- rlm(z ~ cbind(time(z), time(z)^2))

# Cubic trend: center time to reduce collinearity
x <- time(z) - 1990
fit_cubic <- rlm(z ~ cbind(x, x^2, x^3))

# Plot comparison
plot(z, type = "b", pch = 16, col = "black",
     main = "Trend Comparison for Snatching Data (1988–1993)",
     xlab = "Time", ylab = "Number of Snatchings")

lines(as.vector(time(z)), fitted(fit_lin), col = "blue", lty = 1, lwd = 2)
lines(as.vector(time(z)), fitted(fit_quad), col = "red", lty = 2, lwd = 2)
lines(as.vector(time(z)), fitted(fit_cubic), col = "darkgreen", lty = 3, lwd = 2)

legend("topright",
       legend = c("Linear Trend", "Quadratic Trend", "Cubic Trend"),
       col = c("blue", "red", "darkgreen"),
       lty = 1:3, lwd = 2, bty = "n")
#choose quadratic
trendfit<-fit_quad

#seasonality
plot(trendfit$resid,
     type = "b", pch = 16,
     main = "Detrended Series (Residuals from Quadratic Trend)",
     xlab = "Time", ylab = "Residuals")
abline(h = 0, col = "blue")

detrended <- ts(trendfit$resid, frequency = 12, start = start(z))

# Step 1: STL decomposition
h <- stl(z, s.window = "periodic")

# Step 2: Plot decomposition
plot(h,lwd=2)
title("STL Decomposition of Purse Snatching Series")

# Step 3: Monthplot of the seasonal component only
monthplot(h,
          ylab = "Seasonal Component",
          main = "Monthplot of Seasonal Component from STL")

# Create a month vector from the time series
v2 <- outer(cycle(z), 1:12, "==") * 1
v1 <- v2
apply(v1, 2, sum)
months <- c("January","February","March","April","May","June","July","August",
            "September","October","November","December")

par(mfrow = c(3, 4), oma = c(6, 0, 6, 0))

for (i in 1:12){
  yvals <- z[v1[,i]==1]
  n <- length(yvals)
  xvals <- 1:n
  
  if(i <= 6){
    yvals <- c(NA, yvals)
    xvals <- 1:(n+1)
  } else {
    yvals <- c(yvals, NA)
    xvals <- 1:(n+1)
  }
  
  plot(xvals, yvals,
       type = "b", pch = 16, axes = FALSE,
       ylim = range(z), xlim = c(1, max(xvals)),
       main = months[i],
       ylab = "Snatchings", xlab = "")
  
  abline(rlm(yvals ~ xvals))
  
  box()
  axis(2)
  axis(1, at = xvals, labels = 1988:(1988 + length(xvals) - 1))
}
mtext("Seasonal Subseries Plots of Snatching Data",
      side = 3, line = 2, outer = TRUE, cex = 1.5)

detrended <- ts(trendfit$residuals, frequency = 12, start = c(1988, 1))
for (i in 1:12){
  yvals <- detrended[v1[,i] == 1]
  n <- length(yvals)
  xvals <- 1:n
  
  # Padding with NA to align (same logic as before)
  if(i <= 6){
    yvals <- c(NA, yvals)
    xvals <- 1:(n+1)
  } else {
    yvals <- c(yvals, NA)
    xvals <- 1:(n+1)
  }
  
  plot(xvals, yvals,
       type = "b", pch = 16, axes = FALSE,
       ylim = range(detrended), xlim = c(1, max(xvals)),
       main = months[i],
       ylab = "Detrended Snatchings", xlab = "")
  
  abline(rlm(yvals ~ xvals))  # robust linear trend on residuals
  box()
  axis(2)
  axis(1, at = xvals, labels = 1988 + (0:(length(xvals)-1)))
}

# Add overall title
mtext("Seasonal Subseries of Detrended Snatchings", side = 3, line = 2, outer = TRUE, cex = 1.5)
mtext("Linear trend removed ", side = 1, line = 2, outer = TRUE)

# Define Factor() to create dummy variable matrix (excluding the baseline level)
Factor <- function(a) {
  code <- factor(a)
  nlev <- length(levels(code))
  x <- matrix(0, length(a), nlev - 1)
  dimnames(x) <- list(NULL, levels(code)[2:nlev])
  for (i in 2:nlev) {
    x[code == levels(code)[i], i - 1] <- 1
  }
  return(x)
}
#build new season vector
v <- Factor(c(
  rep(c(2,2,1,1,1,1,1,1,1,1,2,2), 5),
  2,2,1,1,1,1,1,1,1
))

# Fit model with quadratic trend and one seasonal group (v2) plus interaction
fits <- rlm(z ~ cbind(time(z), v, v[,1] * time(z)))

# Plot fitted trend + season vs raw data
dev.off()
plot(z, xlab = "Month", ylab = "Snatching data", type = "b", pch = 16,
     main = "Monthly Snatching in Hyde Park",
     sub = "Quadratic trend with seasonal group (Nov–Feb) + interaction")
lines(as.vector(time(z)), fitted(fits), lty = 2, col = "red")
legend("topright",
       legend = "Fitted Values (trend + season)",
       col = "red",
       lty = 2,
       lwd = 2,
       bty = "n")

# Seasonal decomposition components
season <- c(
  rep(c(2,2,1,1,1,1,1,1,1,1,2,2), 5),
  2,2,1,1,1,1,1,1,1
)

# Coefficients
b <- fits$coef

# Only one group effect (Group 2)
ss <- rep(0, 2)
ss[2] <- b[3]

# Build seasonal component with interaction term for Group 2
seasonal <- ss[season] +
  ifelse(season == 2, b[4], 0) * time(z)


# Plot decomposition: original, trend, seasonal, irregular
par(mfrow = c(4, 1), oma = c(0, 0, 3, 0))

# Original series
plot(1:length(z), z, type = "b", pch = 16, axes = FALSE,
     ylab = "Snatching", xlab = "")
box()
axis(2)
axis(1, at = seq(1, length(z), 12), labels = 1988:1993)

# Trend
trend <- fitted(trendfit)
plot(1:length(z), trend, type = "l", axes = FALSE,
     ylab = "Trend", xlab = "")
box()
axis(2)
axis(1, at = seq(1, length(z), 12), labels = 1988:1993)

# Seasonal
plot(1:length(z), seasonal, type = "p", pch = 16, axes = FALSE,
     ylab = "Seasonal", xlab = "")
abline(h = 0)
segments(1:length(z), 0, 1:length(z), seasonal)
box()
axis(2)
axis(1, at = seq(1, length(z), 12), labels = 1988:1993)

# Irregular
plot(1:length(z), fits$resid, type = "p", pch = 16, axes = FALSE,
     ylab = "Irregular", xlab = "")
abline(h = 0)
segments(1:length(z), 0, 1:length(z), fits$resid)
box()
axis(2)
axis(1, at = seq(1, length(z), 12), labels = 1988:1993)

# Title
mtext("Decomposition of Snatching Data Using Grouped Seasonal Model",
      side = 3, outer = TRUE, cex = 1.5)

summary(fits)
dev.off()
plot(trendfit$resid,xlab="Month",ylab="residual if trendfit",type='b',pch=16,
     main="Detrended Monthly Snatching in Hyde Park")
plot(fits$residuals,
     xlab = "Month",
     ylab = "Residual of fits",
     type = "b",
     pch = 16,
     main = "Detrended and Deseasonalised Monthly Snatching in Hyde Park",
     sub = 'Quadratic trend and grouped seasonal')
abline(h=0)

# Box-Jenkins Identification for the irregular series
Ident<-function(x, data = NULL, sub = NULL)
{
  u <- acf(x, plot = F)
  v <- acf(x, type = "partial", plot = F)
  c1 <- 2./sqrt(u$n.used)
  oldpars <- par()
  oldpars$pin <- c(7.256, 5.216)
  par(mfrow = c(2., 1.), oma = c(3., 0., 4., 0.))
  plot(u$lag, u$acf, type = "b", lty = 3., ylab = "ACF", xlab = "Lag",
       ylim = range(u$acf, -c1, c1), main = 
         "Autocorrelation Function")
  segments(u$lag, 0., u$lag, u$acf)
  box(lty = 1.)
  abline(h = c(-c1, c1), lty = 2.)
  abline(h = 0., lty = 1.)
  plot(v$lag, v$acf, type = "b", lty = 2., ylab = "PACF", xlab = "Lag",
       ylim = range(v$acf, -c1, c1), main = 
         "Partial Autocorrelation Function")
  segments(v$lag, 0., v$lag, v$acf)
  box(lty = 1.)
  abline(h = c(-c1, c1), lty = 2.)
  abline(h = 0., lty = 1.)
  if(is.null(data))
    mtext("Box-Jenkins ARMA Model Identification", 3., 2., outer
          = T, cex = 1.5)
  else mtext(paste("Box-Jenkins ARMA Model Identification for", data),
             3., 2., outer = T, cex = 1.5)
  if(!is.null(sub))
    mtext(paste(sub), 1., 2., outer = T)
  invisible()
}
Ident(fits$resid)


.startval<-function (y, Z, tol)
{type <- FALSE     # hard-coded, TRUE as proposed in [MarT82b], [StoDut87]
  #             FALSE as proposed in [MarZ78]
  ar.ls <- lm.fit(Z, y, tol)
  phi <- ar.ls$coefficients
  if (type) {
    s <- sqrt(sum((ar.ls$residuals)^2)/ar.ls$df.residual)
  } else {
    s <- mad(ar.ls$residuals)
  }
  return(list(phi=phi, s=s))
}
psiHuber<-function (t, k=1.345, rho=FALSE) 
{
  if (!rho) {
    pmin(k, pmax(-k, t))
  } else {
    r <- abs(t)
    i <- r < k
    r[i] <- r[i]^2/2
    r[!i] <- k * (r[!i] - k/2)
    r
  }
}
.psi<-function (type)
{
 
  switch(type, 
         Ident = get("psiLS", mode="function"), 
         Huber = get("psiHuber", mode="function"), 
         Tukey = get("psiTukey", mode="function"),
         Hampel = get("psiHampel", mode="function")) 
}
.weights<-function (r, s, u, v, psi1, ...)
{psi <- .psi(psi1)
  n <- length(r)
  w <- rep(NA, n)
  for (i in 1:n) {
    if (r[i] == 0) {
      if (u[i] != 0) {
        w[i] <- v[i]/u[i]
      } else {
        w[i] <- 1
      }
    } else if (u[i] != 0) {
      dummy <- r[i]/s
      w[i] <- v[i]*psi(dummy/u[i], ...)/dummy
    } else if (psi1 == "Ident") {
      w[i] <- 1
    } else {
      w[i] <- 0
    }
  }
  return(w)
}
.Weights<-function (p, Z, invCp, type, psi2, c)
{
 
  psi <- .psi(psi2)
  d <- sqrt(diag(Z%*%invCp%*%t(Z))/p)
  if (psi2 == "Huber") {
    v <- psi(d, k=c)/d
  } else if (psi2 == "Tukey") {
    v <- psi(d, c=c)/d
  } else if (psi2 == "Ident") {
    v <- rep(1, nrow(Z)) 
  } else {
    warning("error in function '.Weights': psi function ", psi2, 
            "not defined \n")
  }
  if (type=="Mallows") {
    u <- rep(1, length(v))
  } else if (type=="Schweppe") { 
    u <- v
  } else {
    warning("error in function '.Weights': wrong GM-estimates type \n")
  }
  return(list(u=u, v=v))
}
.invCp<-function (p, s, Phi)
{
  
  if (p > 1) {
    M1 <- matrix(rep(rev(s), p), p)
    M2 <- cbind(rep(0, (p-1)), (-1)*Phi)
    M2 <- rbind(M2, rep(0, p))
    diag(M2) <- 1
    Ap <- (1/M1)*M2
    invCp <- t(Ap)%*%Ap 
  } else { 
    invCp <- 1/s[1]^2
  }
  return(invCp)
}
.ARmodel<-function (x, p)
{
  
  n <- length(x)
  y <- x[(p+1):n]
  Z <- embed(x, p)[1:(n-p), ]
  return(list(y=y, Z=as.matrix(Z)))
}
.BH<-function (k=1.345) 
{-2*k*dnorm(k) + 2*pnorm(k) + 2*k^2*(1-pnorm(k)) -1
}
.IWLS<-function (y, Z, phi.ini, s.ini, u, v, psi1, niter, tol, ...)
{stop1 <- sqrt(diag(solve(t(Z)%*%Z)))
  B <- NA
  phi <- phi.ini
  s.new <- s.ini
  iter <- 0
  psi <- .psi(psi1)
  if (psi1=="Huber") { 
    B <- .BH(...)
  } else { 
    B <- .BB(...)
  }
  
  while (iter < niter) {
    iter <- iter + 1
    r <- y - Z%*%phi 
    s.old <- s.new
    s2 <- (1/(B*sum(u*v)))*sum(u*v*(psi(r/(u*s.old), ...))^2)*s.old^2
    s.new <- sqrt(s2)
    w <- .weights(r, s.new, u, v, psi1, ...)
    ar.wls <- lm.wfit(Z, y, w, tol)
    tau <- ar.wls$coefficients - phi
    omega <- 1     # hard-coded, 0 < omega < 2, as proposed in [StoDut87]
    phi <- phi + omega*tau
    stop2 <- abs(c(s.old - s.new, omega*tau)) < tol*s.new*c(1, stop1)
    if (sum(stop2) == length(stop2)) break
  }
  return(list(phi=phi, s=s.new, w=w, B=B, niter=iter))
}
arGM<-function (x, order=1, 
          chr=1.5, iterh=maxiter, cbr=5.0, iterb=maxiter, 
          psi2="Tukey", c=4.0, type="Mallows", 
          k=1.5, maxiter=100, tol=1e-08, equal.LS=FALSE,...) 
{
  
  s <- c()
  Phi <- matrix()
  w <- NA
  BH <- BB <- NA
  niterh <- niterb <- niter.testing <- NA
  
  ##  Centering:
  
  ## x.huber <- HuberM(x, ...)     # as proposed in [StoDut87]
  x.huber <- hubers(x)     # as proposed in [MarZ78], [Mart80]
  x <- x - x.huber$mu
  sx <- x.huber$s
  
  ##  Main:
  
  for (p in 1:order) {
    ARmodel <- .ARmodel(x, p)
    y <- ARmodel$y
    Z <- ARmodel$Z
    invCp <- .invCp(p, c(sx, s), Phi)
    Weights <- .Weights(p, Z, invCp, type, psi2, c)
    u <- Weights$u
    v <- Weights$v
    startval <- .startval(y, Z, tol)
    phi <- startval$phi
    s[p] <- startval$s
    if (equal.LS) {     # for testing purpose only
      psi1 <- "Ident"
      niter <- maxiter
      IWLS <- .IWLS(y, Z, phi, s[p], u, v, psi1, niter, tol)
      phi <- IWLS$phi
      w <- IWLS$w
      niter.testing <- IWLS$niter
    } else {
      if ((iterh > 0) & (is.numeric(iterh))) {
        psi1 <- "Huber"
        niter <- iterh
        IWLS <- .IWLS(y, Z, phi, s[p], u, v, psi1, niter, tol, k=chr)
        phi <- IWLS$phi
        s[p] <- IWLS$s
        w <- IWLS$w
        BH <- IWLS$B
        niterh <- IWLS$niter
      }
      if ((iterb > 0) & (is.numeric(iterb))) {
        psi1 <- "Tukey"
        niter <- iterb
        IWLS <- .IWLS(y, Z, phi, s[p], u, v, psi1, niter, tol, c=cbr)
        phi <- IWLS$phi
        s[p] <- IWLS$s
        w <- IWLS$w
        BB <- IWLS$B
        niterb <- IWLS$niter
      }
    }
    if (p > 1) {
      Phi <- cbind(rep(0, (p-1)), Phi)
      Phi <- rbind(phi, Phi)
    } else {
      Phi <- as.matrix(phi)
    }
  }
  Cx <- solve(invCp)
  
  return(list(ar=phi, sinnov=s, Cx=Cx, mu=x.huber$mu, sx=sx, u=u, v=v, w=w, 
              BH=BH, BB=BB, 
              niterh=niterh, niterb=niterb, niter.testing=niter.testing))
  
}
Raic<-function(x, order = 10.)
{
  n <- length(x)
  z <- matrix(0., n - 1., order)
  for(k in 1.:order)
    z[(order + 1. - k):(n - 1.), k] <- x[1.:(n - order - 1. + k)]
  res <- matrix(0., n - 1., order)
  mu <- rep(0., order)
  svec <- mu
  ar <- matrix(0., order, order)
  y <- arGM(x, order, chr = 1.345, psi2="Huber", iterb=0)
  mu[order] <- y$mu
  svec[order] <- y$sinnov[order]
  ar[, order] <- y$ar
  s <- y$sinnov[order]
  res[, order] <- x[2.:n] - rep(y$mu, n - 1.) - z[, order:1.] %*% y$ar
  for(k in 1.:(order - 1.)) {
    y <- arGM(x, k,  chr = 1.345, psi2="Huber", iterb=0)
    mu[k] <- y$mu
    svec[k] <- y$sinnov[k]
    ar[1.:k, k] <- y$ar
    s <- y$sinnov[k]
    res[, k] <- x[2.:n] - rep(y$mu, n - 1.) - as.matrix(z[, order:
                                                            (order - k + 1.)]) %*% as.vector(y$ar)
  }
  list(k = c(1.:order), resid = res, coef = ar, mu = mu,
       sd = svec)
}
fin <- Raic(fits$resid)

res <- fin$resid[,2]
fv <- fits$resid[-2] - res

Ident(res)


#dignostics plot
par(mfcol=c(2,2),oma = c(6,0,6,0))

plot(1:(length(z)-1),res,type="b",
     main="Residual Series",
     xlab="Time",
     ylab="Residuals")

qqnorm(res,
       main="Quantile-Quantile Plot",
       xlab="Gaussian Quantiles",
       ylab="Residuals")
qqline(res);abline(h=0,v=0,col="grey")

plot(fv,res,
     main="Residual Plot",
     xlab="Lagged Irregular Values",
     ylab="Residuals")

plot(fv,abs(res),
     main="Absolute Residual Plot",
     xlab="Lagged Irregular Values",
     ylab="Absolute Residuals")
lines(lowess(fv,abs(res)))

mtext("Diagnostics for the Snatching data",
      side=3,line=2,outer=T,cex=1.5)
mtext("Markov model fitted after removing quadratic trend and seasonal trend 
      ",
      side=1,line=2,outer=T,cex=1.5)
fin$coef[1,1]
fin$coef[2,2]
summary(fits)
