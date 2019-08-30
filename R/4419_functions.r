
inspect <- function(data, x, ...)UseMethod("inspect")
inspect.tbl_df <- function(data, x, ...){
  tmp <- data[[as.character(x)]]
  var.lab <- attr(tmp, "label")
  if(is.null(var.lab)){var.lab <- "No Label Found"}
  val.labs <- attr(tmp, "labels")
  if(is.null(val.labs)){val.labs <- sort(unique(tmp))}
  tab <- cbind(freq = table(tmp), prop = round(table(tmp)/sum(table(tmp), na.rm=T), 3))
  out <- list(variable_label = var.lab, value_labels=t(t(val.labs)), freq_dist = tab)
  return(out)
}
inspect.data.frame <- function(data, x, ...){
  var.lab <- NULL
  var.lab <- attr(data, "var.label")[which(names(data) == x)]
  if(is.null(var.lab) & "label" %in% names(attributes(data[[x]]))){
    var.lab <- attr(data[[x]], "label")
  }
  if(is.null(var.lab)){var.lab <- "No Label Found"}
  if("labels" %in% names(attributes(data[[x]])))
    val.labs <- attr(data[[x]], "labels")
  else{
  val.labs <- {if(!is.null(levels(data[[x]]))){levels(data[[x]])}
    else {sort(unique(data[[x]]))}
  }
  }
  tab <- cbind(freq = table(data[[x]]), prop = round(table(data[[x]])/sum(table(data[[x]]), na.rm=T), 3))
  out <- list(variable_label = var.lab, value_labels=t(t(val.labs)), freq_dist = tab)
  return(out)
}

## concordant, discordant, tau.b, tau.c, ord.somers.d, ord.gamma come from the ryouready package
## Phi and V come from the DescTools package
concordant <- function (x) {
    x <- matrix(as.numeric(x), dim(x))
    mat.lr <- function(r, c) {
        lr <- x[(r.x > r) & (c.x > c)]
        sum(lr)
    }
    r.x <- row(x)
    c.x <- col(x)
    sum(x * mapply(mat.lr, r = r.x, c = c.x))
}
discordant <- function(x){
    x <- matrix(as.numeric(x), dim(x))
    mat.ll <- function(r, c) {
        ll <- x[(r.x > r) & (c.x < c)]
        sum(ll)
    }
    r.x <- row(x)
    c.x <- col(x)
    sum(x * mapply(mat.ll, r = r.x, c = c.x))
}

tau.b <- function (x) {
    x <- matrix(as.numeric(x), dim(x))
    c <- concordant(x)
    d <- discordant(x)
    n <- sum(x)
    SumR <- rowSums(x)
    SumC <- colSums(x)
    tau.b <- (2 * (c - d))/sqrt(((n^2) - (sum(SumR^2))) * ((n^2) -
        (sum(SumC^2))))
    tau.b
}

ord.gamma <- function(x){
    x <- matrix(as.numeric(x), dim(x))
    c <- concordant(x)
    d <- discordant(x)
    gamma <- (c - d)/(c + d)
    class(gamma) <- "ord.gamma"
    gamma
}

ord.somers.d <- function(x){
    x <- matrix(as.numeric(x), dim(x))
    c <- concordant(x)
    d <- discordant(x)
    n <- sum(x)
    SumR <- rowSums(x)
    SumC <- colSums(x)
    sd.cr <- (2 * (c - d))/((n^2) - (sum(SumR^2)))
    sd.rc <- (2 * (c - d))/((n^2) - (sum(SumC^2)))
    sd.s <- (2 * (c - d))/((n^2) - (((sum(SumR^2)) + (sum(SumC^2)))/2))
    res <- list(sd.cr, sd.rc, sd.s)
    names(res) <- c("sd.cr", "sd.rc", "sd.symmetric")
    class(res) <- "ord.somersd"
    res
}

lambda <- function(x){
  wmax <- apply(x, 2, which.max)
  wgmax <- which.max(rowSums(x))
  nullcc <- rowSums(x)[wgmax]
  nullerr <- sum(rowSums(x)[-wgmax])
  corrpred <- x[cbind(wmax, 1:ncol(x))]
  errpred <- colSums(x) - corrpred
  E1 <- nullerr
  E2 <- sum(errpred)
  (E1-E2)/E1
}

phi <- function(x){
    num <- prod(diag(x))- (x[2,1]*x[1,2])
    denom <- sqrt(prod(c(colSums(x), rowSums(x))))
    num/denom
}

V <- function(x){
  if(all(dim(x) == 2)){
    num <- prod(diag(x))- (x[2,1]*x[1,2])
    denom <- sqrt(prod(c(colSums(x), rowSums(x))))
    num/denom
  }
  else{
  chi2 <- chisq.test(x, correct=FALSE)$statistic
  sqrt(chi2/(sum(c(x)) * (min(nrow(x), ncol(x)) -1)))
  }
}

simtable <- function(x,y, n=1000, stat=NULL){
  out <- lapply(1:n, function(i)table(x, sample(y, length(y), replace=F)))
  if(is.null(stat)){
    return(out)
  }
  else{
    sapply(out, stat)
  }

}

simrho <- function(x,y, n=1000){
  rho0 <- cor(x,y, use="pair", method="spearman")
  simrho <- sapply(1:n, function(i)cor(x, sample(y, length(y), replace=F), use="pair", method="spearman"))
  pv <- {if(rho0 >= 0)mean(simrho > rho0)
      else mean(simrho < rho0)}
  return(list(rho0 = rho0, simrho = simrho, pv = pv))
}

makeStats <- function(x,y, chisq=FALSE, phi=FALSE, cramersV=FALSE, lambda=FALSE,
   gamma=FALSE, d=FALSE, taub=FALSE, rho=FALSE, n=1000){

  tabs <- simtable(x,y,n)
  tab <- table(x,y)
allStats <- NULL
if(chisq){
  stat0 <- do.call('chisq.test', list(x=tab, correct=FALSE))$statistic
  stats <- sapply(tabs, function(x)chisq.test(x, correct=FALSE)$statistic)
  pv <- mean(stats > stat0)
  allStats <- rbind(allStats, c(stat0, pv))[,,drop=F]
  rownames(allStats)[nrow(allStats)] <- "Chi-squared"
}
if(phi){
  stat0 <- do.call('phi', list(x=tab))
  stats <- sapply(tabs, function(x)phi(x))
  pv <- mean(stats > stat0)
  allStats <- rbind(allStats, c(stat0, pv))[,,drop=F]
  rownames(allStats)[nrow(allStats)] <- "Phi"
}
if(cramersV){
  stat0 <- do.call('V', list(x=tab))
  stats <- sapply(tabs, function(x)V(x))
  pv <- mean(stats > stat0)
  allStats <- rbind(allStats, c(stat0, pv))[,,drop=F]
  rownames(allStats)[nrow(allStats)] <- "Cramers V"
}
if(lambda){
  stat0 <- do.call('lambda', list(x=tab))
  stats <- sapply(tabs, function(x)lambda(x))
  pv <- mean(stats > stat0)
  allStats <- rbind(allStats, c(stat0, pv))[,,drop=F]
  rownames(allStats)[nrow(allStats)] <- "Lambda"
}
if(gamma){
  stat0 <- do.call('ord.gamma', list(x=tab))
  stats <- sapply(tabs, function(x)lambda(x))
  pv <- {if(stat0 >= 0)mean(stats > stat0)
      else mean(stats < stat0)}
    allStats <- rbind(allStats, c(stat0, pv))[,,drop=F]
  rownames(allStats)[nrow(allStats)] <- "Kruskal-Goodman Gamma"
}
if(d){
  stat0 <- do.call('ord.somers.d', list(x=tab))$sd.symmetric
  stats <- sapply(tabs, function(x)ord.somers.d(x)$sd.symmetric)
  pv <- {if(stat0 >= 0)mean(stats > stat0)
      else mean(stats < stat0)}
  allStats <- rbind(allStats, c(stat0, pv))[,,drop=F]
  rownames(allStats)[nrow(allStats)] <- "Somers D"
}
if(taub){
  stat0 <- do.call('tau.b', list(x=tab))
  stats <- sapply(tabs, function(x)tau.b(x))
  pv <- {if(stat0 >= 0)mean(stats > stat0)
      else mean(stats < stat0)}
  allStats <- rbind(allStats, c(stat0, pv))[,,drop=F]
  rownames(allStats)[nrow(allStats)] <- "Tau-b"
}
if(rho){
  x2 <-as.numeric(x)
  y2 <- as.numeric(y)
  r <- simrho(x2,y2,n)
  allStats <- rbind(allStats, c(r$rho0, r$pv))[,,drop=F]
  rownames(allStats)[nrow(allStats)] <- "Spearmans Rho"
}
if(!is.null(allStats)){
  colnames(allStats) <- c("statistic", "p-value")
  w <- which(allStats[,1] == 0 & allStats[,2] == 0)
  if(length(w) > 0){
    allStats[w,2] <- 1.000
  }
  allStats <- round(allStats, 4)
}
return(allStats)
}

plotStdRes <- function(x, col=RColorBrewer::brewer.pal(10, "RdBu"), ...){
  x2 <- chisq.test(x)
  res <- x2$stdres
  minres <- ifelse(min(c(res)) > -3.5, -3.5, min(c(res)))
  maxres <- ifelse(max(c(res)) <3.5, 3.5, max(c(res)))
  if(maxres > abs(minres)){
    minres <- -maxres
  }
  if(maxres < abs(minres)){
    maxres <- -minres
  }
  lattice::levelplot(res, col.regions=col, cuts=10, at=c(minres, -3, -2, -1, 0, 1, 2, 3, maxres), ...)
}

pmc <- function(obj, col=RColorBrewer::brewer.pal(6, "RdBu"), ...){
  X <- model.matrix(obj)[,-1]
  R <- cor(X)
  diag(R) <- 0
  return(lattice::levelplot(R, at=c(-1, -.75, -.5, 0, .5, .75, 1), col.regions=col, ..., xlab="", ylab="", scales=list(x=list(rot=90))))
}

plotCIgroup <- function(form, data, includeOverall=TRUE, ...){
    cfun <- function(x, ...){tmp <- confidenceInterval(x, ...); data.frame(y = tmp[1], ymin = tmp[2], ymax=tmp[3])}
    dot.args <- as.list(match.call(expand.dots = FALSE)$`...`)
    mf <- model.frame(form, data)
    nums <- sort(unique(as.numeric(mf[[2]])))
    if(includeOverall){
      mf2 <- mf
      mf2[[2]] <- factor(max(nums)+1, levels=c(nums, max(nums)+1), labels=c(levels(mf[[2]]), "Overall"))
      levels(mf[[2]]) <- c(levels(mf[[2]]), "Overall")
      mf <- rbind(mf, mf2)
    }
    
    ggplot(mf, aes_string(x=names(mf)[2], y=names(mf)[1])) + stat_summary(fun.data=cfun, fun.args = dot.args)
}

searchVarLabels <- function(dat, str) UseMethod("searchVarLabels")
searchVarLabels.data.frame <-
function (dat, str)
{
    vlat <- NULL
    if ("var.labels" %in% names(attributes(dat))) {
        vlat <- "var.labels"
      ind <- sort(union(grep(str, attr(dat, vlat), ignore.case = T), grep(str, names(dat), ignore.case = T)))
      labs <- attr(dat, vlat)
    }
    if ("variable.labels" %in% names(attributes(dat))) {
        vlat <- "variable.labels"
      ind <- sort(union(grep(str, attr(dat, vlat), ignore.case = T), grep(str, names(dat), ignore.case = T)))
      labs <- attr(dat, vlat)
    }
    natt <- sapply(1:ncol(dat), function(i)names(attributes(dat[[i]])))
    natt1 <- unique(c(unlist(natt)))
    if("label" %in% natt1){
      haslabs <- sapply(natt, function(x)"label" %in% x)
      labs <- vlat <- sapply(1:length(haslabs), function(x)ifelse(haslabs[x], attr(dat[[x]], "label"), ""))
      ind <- sort(union(grep(str, labs, ignore.case = T), grep(str, names(dat), ignore.case = T)))
  }
    if(is.null(vlat))stop("No Labels to Search")
    vldf <- data.frame(ind = ind, label = labs[ind])
    rownames(vldf) <- names(dat)[ind]
    vldf
}
searchVarLabels.tbl_df <-
function (dat, str)
{
    vlat <- unlist(sapply(1:ncol(dat), function(i)attr(dat[[i]], "label")))
    ind <- sort(union(grep(str, vlat, ignore.case = T), grep(str, names(dat), ignore.case = T)))
    vldf <- data.frame(ind = ind, label = vlat[ind])
    rownames(vldf) <- names(dat)[ind]
    vldf
}

freqDist <- function(x){
  tab <- table(x)
  ntab <- names(tab)
  pct <- tab/sum(tab)*100
  cpct <- cumsum(pct)
  tab <- c(tab, sum(tab))
  names(tab) <- c(ntab, "Total")
  pct <- c(pct, 100)
  cpct <- c(cpct, 100)
  maxnum <- max(nchar(tab))
  fr <- sprintf(paste0("%", maxnum, ".0f"), tab)
  pc <- sprintf("%6.2f", pct)
  cp <- sprintf("%6.2f", cpct)
  cp[length(cp)] <- ""
  out <- cbind(fr, pc, cp)
  rownames(out) <- names(tab)
  colnames(out) <- c("Freq", "  %   ", " Cu % ")
  noquote(out)
}

histDiscrete <- function(x, data, ...){
  m <- min(data[[x]], na.rm=TRUE)    
  ggplot(data, aes_string(x=x)) + geom_histogram(binwidth=1, center=m, color="black", fill="gray75")
}

unalike <- function(x){
  o <- outer(x, x, "!=")
  mean(c(o[lower.tri(o)]), na.rm=T)
}

sig.cor <- function(x,y, method=c("z", "t", "sim"), n.sim = 1000, two.sided=TRUE, ...){
meth <- match.arg(method)
r <- cor(x,y, use="pairwise.complete.obs", ...)
n <- sum(!is.na(x)*!is.na(y))
if(meth == "z"){
  z <- .5*log((1+r)/(1-r))
  sez <- 1/sqrt(n-3)
  pv <- (2^two.sided)*pnorm(abs(z), 0, sez, lower.tail=F)
}
if(meth == "t"){
  tstat <- r*sqrt((n-2)/(1-r^2))
  pv <- (2^two.sided)*pt(abs(tstat), n-2, lower.tail=F)
}
if(meth == "sim"){
  xmat <- sapply(1:n.sim, function(z)sample(x, length(x), replace=F))
  r0 <- c(cor(y, xmat))
  pv <- {if(two.sided){
    mean(r0 < (-abs(r))) + mean(r0 > abs(r))
  }
  else{
    if(r > 0){
      mean(r > r0)
    }
    else{
      mean(r < r0)
    }
  }}
}
return(list(r=r, p = pv))
}


pwCorrMat <- function(X, method=c("z", "t", "sim"), ...){
  meth <- match.arg(method)
  out <- p.out <- diag(ncol(X))
  for(i in 1:(ncol(X)-1)){
    for(j in i:ncol(X)){
      f <- sig.cor(X[,i], X[,j], method=meth, ...)
      out[i,j] <- out[j,i] <- f$r
      p.out[i,j] <- p.out[j,i] <- f$p
    }
  }
  outSig <- matrix(sprintf("%.3f", out), ncol=ncol(X))
  outSig[which(p.out > .05, arr.ind=T)] <- ""
  diag(outSig) <- ""
  outSig[upper.tri(outSig)] <- ""
  colnames(outSig) <- colnames(out) <- rownames(outSig) <- rownames(out) <- colnames(p.out) <- rownames(p.out) <- colnames(X)
  ret <- list(rSig=outSig, r=out, p = p.out )
  class(ret) <- "pwc"
  return(ret)
}

print.pwc <- function(x, ...){
cat("All Correlations\n")
 print(noquote(round(x$r,3)))
cat("\nOnly Significant Correlations\n")
 print(noquote(x$rSig))
}


GKGamma <- function (x, y = NULL, conf.level = NA, ...){
## Function taken from DescTools v0.99.22
    if (!is.null(y))
        tab <- table(x, y, ...)
    else tab <- as.table(x)
    x <- ConDisPairs(tab)
    psi <- 2 * (x$D * x$pi.c - x$C * x$pi.d)/(x$C + x$D)^2
    sigma2 <- sum(tab * psi^2) - sum(tab * psi)^2
    gamma <- (x$C - x$D)/(x$C + x$D)
    if (is.na(conf.level)) {
        result <- gamma
    }
    else {
        pr2 <- 1 - (1 - conf.level)/2
        ci <- qnorm(pr2) * sqrt(sigma2) * c(-1, 1) + gamma
        result <- c(gamma = gamma, lwr.ci = max(ci[1], -1), ups.ci = min(ci[2],
            1))
    }
    class(result) <- "gkg"
    return(result)
}

confidenceInterval <- function (x, confidence = 0.95,  na.rm = TRUE, distr=c("normal", "t")){
    distr <- match.arg(distr)
    nobs <- sum(!is.na(x))
    est <- mean(x, na.rm = na.rm)
    stderr <- sd(x, na.rm = na.rm)/sqrt(nobs)
    alpha <- 1-confidence
    if(distr == "t"){
      ci.low <- est + qt(alpha/2, nobs - 1) * stderr
      ci.high <- est - qt(alpha/2, nobs - 1) * stderr
    }
    else{
      ci.low <- est + qnorm(alpha/2) * stderr
      ci.high <- est - qnorm(alpha/2) * stderr
    }
    retval <- c(Estimate = est, `CI lower` = ci.low, `CI upper` = ci.high,
        `Std. Error` = stderr)
    retval
}

print.gkg <- function(x, ...){
  if(class(x) != "gkg")stop("Object must be of class gkg\n")
  if(length(x) == 1){
    cat("Goodman-Kruskal's Gamma = ", round(x,3), "\n", sep="")
  }
  if(length(x) == 3){
    cat("Goodman-Kruskal's Gamma = ", round(x,3), ", 95% CI (", round(x[2], 3), ", ", round(x[3],3),  ")\n", sep="")
  }
}
print.ktb <- function(x, ...){
  if(class(x) != "ktb")stop("Object must be of class ktb\n")
  cat("Kendall's Tau-b = ", round(x,3), "\n", sep="")
}

barplotStats <- function(x, y, data, stat="sum", ...){
  dot.args <- as.list(match.call(expand.dots = FALSE)$`...`)
  ggplot(data, aes_string(x=x, y=y)) + stat_summary(fun.y = stat, geom="bar", fun.args=dot.args)
}

sumStats <- function(data, vars, byvar=NULL, convertFactors=FALSE){
  if(is.null(byvar)){
    X <- data[,vars, drop=FALSE]
    if(convertFactors){
      for(i in 1:ncol(X)){
        if(is.factor(X[[i]]))X[[i]] <- as.numeric(X[[i]])
      }
    }
    means <- colMeans(X, na.rm=T)
    sds <- apply(X, 2, sd, na.rm=T)
    qtiles <- t(apply(X, 2, quantile, probs = c(0,.25,.5,.75,1), na.rm=TRUE))[,,drop=FALSE]
    iqr <- qtiles[,4]-qtiles[,2]
    n <- apply(X, 2, function(x)sum(!is.na(x)))
    na <- apply(X, 2, function(x)sum(is.na(x)))
    out <- cbind(means, sds, iqr, qtiles, n, na)
    colnames(out) <- c("Mean", "SD", "IQR", "0%", "25%", "50%", "75%", "100%", "n", "NA")
  }
  else{
    unvals <- unique(na.omit(data[[byvar]]))
    out <- vector(mode="list", length=length(unvals))
    for(i in 1:length(unvals)){
      X <- data[which(data[[byvar]] == unvals[i]),vars, drop=FALSE]
      if(convertFactors){
        for(i in 1:ncol(X)){
          if(is.factor(X[[i]]))X[[i]] <- as.numeric(X[[i]])
        }
      }
      means <- colMeans(X, na.rm=T)
      sds <- apply(X, 2, sd, na.rm=T)
      qtiles <- t(apply(X, 2, quantile, probs = c(0,.25,.5,.75,1), na.rm=TRUE))[,,drop=FALSE]
      iqr <- qtiles[,4]-qtiles[,2]
      n <- apply(X, 2, function(x)sum(!is.na(x)))
      na <- apply(X, 2, function(x)sum(is.na(x)))
      out[[i]] <- cbind(means, sds, iqr, qtiles, n, na)
      colnames(out[[i]]) <- c("Mean", "SD", "IQR", "0%", "25%", "50%", "75%", "100%", "n", "NA")
      names(out)[[i]] <- paste(byvar, " = ", unvals[i], sep="")
    }
  }
  out
}

histNorm <- function(x, data, normCurve=TRUE, densCurve=FALSE, bins=30){
  s <- seq(min(data[[x]], na.rm=TRUE), max(data[[x]], na.rm=T), length=100)
  dn <- dnorm(s, mean(data[[x]], na.rm=T), sd(data[[x]], na.rm=TRUE))
  tmp <- data[[x]]
  dens <- density(na.omit(tmp))
  dens <- data.frame(x=dens$x, y=dens$y)
  g <- ggplot(data, aes_string(x=x)) + geom_histogram(aes(y=stat(density)), bins=bins, color="black", fill="gray70")
  if(normCurve){
    g <- g+stat_function(fun=dnorm, n=101, args=list(mean=mean(data[[x]], na.rm=TRUE), sd=sd(data[[x]], na.rm=TRUE)), color="blue")
  }
  if(densCurve){
    g <- g+geom_line(data=dens, aes_string(x="x",y="y"), color="red")
  }
  g
}

propci <- function(x, n=NULL, conf.level=.95){
  tailprob <- (1-conf.level)/2
  ltp <- tailprob
  utp <- 1-tailprob
  nobs <- ifelse(is.null(n), sum(!is.na(x)), n)
    if(length(x) > 1){
    p <- mean(x, na.rm=T)
    nx <- sum(x == 1, na.rm=TRUE)
  }
  else{
    nx <- x
    p <- nx/nobs
  }
  zcrit <- abs(qnorm((1-conf.level)/2))
  norm.ci <- p + c(-1,1)*zcrit * sqrt((p*(1-p))/nobs)
  binom.ci <- c(qbeta(ltp, nx, nobs - nx + 1), qbeta(utp, nx + 1, nobs - nx))
  out <- rbind(norm.ci, binom.ci)
  colnames(out) <- c("Lower", "Upper")
  rownames(out) <- c("Normal Approx", "Exact")
  return(out)
}

tTest <- function(x,y, data, ...){
  formula <- as.formula(paste(y, x, sep=" ~ "))
  tmp <- get_all_vars(formula, data)
  tmp <- na.omit(tmp)
  if(is.factor(tmp[[x]]))tmp[[x]] <- droplevels(tmp[[x]])
  g <- vector(mode="list", length=3)
  ng <- levels(tmp[[x]])
  if(is.null(ng)){
    ng <- sort(unique(ng))
  }
  tt <- t.test(formula, data, ...)
  names(g) <- c(ng, "Difference")
  g[[1]]$mean <- mean(tmp[[y]][which(tmp[[x]] == ng[1])], na.rm=TRUE)
  g[[1]]$n <- sum(tmp[[x]] == ng[1])
  g[[1]]$se <- sd(tmp[[y]][which(tmp[[x]] == ng[1])], na.rm=TRUE)
  g[[2]]$mean <- mean(tmp[[y]][which(tmp[[x]] == ng[2])], na.rm=TRUE)
  g[[2]]$n <- sum(tmp[[x]] == ng[2])
  g[[2]]$se <- sd(tmp[[y]][which(tmp[[x]] == ng[2])], na.rm=TRUE)
  g[[3]]$mean <- g[[1]]$mean - g[[2]]$mean
  g[[3]]$n <- g[[1]]$n+g[[2]]$n
  g[[3]]$se <- tt$stderr
  res <- list(sum=do.call(rbind, g), tt=tt)
  class(res) <- "tTest"
  print(res)
}

print.tTest <- function(x, ...){
  cat("Summary:\n")
  print(x$sum)
  cat("\n")
  print(x$tt)
}
