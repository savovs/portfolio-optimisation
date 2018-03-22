library(PortfolioAnalytics)
library(GA)
library(quantmod)
library(timeSeries)
library(sm)
library(corrplot)

# Some stock symbols inspired by Vanguard Fund
# https://www.vanguardinvestor.co.uk/investments/vanguard-global-balanced-gbp-accumulation-shares?intcmpgn=blendedglobal_globalbalancedfund_fund_link
symbols <- c("AAPL", "GOOGL", "MSFT", "JNJ", "WMT", "CVX", "PG", "WFC", "INTC", "UPS", "PNC", "BAC", "CSCO")
getSymbols(symbols, src = "yahoo", from = "2013-01-01", to = "2018-01-01")

returns <- lapply(
  symbols,
  function(s) monthlyReturn(eval(parse(text = s)))
)

returns <- do.call(cbind, returns)
returns <- na.omit(returns, 0)
colnames(returns) <- symbols

# Train / test split
trainReturns <- head(returns, nrow(returns) - 12)
testReturns <- tail(returns, 12)

nStocks <- ncol(trainReturns)
meanReturns <- colMeans(trainReturns)
covarianceMatrix <- cov(trainReturns)
volatility <- sqrt(diag(covarianceMatrix))

# Explore the Data
plot(volatility, meanReturns, type = "n", panel.first = grid(),
     xlab = "Std. dev. monthly returns", ylab = "Average monthly returns")

text(volatility, meanReturns, names(meanReturns), col = .colorwheelPalette(10), font = 2)

corrplot(cor(returns), method = "color", title = "Correlation of Monthly Returns",  mar = c(0,0,1,0))


# Set up and run GA
normaliseWeights <- function(w) {
  drop(w / sum(w))
}

expectedReturn <- function(w) {
  sum(normaliseWeights(w) * meanReturns)
}

portfolioVariance <- function(w) {
  w <- normaliseWeights(w)
  drop(w %*% covarianceMatrix %*% w)
}

fitness1 <- function(w) expectedReturn(w) ^ 2 / portfolioVariance(w)
fitness2 <- function(w) {
  ret <- expectedReturn(w) - 0.01

  penalty <- if(ret < 0) ret ^ 2 else 0
  
  return(-(portfolioVariance(w) + penalty))
}

GA <- ga(
  type = "real-valued",
  fitness = fitness1,
  min = rep(0, nStocks),
  max = rep(1, nStocks),
  names = symbols,
  maxiter = 1000,
  run = 200,
  optim = TRUE,
  popSize = 50
)

gaWeights <- normaliseWeights(GA@solution)
gaExpectedReturn <- expectedReturn(gaWeights)
gaVariance <- portfolioVariance(gaWeights)

barplot(gaWeights, xlab = "Stocks", ylab = "Portfolio weights",
        cex.names = 0.7, col = .colorwheelPalette(10))

plot(GA)

# Generate random weights for comparison
randomNumbers <- runif(length(symbols))
randomSum <- sum(randomNumbers)
randomWeights <- sapply(randomNumbers, function(x) x / randomSum)

#Generate weights using PortfolioAnalytics https://cran.r-project.org/web/packages/PortfolioAnalytics/vignettes/portfolio_vignette.pdf
port <- portfolio.spec(assets = symbols)
port <- add.constraint(port, type = "box", min = 0.05, max = 0.8)
port <- add.constraint(portfolio = port, type = "full_investment")
solver <- random_portfolios(port, permutations = 5000, rp_method = "sample")

# Get minimum variance portfolio and optimization
minvar.port <- add.objective(port, type = "risk", name = "var")
minvar.opt <- optimize.portfolio(trainReturns, minvar.port, optimize_method = "random", rp = solver)

# Maximum return portfolio and optimization
maxret.port <- add.objective(port, type = "return", name = "mean")
maxret.opt <- optimize.portfolio(trainReturns, maxret.port, optimize_method = "random", rp = solver)

# Get blended weights between minimized variance / maximized return weights
combinedOptimizations <- combine.optimizations(list(minVar = minvar.opt, maxRet = maxret.opt))
paWeights <- colMeans(extractWeights(combinedOptimizations))


# Compare weights from GA, PortfolioAnalytics and random
gaReturns <- sweep(testReturns, 2, gaWeights, '*')
paReturns <- sweep(testReturns, 2, paWeights, '*')
randomReturns <- sweep(testReturns, 2, randomWeights, '*')

rmse <- function(error) sqrt(mean(error ^ 2))

errors <- list(
  random = rmse(gaReturns - randomReturns),
  pa = rmse(gaReturns - paReturns)
)

print("RMSE between GA returns and two other returns")
print(errors)

sums <- c(sum(gaReturns), sum(paReturns), sum(randomReturns))

barplot(
  sums,
  main = sprintf("Sum of returns of portfolios, std: %#.03f", sd(sums)),
  xlab = "Methods",
  ylab = "Sum of Returns",
  col = 2:5
)

legend("bottomright", c("Genetic Algorithm", "PortfolioAnalytics", "Random"), fill = 2:5)



