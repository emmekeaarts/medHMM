## Sensibility for Chi-Square prior

library(extraDistr)

# Emiss:
x <- rinvchisq(n = 1e6, nu = 50, tau = 10)

plot(density(x))

median(x)
sd(x)
range(x)
IQR(x)

# Dwell:
x <- rinvchisq(n = 1e6, nu = 3, tau = 1)

plot(density(x))

median(x)
sd(x)
range(x)
IQR(x)


# normal sd emiss:
x <- extraDistr::rinvgamma(n = 1e6, alpha = 0.01, beta = 0.01)

x <- MCMCpack::rinvgamma(n = 1e6, shape = 10, scale = 1)

# x <- rgamma(n = 1e6, shape = 0.1, rate = 0.1)
# x <- rgamma(n = 1e6, shape = 0.1, scale = 0.1)

plot(density(x))

median(x)
sd(x)
range(x)
IQR(x)




##

# Define the means and variances
means <- c(10, 17.38, 24.76)
means <- c(10,31.47,52.94)
variances <- c(30, 30, 30)
sds <- sqrt(variances)  # Standard deviations

# Set up the x range for plotting
x_range <- seq(min(means) - 3 * max(sds), max(means) + 3 * max(sds), length.out = 1000)

# Set up the plotting area
plot(x_range, dnorm(x_range, mean = means[1], sd = sds[1]), type = "n",
     xlab = "X", ylab = "Density", main = "Density Plot of Three Normal Distributions")

# Add the normal distribution curves
colors <- c("red", "blue", "green")

for (i in 1:3) {
    curve(dnorm(x, mean = means[i], sd = sds[i]), add = TRUE, col = colors[i], lwd = 2)
}


