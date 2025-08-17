library(ggplot2)
library(tidyr)

r2 <- 0.1
K <- 0.01
N <- 1:10
p <- 0.4
# For binomial:
# n_0 * p /(1-p) = N, so n_0 = (1-p) * N / p
n0 <- (1-p)*N / p
# For poisson, we solve
# lambda * p /(1-exp(-lambda*p)) = N
lambda <- sapply(N, function(N) 
  optimise(function(lambda) (lambda*p/(1-exp(-lambda*p)) - N)^2, 
         c(0.001, 100))$minimum)
result <- matrix(nrow=length(N), ncol=3)
for (n in N) {
  result[n, 1] <- 100*risk_reduction_lowest(r2, K, n)
  result[n, 2] <- 100*risk_reduction_lowest_bin(r2, K, n0[n], p)
  result[n, 3] <- 100*risk_reduction_lowest_pois(r2, K, lambda[n] * p)
}

result <- data.frame(result) 
names(result) <- c("Fixed live births", "Binomial", "Poisson")
result$n <- N
result_long <- pivot_longer(result, cols = 1:3)
result_long <- result_long |>
  subset(n > 1)

result_long$name <- factor(result_long$name, 
                           levels = c("Fixed live births", "Binomial", "Poisson"))


ggplot(data=result_long, aes(n, value, col=name)) +
  geom_point(size=3) +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.title = element_text(size=20),
        axis.text = element_text(size=20),
        legend.text = element_text(size=20),
        legend.position = "bottom") +
  labs(x = "Expected number of live births",
       y = "Relative risk reduction") +
  ylim(c(0, max(result_long$value)))

# ggsave("~/Desktop/figure.png", 
#        width=8, height=5)
