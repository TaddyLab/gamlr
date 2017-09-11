beta <- 1:100
y <- log ( 1 + 2* abs(beta))
plot(beta, y)

plot(beta, abs(beta))

plot(beta, beta^2)

y <- 1/( 1 + 2* abs(beta))
plot(beta, y)

y <- (-2)/( 1 + 2* abs(beta))^2
plot(beta, y)

beta <- seq(0, 10, length.out = 1000)
w <- -log(0.9*dnorm(beta, 0, 0.00001) + 0.025*dnorm(beta, 0, 0.5) + 0.025*dnorm(beta, 0, 1) 
         + 0.025*dnorm(beta, 0, 2) + 0.025*dnorm(beta, 0, 5))
plot(beta, w)


plot(beta, beta^2)

