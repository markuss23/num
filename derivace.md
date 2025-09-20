<img width="591" height="491" alt="image" src="https://github.com/user-attachments/assets/b4c957b4-738c-41e1-9dae-472e75e7b214" />


```r

# 00_settings.R -----------------------------------------------------------
rm(list=ls()); cat("\014")
set.seed(123)
options(digits = 8)

# Pomocné utilitky
check_tol <- function(x) format(x, scientific=TRUE, digits=3)
timer <- function(expr) { t0 <- proc.time(); on.exit(print(proc.time()-t0)); force(expr) }

# 01_problem.R ------------------------------------------------------------
# Definuj testovací funkci a „ground truth“ (když jde)
f  <- function(x) sin(x)
df <- function(x) cos(x)           # přesná derivace, když máme
Fint <- function(a,b) -cos(b)+cos(a) # přesný integrál ∫ sin = -cos

# 02_methods.R ------------------------------------------------------------
# Derivace
ForwardDiff <- function(f,x,h=1e-6) (f(x+h)-f(x))/h
CentralDiff <- function(f,x,h=1e-6) (f(x+h)-f(x-h))/(2*h)

# Integrály
MidpointRule <- function(f,a,b,n=1000) {
  h <- (b-a)/n
  h*sum(f(a + h*(1:n) - h/2))
}
SimpsonRule <- function(f,a,b,n=50){ # n = počet dílků
  h <- (b-a)/n
  hpul <- h/2
  suma <- f(a)+f(b)
  xs <- h*(1:n) + a - hpul
  suma <- suma + 4*sum(f(xs))
  if(n>1){ xl <- (xs + hpul)[-n]; suma <- suma + 2*sum(f(xl)) }
  suma*h/6
}

# Kořeny
BisecRoot <- function(f,a,b,tol=1e-8, maxit=1e5){
  fa <- f(a); fb <- f(b); if(fa*fb >= 0) stop("Need sign change on [a,b]")
  it <- 0
  repeat{
    c <- (a+b)/2; fc <- f(c); it <- it+1
    if(abs(fc)<tol || (b-a)/2 < tol || it>=maxit) return(c)
    if(fa*fc < 0){ b <- c; fb <- fc } else { a <- c; fa <- fc }
  }
}
NewtonRoot <- function(f,x0,tol=1e-8, maxit=1e5){
  h <- sqrt(tol); x <- x0
  for(it in 1:maxit){
    d <- (f(x+h)-f(x-h))/(2*h)
    if(abs(d) < .Machine$double.eps) stop("Derivative ~0")
    dx <- f(x)/d
    x_new <- x - dx
    if(abs(dx) < tol) return(x_new)
    x <- x_new
  }
  stop("No convergence")
}

# ODE
EulerStep <- function(f,x,y,h) y + h*f(x,y)
RK4 <- function(f,x,y,h){
  hh <- 0.5*h
  k1 <- f(x,y)
  k2 <- f(x+hh, y+hh*k1)
  k3 <- f(x+hh, y+hh*k2)
  k4 <- f(x+h , y+h *k3)
  y + h*(k1 + 2*(k2+k3) + k4)/6
}

# Interpolace
Lagrange <- function(xa,x,y){
  n <- length(x); S <- 0
  for(i in 1:n){
    L <- 1
    for(j in 1:n) if(j!=i) L <- L * (xa - x[j])/(x[i]-x[j])
    S <- S + y[i]*L
  }
  S
}

# Polynomy
Horner <- function(a,x){
  n <- length(a); y <- a[n]
  for(i in (n-1):1) y <- y*x + a[i]
  y
}

# 03_experiments.R --------------------------------------------------------
# A) Derivace test (konvergenční studie)
x0 <- 1
for(h in c(1e-1,1e-2,1e-3,1e-4,1e-5)){
  num <- CentralDiff(f, x0, h)
  err <- abs(num - df(x0))
  cat("h =", h, " | derivace ~", num, " | chyba =", check_tol(err), "\n")
}

# B) Integrál test
a <- 0; b <- pi
I_true <- Fint(a,b)
for(n in c(10, 50, 200, 1000)){
  I_mid <- MidpointRule(f,a,b,n)
  I_simp <- SimpsonRule(f,a,b,n)
  cat("n =", n, "| Mid:", check_tol(abs(I_mid-I_true)), "| Simp:", check_tol(abs(I_simp-I_true)), "\n")
}

# C) Kořen
g <- function(x) cos(x)-x
root_b <- BisecRoot(g, 0, 1)
root_n <- NewtonRoot(g, 0.5)
cat("Root bisekce:", root_b, " | Newton:", root_n, "\n")

# D) ODE y' = y, y(0)=1 na [0,2]
fode <- function(x,y) y
h <- 0.1; xs <- seq(0,2,h); y_e <- y_rk <- numeric(length(xs))
y_e[1] <- y_rk[1] <- 1
for(i in 1:(length(xs)-1)){
  y_e[i+1]  <- EulerStep(fode, xs[i], y_e[i], h)
  y_rk[i+1] <- RK4(fode, xs[i], y_rk[i], h)
}
# rychlá vizuální kontrola
plot(xs, y_e, type="l", lwd=2, col="red", main="ODE y'=y")
lines(xs, y_rk, lwd=2, col="darkgreen")
lines(xs, exp(xs), lwd=2, lty=2)
legend("topleft", c("Euler","RK4","exp"), col=c("red","darkgreen","black"), lwd=2, lty=c(1,1,2))

# E) Interpolace ukázka
xi <- c(0,1,2); yi <- c(0,1,4)
xx <- seq(0,2,0.05); yy <- sapply(xx, function(z) Lagrange(z, xi, yi))
plot(xi, yi, pch=19, col="red", main="Lagrange"); lines(xx, yy, lwd=2)

# F) Polynomy – Horner
a <- c(1,2,3) # 1 + 2x + 3x^2
Horner(a, 2)
```


<img width="654" height="725" alt="image" src="https://github.com/user-attachments/assets/e15a4228-67f2-468d-859e-64349330e622" />


