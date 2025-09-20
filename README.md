# Numerické metody

## Derivace

### Dopředná

```{r}
ForwardDiff <- function(f,x,h=1e-6) {
	return((f(x+h)-f(x))/h)
}
```

```r

x <- seq(-5, 5, 0.1)
y <- sapply(x, function(x0) ForwardDiff(f, x0))

plot(x, y, type="l", lwd=2 ,col="red")
lines(x, 2*x, col="blue", lwd=2, lty=2)   # přesná derivace

```

Spočte derivaci spojité funkce $f(x)$ v bodě $x$ s krokem $h$.

### Centrální diference

```{r}
CentralDiff <- function(f,x,h=1e-6) {
	return((f(x+h)-f(x-h))/(2*h))
}
```

Obdobně jako dopředná, jen přesnější.

## Integrály

### Pravidlo středního bodu

```{r}
MidpointRule <- function(f,a,b,n=1e6) {
	h <- (b-a)/n
	return(h*sum(f(a+h*(1:n)-h/2)))
}
```
<img width="741" height="598" alt="image" src="https://github.com/user-attachments/assets/76a46137-0ccf-45d2-8935-31022107fe54" />


```r
f <- function(x) x^2

a <- 0
b <- 1
n <- 5                        # počet obdélníků
h <- (b-a)/n

# Midpoint hodnota
MidpointRule(f,a,b,n)

# Vykreslení
x <- seq(a, b, 0.01)
plot(x, f(x), type="l", lwd=2, col="blue", main="Midpoint Rule", ylab="f(x)")

# nakreslíme obdélníky
for(i in 1:n){
  mid <- a + (i-0.5)*h
  rect(a+(i-1)*h, 0, a+i*h, f(mid), border="red", col=rgb(1,0,0,0.3))
}
```


Spočte integrál $f(x)$ na intervalu $(a,b)$.

### Simpsonovo pravidlo

```{r}
SimpsonRule <- function(f,a,b,n=1){
  h <- (b-a)/n
  hpul <- h/2
  suma <- f(a) + f(b)
  xs <- h*(1:n) + a - hpul
  suma <- suma + 4*sum(f(xs))
  if(n > 1){
    xl <- (xs + hpul)[-n]
    suma <- suma + 2*sum(f(xl))
  }
  return(suma*h/6)
}
```

Obdobně jako pravidlo středního bodu.

### Monte Carlo

```{r}
GeomMethod <- function(f,a,b,h,n=1e6){
	x <- runif(n,a,b)
	y <- runif(n,0,h)
	return(h*(b-a)*sum(f(x) > y)/n)
}
```

Obdobně jako pravidlo středního bodu.

## Soustavy rovnic

### Gaussova eliminace

```{r}
GaussPivot <- function(A,b) {
	Ab <- cbind(A, b)
	n <- length(b)
	for(k in 1:(n-1)) {
		pivot <- which.max(abs(Ab[k:n,k])) + k - 1
		if(pivot != k){
			j <- k:(n+1)
			pom <- Ab[k,j]
			Ab[k,j] <- Ab[pivot,j]
			Ab[pivot,j] <- pom
		}
		for (i in (k+1):n) {
			j <- (k+1):(n+1)
			nasobek <- Ab[i,k] / Ab[k,k]
			Ab[i,j] <- Ab[i,j] - nasobek*Ab[k,j]
		}
	}
	x <- b
	x[n] <- Ab[n,n+1] / Ab[n,n]
	for(i in (n-1):1){
		j <- (i+1):n
		x[i] <- (Ab[i,n+1]-sum(Ab[i,j]*x[j]))/Ab[i,i]
	}
	return(x)
}
```

<img width="749" height="386" alt="image" src="https://github.com/user-attachments/assets/47938b17-4755-47c3-b988-56e285d992d0" />


```r

A <- matrix(c(2,1,1,-1), nrow=2, byrow=TRUE)
b <- c(5,1)

GaussPivot(A,b)   # vlastní funkce
solve(A,b)        # vestavěná funkce v R

```

Řeší soustavu rovnic zadanou ve tvaru $Ax = b$.

## Interpolace

<!--
## Newton
-->

### Lagrange

```{r}
Lagrange <- function(xa,x,y){
  n <- length(x)
  suma <- 0
  for(i in 1:n){
    nasobic <- 1
    for(j in 1:n){
      if(j != i) nasobic <- nasobic*(xa-x[j])/(x[i]-x[j])
    }
    suma <- suma + y[i]*nasobic
  }
  return(suma)
}
```

<img width="838" height="618" alt="image" src="https://github.com/user-attachments/assets/41c9cbb5-0663-439a-af08-b405bbb7afcd" />

```r

x <- c(0,1,2)
y <- c(0,1,4)

# interpolace v bodě 1.5
Lagrange(1.5, x, y)

# jemná síť pro hladký graf
xx <- seq(0,2,0.1)
yy <- sapply(xx, function(z) Lagrange(z, x, y))

plot(x, y, col="red", pch=19, xlim=c(0,2), ylim=c(0,4), main="Lagrangeova interpolace")
lines(xx, yy, col="blue", lwd=2)
legend("topleft", legend=c("body","interpolace"), col=c("red","blue"), pch=c(19,NA), lty=c(NA,1))


```

Spočte hodnotu interpolačního polynomu v bodě $x_a$ pro body zadané vektory $\vec{x}$, $\vec{y}$.

## Aproximace

### LSA

```{r}
LSA <- function(x,y,n){
  X <- matrix(1, nrow=n, ncol=length(x))
  for(i in 2:n) X[i, ] <- X[i-1, ]*x
  return(c(solve(X%*%t(X), X%*%y)))
}
```

<img width="773" height="1015" alt="image" src="https://github.com/user-attachments/assets/cf54da82-7459-45a8-af64-a2ecfb1add9d" />


```r
set.seed(1)
x <- seq(-3, 3, 0.5)
y <- x^2 + rnorm(length(x), sd=2)   # x^2 + šum

# fit kvadratického polynomu (n=3 → a0 + a1*x + a2*x^2)
coef <- LSA(x,y,n=3)
coef
xx <- seq(-3, 3, 0.1)
yy <- coef[1] + coef[2]*xx + coef[3]*xx^2

plot(x, y, col="red", pch=19, main="Aproximace dat polynomem")
lines(xx, yy, col="blue", lwd=2)
legend("topleft", legend=c("data","aproximace"), col=c("red","blue"), pch=c(19,NA), lty=c(NA,1))


```

Spočte koeficienty aproximačního polynomu stupně $n-1$ pro body zadané vektory $\vec{x}$, $\vec{y}$.

## Diferenciální rovnice

### Eulerova metoda

```{r}
EulerStep <- function(f,x,y,h) {
	return(y+h*f(x,y))
}
```
<img width="773" height="1178" alt="image" src="https://github.com/user-attachments/assets/9c66eb6c-a823-4c6e-9c3b-c579bbeae344" />

```r
f <- function(x,y) y

# parametry
h <- 0.1
x <- seq(0, 2, h)
y <- numeric(length(x))
y[1] <- 1   # počáteční podmínka

# Eulerův krok
for(i in 1:(length(x)-1)){
  y[i+1] <- EulerStep(f, x[i], y[i], h)
}

# graf
plot(x,y, type="l", col="red", lwd=2, main="Eulerova metoda")
lines(x, exp(x), col="blue", lwd=2, lty=2)
legend("topleft", legend=c("Euler","přesné e^x"), col=c("red","blue"), lwd=2, lty=c(1,2))



```

Řeší jeden krok obyčejné diferenciální rovnice zadané ve formě $f(x,y)$ s krokem $h$.

### RK4

```{r}
RK4 <- function(f,x,y,h){
  hhalf <- 0.5*h
  k1 <- f(x,y)
  k2 <- f(x+hhalf, y+hhalf*k1)
  k3 <- f(x+hhalf, y+hhalf*k2)
  k4 <- f(x+h, y+h*k3)
  return(y+h*(k1+2*(k2+k3)+k4)/6)
}
```

Obdobně jako Eulerova metoda, jen přesnější.

## Hledání kořenů

### Bisekce

```{r}
BisecRoot <- function(f,a,b){
	fa <- f(a)
	fb <- f(b)
	if(fa*fb<0) {
		repeat {
			c <- (a+b)/2
			if(c==a | c==b) return(c)
			fc <- f(c)
			if(fa*fc<0) {
				b <- c
				fb <- fc
			} else {
				a <- c
				fa <- fc
			}
		}
	} else {
		stop("f(a)*f(b) < 0 not satisfied")
	}
}
```
<img width="735" height="1001" alt="image" src="https://github.com/user-attachments/assets/fa6f24ec-ffd3-4616-a058-291a22404828" />

```r

f <- function(x) cos(x) - x

# zvolíme interval, kde se mění znaménko
a <- 0; b <- 1
root_bisec <- BisecRoot(f, a, b)
root_bisec


xx <- seq(0,1,0.001)
plot(xx, f(xx), type="l", lwd=2, main="f(x)=cos x - x")
abline(h=0, lty=2)
abline(v=root_bisec, col="red")

```

Nalezne kořen spojité funkce $f(x)$ na intervalu $(a,b)$. $f(a)$ a $f(b)$ se musí lišit znaménkem.

### Newtonova metoda

```{r}
NewtonRoot <- function(f,x0,tol=1e-6) {
	x <- x0
    h <- tol
	repeat{
        dx <- 2*h*f(x)/(f(x+h)-f(x-h))
		if(abs(dx) < tol) return(x)
		x <- x - dx
	}
}
```

<img width="749" height="1203" alt="image" src="https://github.com/user-attachments/assets/6af47cb3-e829-41e7-88b8-fcf96a68c627" />

```r
f <- function(x) cos(x) - x

root_newton <- NewtonRoot(f, x0 = 0.5)   # rozumný start v okolí kořene
root_newton

```

Nalezne kořen spojité funkce $f(x)$ při počátečním odhadu $x_0$.

## Polynomy

### Horner

```{r}
Horner <- function(a,x){
	n <- length(a)
	y <- a[n]
	for(i in (n-1):1) y <- y*x+a[i]
	return(y)
}
```
<img width="786" height="658" alt="image" src="https://github.com/user-attachments/assets/84c77485-69fd-4dcf-ab99-899e4bbaec04" />

```r
a <- c(1,2,3)   # 1 + 2x + 3x^2
Horner(a, 2)    # výpočet v bodě x=2
```

Efektivě počítá hodnotu polynomu v bodě $x$, jež je určen vektorem koeficientů $\vec{a}$.

### Čebyšev

```{r}
ChebyCoef <- function(n){
	a0 <- numeric(n)
	a0[1] <- 1
	if(n==1) return(a0)
	a1 <- numeric(n)
	a1[2] <- 1
	if(n==2) return(a1)
	for(i in 3:n){
		a <- 2*c(0, a1[-n])-a0
		a0 <- a1
		a1 <- a
	}
	return(a)
}
```
<img width="818" height="1174" alt="image" src="https://github.com/user-attachments/assets/69b54d7a-b9ec-4254-9803-b488bafd3c7c" />


```r
ChebyCoef(1)  # T0(x) = 1
ChebyCoef(2)  # T1(x) = x
ChebyCoef(3)  # T2(x) = 2x^2 - 1
ChebyCoef(4)  # T3(x) = 4x^3 - 3x

x <- seq(-1,1,0.01)

plot(x, Horner(ChebyCoef(2), x), type="l", col="red", lwd=2, ylim=c(-1,1), main="Čebyševovy polynomy")
lines(x, Horner(ChebyCoef(3), x), col="blue", lwd=2)
lines(x, Horner(ChebyCoef(4), x), col="green", lwd=2)
legend("topright", legend=c("T1","T2","T3"), col=c("red","blue","green"), lwd=2)
```

Vrací koeficienty Čebyševova polynomu stupně $n-1$.

#### -----

Derivace

- Dopředná – když jsem na začátku intervalu, potřebuji rychlý odhad.

- Centrální – když mám hodnoty vlevo i vpravo, přesnější.

Integrály

- Střední bod – rychlý přibližný výpočet plochy pod křivkou.

- Simpson – přesnější výpočet integrálu, lepší pro hladké funkce.

- Monte Carlo – když je integrál složitý nebo ve více rozměrech.

Soustavy rovnic

- Gaussova eliminace – klasická metoda na řešení soustavy Ax=b.

Interpolace

- Lagrange – spočítá hodnotu mezi známými body.

Aproximace

- LSA (nejmenší čtverce) – když mám data se šumem, najdu polynom, který se nejlépe hodí.

Diferenciální rovnice

- Eulerova metoda – jednoduchá, hrubý odhad řešení.

- RK4 – přesnější než Euler, používá se skoro vždy.

Hledání kořenů

- Bisekce – pomalá, ale jistá, když se funkce na intervalu mění ze + na –.

- Newtonova metoda – rychlá, když mám dobrý start.

Polynomy

- Horner – rychlé a jednoduché vyčíslení polynomu.

- Čebyšev – speciální polynomy na lepší aproximaci funkcí.
--------


Derivace = sklon.

Integrály = plocha.

Soustavy rovnic = řešení lineárních úloh.

Interpolace = hodnota mezi body.

Aproximace = fit k datům.

Diferenciální rovnice = vývoj v čase.

Kořeny = kde je nula.

Polynomy = rychlý výpočet/lepší aproximace.


---- 

zakladni prace s grafy a funkcemi
```r
f <- function(x) {
  return(x^2)
}

x = seq(-5, 5, 1)
y = f(x)

plot(x,y)


g = function(x) x^2

lines(x, g(x))
```


Lorenzův atraktor


<img width="614" height="1186" alt="image" src="https://github.com/user-attachments/assets/0b85cb9c-e82a-4f23-8dee-b9c57223db97" />


```r

# Runge-Kutta 4. řádu pro systém
RK4_step <- function(f, t, y, h) {
  k1 <- f(t, y)
  k2 <- f(t + h/2, y + h/2*k1)
  k3 <- f(t + h/2, y + h/2*k2)
  k4 <- f(t + h, y + h*k3)
  return(y + h*(k1 + 2*k2 + 2*k3 + k4)/6)
}

# Lorenzův systém
lorenz <- function(t, state, sigma=10, rho=28, beta=8/3) {
  x <- state[1]; y <- state[2]; z <- state[3]
  dx <- sigma * (y - x)
  dy <- x * (rho - z) - y
  dz <- x*y - beta*z
  return(c(dx, dy, dz))
}
# parametry
h <- 0.01
n <- 10000
t <- numeric(n)
out <- matrix(0, nrow=n, ncol=3)
out[1,] <- c(1, 1, 1)   # počáteční stav

for(i in 1:(n-1)) {
  out[i+1,] <- RK4_step(function(t,y) lorenz(t,y), t[i], out[i,], h)
  t[i+1] <- t[i] + h
}
library(rgl)  # balíček na 3D graf

par(mfrow=c(1,3))
plot(out[,1], out[,2], type="l", col="red", main="x-y")
plot(out[,1], out[,3], type="l", col="blue", main="x-z")
plot(out[,2], out[,3], type="l", col="green", main="y-z")

```
