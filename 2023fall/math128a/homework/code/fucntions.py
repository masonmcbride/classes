from math import log, exp

epsilon = 2e-4

def newton(fun, derfun, p0, n):
    p_next = p0 - (fun(p0)/derfun(0))
    if n == 0 or abs(p_next - p0) < epsilon:
        return p_next
    else:
        return newton(fun, derfun, p_next, n-1)

def secant(f, p0, p1):
    for _ in range(100):
        if f(p0) == f(p1):
            return p1, _
        
        p0, p1 = p1, p1 - f(p1)*(p1 - p0)/(f(p1) - f(p0))
        if abs(f(p1)) < epsilon:
            return p1, _
    

def fpi(g, p0) -> float:
    p = p0
    for _ in range(100):
        p_next = g(p)
        if abs(p_next - p) < epsilon:
            return p_next, _
        p = p_next
    return p, 100

def p(n): return 1./(n*n)
    
if __name__ == '__main__':
    def f(x): return exp(6*x) + 3*(log(2)*log(2))*exp(2*x) - (log(8))*exp(4*x) - (log(2)*log(2)*log(2))
    def fprime(x): return 6*exp(6*x) + 6*(log(2)*log(2))*exp(2*x) - 4*log(8)*exp(4*x)
    def g(x): return x - f(x)/fprime(x)
    
    p = newton(fun=f, derfun=fprime, p0=0, n=1000)
    print(f"newtwon {p=}")
    p, i = fpi(g=g, p0=0)
    print(f"fpi {p=}, {i=}")

