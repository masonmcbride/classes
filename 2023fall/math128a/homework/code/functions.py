import math
from functools import reduce
from operator import mul
import numpy as np
#from sympy import *
from pprint import pprint

epsilon = 1e-4

def bisection(f, a, b): 
    if np.sign(f(a)) == np.sign(f(b)):
        raise Exception(
         "The scalars a and b do not bound a root")
        
    m = (a + b)/2
    
    if np.abs(f(m)) < epsilon:
        # stopping condition, report m as root
        return m
    elif np.sign(f(a)) == np.sign(f(m)):
        return bisection(f, m, b)
    elif np.sign(f(b)) == np.sign(f(m)):
        return bisection(f, a, m)

def newton(fun, derfun, p0, N):
    p_next = p0 - (fun(p0)/derfun(p0))
    if N == 0 or abs(p_next - p0) < epsilon:
        return p_next, N
    else:
        return newton(fun, derfun, p_next, N-1)

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

def steffensen(g, p0, N):
    p1 = g(p0)
    p2 = g(p1) 
    #p = p0 - ((p1-p0)*(p1-p0) / (p2 - 2*p1 + p0)) these lines should be the same
    p = p0 - ((p1-p0)*(p1-p0) / (p2 - p1)-(p1 - p0))
    print(f"{p=} for {N=}")
    if N == 0 or abs(p - p0) < epsilon:
        return p
    else:
        return steffensen(g, p, N - 1)

def muller(f, p0, p1, p2, N):
    h1, h2 = p1 - p0, p2 - p1
    d1 = (f(p1) - f(p0))/h1
    d2 = (f(p2) - f(p1))/h2
    d = (d2 - d1)/(h1 + h2)
    for i in range(3, N):
        b = d2 + h2*d
        D = sqrt(b*b - 4*f(p2)*d)
        E = max([b-D, b+D],key=abs)
    
        h = -2*f(p2)/E
        p = p2 + h
    
        if abs(h) < epsilon:
            return p, i
        
        p0, p1, p2 = p1, p2, p

        h1, h2 = p1-p0, p2-p1
        d1 = (f(p1) - f(p0))/h1
        d2 = (f(p2) - f(p1))/h2
        d = (d2 - d1)/(h1 + h2)

def P(xx):
    """Linear interpolation via Lagrange Polynomials"""
    x = {0: -0.5, 1: -0.25, 2: 0.25, 3: 0.5}
    f = {0: 1.93750, 1: 1.33203, 2: 0.800781, 3: 0.687500}
    n = len(x) - 1

    def L(n, k): 
        return reduce(mul, [(xx - x[i])/(x[k] - x[i]) for i in range(n+1) if i != k])

    return sum(f[k]*L(n, k) for k in range(n))

def neville(x, xs, fxs):
    """to find a n degree interpolating poynomial with n+1 of x0,...,xn, 
    perform neville's method"""
    n = len(xs)
    Q = np.zeros((n,n))
    # initalize first column w known values
    for i in range(n):
        Q[i,0] = fxs[i]
    
    for j in range(1,n):
        for i in range(j,n):
            Q[i, j] = (x - xs[i-j])*Q[i,j-1] - (x - xs[i])*Q[i-1,j-1]
            Q[i, j] /= xs[i] - xs[i-j]
    
    return Q

def NDDF(xs, fxs):
    n = len(xs)
    F = np.zeros((n,n))
    # initialize first column w known values
    for i in range(n):
        F[i,0] = fxs[i]
    
    for j in range(1,n):
        for i in range(j,n):
            F[i,j] = F[i,j-1] - F[i-1,j-1]
            F[i,j] /= xs[i] - xs[i-j]

    return F

def hermite(xs, fxs, f1xs):
    n = len(xs)
    Q = np.zeros((2*n+1,2*n+1))
    z = [0 for _ in range(2*n+1)]
    for i in range(n):
        z[2*i] = xs[i]
        z[2*i+1] = xs[i]
        Q[2*i,0] = fxs[i]
        Q[2*i+1,0] = fxs[i]
        Q[2*i+1,1] = f1xs[i]

        if i != 0:
            Q[2*i,1] = Q[2*i,0] - Q[2*i-1,0]
            Q[2*i,1] /= z[2*i] - z[2*i-1]
        
    for j in range(2,2*n+1):
        for i in range(j,2*n+1):
            Q[i,j] = Q[i,j-1] - Q[i-1,i-1]
            Q[i,j] /= z[i] - z[i-j]
    
    return Q
    
def H(Q, x, xs):
    m = len(Q)
    zs = [x for x in xs for _ in (0, 1)]
    return Q[0,0] + sum([Q[i,i]*reduce(mul, [(x - zs[j]) 
                                             for j in range(i)]) 
                                             for i in range(1,m)])

def interpolation(xs, fxs):
    n = len(xs)
    x = symbols('x')

    F = NDDF(xs=xs, fxs=fxs)
    P = F[0,0] + sum(F[i,i]*reduce(mul, (x-xs[j] for j in range(i))) for i in range(1,n))
    return P

def natural_cubic_spline(xs, fxs):
    n = len(xs) - 1
    a = fxs
    h = [0 for _ in range(n)]
    alpha = [0 for _ in range(n)]
    for i in range(n):
        h[i] = xs[i+1] - xs[i] 
    for i in range(n):
        alpha[i] = 3/h[i]*(a[i+1] - a[i]) - 3/h[i-1]*(a[i]-a[i-1])

    l = [1] + [0 for _ in range(1,n+1)] 
    zs = [0 for _ in range(n+1)]
    mu = [0 for _ in range(n)]
    for i in range(1,n):
        l[i] = 2*(xs[i+1] - xs[i-1]) - h[i-1]*mu[i-1]
        mu[i] = h[i]/l[i]
        zs[i] = (alpha[i] - h[i-1]*zs[i-1])/l[i]
    
    l[n] = 1
    zs[n] = 0
    c = [0 for _ in range(n+1)]
    b = [0 for _ in range(n)]
    d = [0 for _ in range(n)] 

    for j in reversed(range(n)):
        c[j] = zs[j] - mu[j]*c[j+1]
        b[j] = (a[j+1] - a[j])/h[j] - h[j]*(c[j+1] + 2*c[j])/3
        d[j] = (c[j+1] - c[j])/(3*h[j])

    return [(a[i], b[i], c[i], d[i]) for i in range(n)]

def S(xj, coeffs): 
    x = symbols('x')
    a, b, c, d = coeffs 
    return a + b*(x-xj) + c*(x-xj)*(x-xj) + d*(x-xj)*(x-xj)*(x-xj)

def clamped_cubic_splines(xs, fxs, f1xs): 
    n = len(xs) - 1
    a = fxs
    FPO, FPN = f1xs[0], f1xs[-1]

    h = [xs[i+1] - xs[i] for i in range(n)]
    alpha = [0 for _ in range(n+1)]
    alpha[0] = 3*(a[1]- a[0])/h[0] - 3*FPO
    alpha[n] = 3*FPN - 3*(a[n]-a[n-1])/h[n-1]
    for i in range(1,n):
        alpha[i] = 3/h[i]*(a[i+1] - a[i]) - 3/h[i-1]*(a[i] - a[i-1])

    l = [2*h[0]] + [0 for _ in range(1,n+1)] 
    zs = [alpha[0]/l[0]] + [0 for _ in range(1,n+1)]
    mu = [0.5] + [0 for _ in range(1,n+1)]
    for i in range(1,n):
        l[i] = 2*(xs[i+1] - xs[i-1]) - h[i-1]*mu[i-1]
        mu[i] = h[i]/l[i]
        zs[i] = (alpha[i] - h[i-1]*zs[i-1])/l[i]
    l[n] = h[n-1]*(2 - mu[n-1])
    zs[n] = (alpha[n] - h[n-1]*zs[n-1])/l[i]
    c = [zs[n] for _ in range(n+1)]
    b = [0 for _ in range(n)]
    d = [0 for _ in range(n)] 

    for j in reversed(range(n)):
        c[j] = zs[j] - mu[j]*c[j+1]
        b[j] = (a[j+1] - a[j])/h[j] - h[j]*(c[j+1] + 2*c[j])/3
        d[j] = (c[j+1] - c[j])/(3*h[j])

    return [(a[i], b[i], c[i], d[i]) for i in range(n)]

### SUB-SECTION 3.6

def bezier_curves(xs, xps, xns):
    n = len(xs) - 1
    assert len(xs) == len(xps) == len(xns)
    a0 = [0] * n
    a1 = [0] * n
    a2 = [0] * n
    a3 = [0] * n
    b0 = [0] * n
    b1 = [0] * n
    b2 = [0] * n
    b3 = [0] * n 
    for i in range(n):
        xi, yi = xs[i]
        xii, yii = xs[i+1]
        xpi, ypi = xps[i]
        xnii, ynii = xns[i+1]
        a0[i] = xi
        b0[i] = yi
        a1[i] = 3*(xpi - xi)
        b1[i] = 3*(ypi - yi)
        a2[i] = 3*(xi + xnii - 2*xpi)
        b2[i] = 3*(yi + ynii - 2*ypi)
        a3[i] = xii - xi + 3*xpi - 3*xnii
        b3[i] = yii - yi + 3*ypi - 3*ynii
    return [ (a0[i], a1[i], a2[i], a3[i], b0[i], b1[i], b2[i], b3[i]) 
            for i in range(n)]

def eulers_method(f, alpha, interval, N) -> list[float]:
    [a, b] = interval
    h = (b - a) / N

    t = {0:a}
    w = {0:alpha}

    for i in range(N):
        w[i+1] = w[i] + h*f(t[i],w[i])
        t[i+1] = t[i] + h 
        
    return list(t.values()),list(w.values())

def midpoint_method(f, alpha, N, interval):
    pass

def modified_eulers_method(f, alpha, N, interval):
    [a, b] = interval
    h = (b - a) / N
    # N = (b-a) * 1/h
    t = {i: (a+h*i) for i in range(N+1)}
    print(f"{t=}")
    assert t[N] == b

    w = {0:alpha}
    print(f"{f(t[0],w[0])=}")
    for i in range(0, N):
        w[i+1] = w[i] + (h/2)*(f(t[i],w[i]) + f(t[i+1],w[i] + h*f(t[i],w[i])))
    
    #return w[-1] # return last prediction
    return list(w.values()),list(t.values())

def Runge_Kutta_FehlBerg(f, alpha, interval, TOL):
    HMIN, HMAX = 0.05, 0.5
    [a, b] = interval
    t = a
    w = alpha
    h = HMAX
    FLAG = 1

    yield (t, w, h)
    while FLAG:
        K1 = h*f(t,w)
        K2 = h*f(t+h/4,w+K1/4)
        K3 = h*f(t+3/8*h,w+3/32*K1+9/32*K2)
        K4 = h*f(t+12/13*h,w+1932/2197*K1-7200/2197*K2+7296/2197*K3)
        K5 = h*f(t+h,w+439/216*K1-8*K2+3680/513*K3-845/4104*K4)
        K6 = h*f(t+h/2,w-8/27*K1+2*K2-3544/2565*K3+1859/4104*K4-11/40*K5)
        
        R = 1/h*abs(K1/360-128/4275*K3-2197/75240*K4+K5/50+2/55*K6)

        if R <= TOL:
            t = t + h # (Aproximation accepted)
            w = w + 25/216*K1 + 1408/2565*K3 + 2197/4104*K4 - K5/5
            yield (t, w, h)
        
        delta = 0.84*pow(TOL/R, 1/4)
        if delta <= 0.1:
            h = 0.1*h
        elif delta >= 4:
            h = 4*h
        else:
            h = delta*h # (calculate new h)
        
        h = min(h, HMAX)
        if t >= b:
            FLAG = 0
        elif t + h > b:
            h = b - t
        elif h < HMIN:
            FLAG = 0
            print("FAIL: Minimum h exceeded")

def Adams_Bashforth_2step(f, alphas, N, interval):
    assert len(alphas) == 2
    [a, b] = interval
    h = (b - a) / N
    t = {i: (a+h*i) for i in range(N+1)}
    assert t[N] == b

    w = dict(enumerate(alphas))
    print(f"{f(t[0],w[0])=}")
    for i in range(len(w)-1, N):
        w[i+1] = w[i] + h/2*(3*f(t[i],w[i]) - f(t[i-1],w[-1]))
    
    #return w[-1] # return last prediction
    return list(w.values()),list(t.values())

def Adams_Bashforth_3step(f, alphas, N, interval):
    assert len(alphas) == 3
    [a, b] = interval
    h = (b - a) / N
    t = {i: (a+h*i) for i in range(N+1)}
    assert t[N] == b

    w = dict(enumerate(alphas))
    print(f"{f(t[0],w[0])=}")
    for i in range(len(w)-1, N):
        w[i+1] = w[i] + h/12*(23*f(t[i],w[i]) - 16*f(t[i-1],w[i-1]) + 5*f(t[i-2],[w-2]))
    
    #return w[-1] # return last prediction
    return list(w.values()),list(t.values())

def Adams_Bashforth_4step(f, alphas, N, interval):
    assert len(alphas) == 4
    [a, b] = interval
    h = (b - a) / N
    t = {i: (a+h*i) for i in range(N+1)}
    assert t[N] == b

    w = dict(enumerate(alphas))
    print(f"{f(t[0],w[0])=}")
    for i in range(len(w)-1, N):
        w[i+1] = w[i] + h/24*(55*f(t[i],w[i]) \
                              - 59*f(t[i-1],w[i-1])\
                                 + 37*f(t[i-2],w[i-2]) \
                                    - 9*f(t[i-3],w[i-3]))
    
    #return w[-1] # return last prediction
    return list(w.values()),list(t.values())

def Adams_Bashforth_5step(f, alphas, N, interval):
    [a, b] = interval
    h = (b - a) / N
    t = {i: (a+h*i) for i in range(N+1)}
    assert t[N] == b

    w = dict(enumerate(alphas))
    print(f"{f(t[0],w[0])=}")
    for i in range(len(w)-1, N):
        w[i+1] = w[i] + h/720*(1901*f(t[i],w[i])\
                               -2274*f(t[i-1],w[i-1])\
                                +2616*f(t[i-2],w[i-2])\
                                -1274*f(t[i-3],w[i-3])\
                                +251*f(t[i-4],w[i-4]))
    
    #return w[-1] # return last prediction
    return list(w.values()),list(t.values())

def Runge_Kutta_for_Systems(equations, initial_conditions, interval, N):
    a, b = interval
    alphas: list[float] = initial_conditions
    t = a
    h = (b-a)/N
    m = len(equations)

    from collections.abc import Mapping
    w: list[float] = alphas
    f: list[Mapping[list[float], float]] = equations
    k = np.zeros((N,m))
    yield tuple([t] + w)
    #NOTE wonder if using itertools.products would break this bc it's not in "order"
    for i in range(N):
        for j in range(m): k[1,j] = h*f[j](t, *w)
        for j in range(m): k[2,j] = h*f[j](t + h/2, *[wi + k[1,j]/2 for wi in w])
        for j in range(m): k[3,j] = h*f[j](t + h/2, *[wi + k[2,j]/2 for wi in w])
        for j in range(m): k[4,j] = h*f[j](t+h, *[wi + k[3,j] for wi in w])
        for j in range(m): w[j] = w[j] + (k[1,j] + 2*k[2,j] + 2*k[3,j] + k[4,j])/6
        t = a + i*h + h
        yield tuple([t] + w)

def Trapezoidal_Algorithm(f, fy, alpha, interval, N, TOL, M=100):
    [a,b] = interval
    h = (b-a)/N
    t = a
    w = alpha

    for i in range(1, N+1):
        k1 = w + h/2*f(t,w)
        w0 = k1
        j = 1
        FLAG = 0
        while not FLAG:
            w = w0 - (w0 - h/2*f(t+h, w0) - k1)/(1 - h/2*fy(t+h, w0))
            if abs(w - w0) < TOL:
                FLAG = 1
            else:
                j = j + 1
                w0 = w
                if j > M:
                    print("max number of iterations exceeded")
                    return None
        t = a + i*h
        yield (t, w)
    
if __name__ == '__main__':
    pass