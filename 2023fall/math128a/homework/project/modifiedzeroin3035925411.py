def IQI(f,  x0, x1, x2):
    """Given three pairs of points (x0, f0) , (x1, f1) , (x2, f2), IQI defines
    a quadratic polynomial x(f) s.t. x(0) =: x3"""
    f0 = f(x0); f1 = f(x1); f2 = f(x2)
    #print(f"f({x0})={f0} f({x1})={f1} f({x2})={f2}")
    a0 = f1*f2/((f0-f1)*(f0-f2))
    a1 = f0*f2/((f1-f0)*(f1-f2))
    a2 = f0*f1/((f2-f0)*(f2-f1))
    return a0*x0 + a1*x1 + a2*x2
    
from numpy import sign

def bisection_step(f,a,m,b) -> tuple:
    if sign(f(a)) == sign(f(m)):
        [a,b] = [m,b]
    elif sign(f(m)) == sign(f(b)):
        [a,b] = [a,m] 
    return [a,b]

def modifiedzeroin(f, interval, params) -> tuple[float, bool]:
    """Returns root, where f(root) = 0 over the given interval."""
    [a,b] = interval
    xs = []
    for i in range(25):
        #assert f(a)*f(b) < 0
        m=(a+b)/2
        x = IQI(f, a, b, m)
        xs.append(x)
        print(f"IQI([{a},{m},{b}]) -> {x}")
        if not (a < x < b):
            [a,b] = bisection_step(f,a,m,b)
        elif not any(abs(x) < 0.5 * abs(prev) for prev in xs[-4:]) and len(xs) > 4:
            [a,b] = bisection_step(f,a,m,b)

    return xs[-1], True

from functools import lru_cache
from math import exp, cos, sqrt
TOL = 1e-15


def test1():
    @lru_cache(maxsize=None)
    def f(x): return x*exp(-x)-2*x+1
    print(f.cache_info())
    interval = [0, 3]
    root = 0.671553094250269
    approx, success = modifiedzeroin(f=f,interval=interval,params=None)
    print(f.cache_info().currsize)
    assert abs(root-approx) < TOL
test1()

def test2():
    @lru_cache(maxsize=None)
    def f(x): return x*cos(x)-2*x*x+3*x-1
    print(f.cache_info())
    interval = [1, 3]
    root = 1.256623322505569
    approx, success = modifiedzeroin(f=f,interval=interval,params=None)
    print(f.cache_info().currsize)
    assert abs(root-approx) < TOL

def test3():
    @lru_cache(maxsize=None)
    def f(x): return x*x*x - 7*x*x + 14*x - 6
    print(f.cache_info())
    interval = [0, 1]
    root = 0.585786437626905
    approx, success = modifiedzeroin(f=f,interval=interval,params=None)
    print(f.cache_info().currsize)
    assert abs(root-approx) < TOL

def test4():
    @lru_cache(maxsize=None)
    def f(x): return sqrt(x) - cos(x) 
    print(f.cache_info())
    interval = [0, 1]
    root = 0.641714370872883
    approx, success = modifiedzeroin(f=f,interval=interval,params=None)
    print(f.cache_info())
    assert abs(root-approx) < TOL

def test5():
    @lru_cache(maxsize=None)
    def f(x): return 2*x*cos(2*x) - (x+1)*(x+1)
    print(f.cache_info())
    interval = [-4, -2]
    root = -2.191308011797247
    approx, success = modifiedzeroin(f=f,interval=interval,params=None)
    print(f.cache_info().currsize)
    assert abs(root-approx) < TOL

if __name__ == '__main__':
    test1()
    test2()
    test3()
    test4()
    test5()