function [root, info] = modifiedzeroin(@func,Int,params)
    a = Int.a
    b = Int.b
    
        m=(a+b)/2
        x = IQI(f, a, b, m)
        xs.append(x)
        if a < x:
            [a,b] = bisection_step(f,a,m,b)
        elif not any(abs(x) < 0.5 * abs(prev) for prev in xs[-4:]) and len(xs) > 4:
            [a,b] = bisection_step(f,a,m,b)

function [a,b] = bisection_step(@func,a,m,b):
    if sign(f(a)) == sign(f(m)):
        [a,b] = [m,b]
    elif sign(f(m)) == sign(f(b)):
        [a,b] = [a,m] 
    return [a,b]

def IQI(f,  x0, x1, x2):
    """Given three pairs of points (x0, f0) , (x1, f1) , (x2, f2), IQI defines
    a quadratic polynomial x(f) s.t. x(0) =: x3"""
    f0 = f(x0); f1 = f(x1); f2 = f(x2)
    #print(f"f({x0})={f0} f({x1})={f1} f({x2})={f2}")
    a0 = f1*f2/((f0-f1)*(f0-f2))
    a1 = f0*f2/((f1-f0)*(f1-f2))
    a2 = f0*f1/((f2-f0)*(f2-f1))
    return a0*x0 + a1*x1 + a2*x2