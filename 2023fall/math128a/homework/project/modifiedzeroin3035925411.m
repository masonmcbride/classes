function [root, info] = modifiedzeroin(@func,Int,params)
    a = Int.a
    b = Int.b
    xs = []
    for i in range(25):
        m=(a+b)/2
        x = IQI(f, a, b, m)
        if not (a < x < b):
            [a,b] = bisection_step(f,a,m,b)
        elif not any(abs(x) < 0.5 * abs(prev) for prev in xs[-4:]) and len(xs) > 4:
            [a,b] = bisection_step(f,a,m,b)

function [a,b] = bisection_step(@func,a,m,b):
    if sign(f(a)) == sign(f(m)):
        [a,b] = [m,b]
    elif sign(f(m)) == sign(f(b)):
        [a,b] = [a,m] 
    return [a,b]