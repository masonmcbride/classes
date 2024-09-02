import math
import numpy as np
import matplotlib.pyplot as plt

def Section5_9Question2b():
    print("\\\ Chapter 5 Section 9. Question 2(b).")
    from functions import Runge_Kutta_for_Systems

    interval = [a:=0, b:=1] # t in [0,1]
    equations = [lambda t,u1,u2: u1 - u2 + 2, 
                 lambda t,u1,u2: -u1 + u2 + 4*t]
    initial_conditions = [-1,
                           0]
    N = 10 
    h = 0.1
    actual = [(t,
               (lambda t: -math.exp(2*t)/2 + t*t + 2*t - 1/2)(t),
               (lambda t: math.exp(2*t)/2 + t*t - 1/2)(t))
              for t in [a + i*h for i in range(N+1)]]


    iterations = list(Runge_Kutta_for_Systems(equations, initial_conditions, interval, N))
    error = [np.array(w) - np.array(y) for w,y in zip(iterations, actual)]
    from pprint import pprint
    pprint(error)


def Section5_9Question4b():
    print("\\\ Chapter 5 Section 9. Question 4(b).")
    from functions import Runge_Kutta_for_Systems

    interval = [a:=1, b:=3] # t in [1,3]
    equations = [lambda t,y,y1,y2: (3*t + t*y1 + t*t*y2)/4, #y
                 lambda t,y,y1,y2: (-3*t + 4*y - t*t*y2)/t, #y'
                 lambda t,y,y1,y2: (-3*t + 4*y - t*y1)/(t*t)] #y''
    initial_conditions = [4,
                          3,
                         equations[2](1,4,3,None)]
    N = 10 
    h = 0.2
    actual = [(t,
               (lambda t: 2*t*t + t + 1/(t*t))(t),
               (lambda t: 4*t + 1 - 2/(t*t*t))(t),
               (lambda t: 4 + 6/(t*t*t*t))(t)) 
               for t in [a + i*h for i in range(N+1)]]

    iterations = list(Runge_Kutta_for_Systems(equations, initial_conditions, interval, N))
    error = [np.array(w) - np.array(y) for w,y in zip(iterations, actual)]
    from pprint import pprint
    pprint(error)

def Section5_9Question5():
    print("\\\ Chapter 5 Section 9. Question 5.")
    from functions import Runge_Kutta_for_Systems
    interval = [a:=0, b:=4]
    k1,k2,k3,k4 = 3, 2e-3, 6e-4, 1/2
    equations = [lambda t,x1,x2: k1*x1 - k2*x1*x2, # x1'
                 lambda t,x1,x2: k3*x1*x2 - k4*x2] # x2'
    initial_conditions = [1000,
                          500]
    h = 0.01
    
    N = int(round((b-a)/h))
    print(f"{N=}")

    iterations = list(Runge_Kutta_for_Systems(equations, initial_conditions, interval, N))
    t = np.arange(0,4+h,h)
    prey = [x1 for _,x1,_ in iterations]
    predator = [x2 for _,_,x2 in iterations]
    # plot that shit
    fig = plt.figure()
    fig.suptitle('Predicting the Population of a Prey and Predator', fontsize=20)
    plt.plot(t, prey, label='Prey', color='b',linewidth=0.5)
    plt.plot(t, predator, label='Predator', color='y', linewidth=1.5)
    plt.xlabel('Time (s.)', fontsize=20)
    plt.ylabel('Population (num.)', fontsize=20)
    plt.legend()
    plt.show()

def Section5_9Question6():
    print("\\\ Chapter 5 Section 9. Question 6.")
    from functions import Runge_Kutta_for_Systems
    interval = [a:=0, b:=4]
    equations = [lambda t,x1,x2: x1*(4 - 3e-4*x1 - 4e-4*x2), # x1'
                 lambda t,x1,x2: x2*(2 - 2e-4*x1 - 1e-4*x2)] # x2'
    initial_conditions = [10000,
                          5000]
    h = 0.01
    
    N = int(round((b-a)/h))
    print(f"{N=}")

    iterations = list(Runge_Kutta_for_Systems(equations, initial_conditions, interval, N))
    t = np.arange(0,4+h,h)
    prey = [x1 for _,x1,_ in iterations]
    predator = [x2 for _,_,x2 in iterations]
    # plot that shit
    fig = plt.figure()
    fig.suptitle('Predicting the Population of a Prey and Predator', fontsize=20)
    plt.plot(t, prey, label='Prey', color='b',linewidth=0.5)
    plt.plot(t, predator, label='Predator', color='y', linewidth=1.5)
    plt.xlabel('Time (s.)', fontsize=20)
    plt.ylabel('Population (num.)', fontsize=20)
    plt.legend()
    plt.show()

def Section5_11Question2b():
    print("\\\ Chapter 5 Section 11. Question 2(b)")
    from functions import eulers_method
    interval = [0, 1]
    def y1(t, y): return -10*y + 10*t + 1
    def y(t): return math.exp(-5*t) + math.exp(t)
    alpha = math.e
    N = 10
    t, w = eulers_method(f=y1,alpha=alpha,interval=interval, N=N)
    print(t)
    print(w[-1])
    print(y(1))

def Section5_11Question6b():
    interval = [0, 1]
    def y1(t, y): return -10*y + 10*t + 1
    alpha = math.e
    N = 10
    pass

def Section5_11Question8b():
    from functions import Trapezoidal_Algorithm
    interval = [0, 1]
    def y1(t, y): return -10*y + 10*t + 1
    def fy(t, y): return -10
    def y(t): return math.exp(-5*t) + math.exp(t)
    alpha = math.e
    N = 10
    solutions = list(Trapezoidal_Algorithm(y1, fy, alpha, interval,N, 1e-5))
    print(solutions)
    print(y(1))

#### TEST BED ####
if __name__ == '__main__':
    Section5_9Question5

