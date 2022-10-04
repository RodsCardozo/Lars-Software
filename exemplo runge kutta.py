import numpy as np
from math import sqrt
import matplotlib.pyplot as plt
import numba
import time
start_time = time.clock()

a=1
b=1
c=1
d=1

# Equations:
@numba.jit()

#du/dt=V(u,t)
def V(u,t):
    x, y, vx, vy = u
    return np.array([vy,vx,a*x+b*y,c*x+d*y])

def rk4(f, u0, t0, tf , n):
    t = np.linspace(t0, tf, n+1)
    u = np.array((n+1)*[u0])
    h = t[1]-t[0]
    for i in range(n):
        k1 = h * f(u[i], t[i])
        k2 = h * f(u[i] + 0.5 * k1, t[i] + 0.5*h)
        k3 = h * f(u[i] + 0.5 * k2, t[i] + 0.5*h)
        k4 = h * f(u[i] + k3, t[i] + h)
        u[i+1] = u[i] + (k1 + 2*(k2 + k3 ) + k4) / 6
    return u, t

u, t  = rk4(V, np.array([1., 0., 1. , 0.]) , 0. , 10. , 100000)
x,y, vx,vy  = u.T
# plt.plot(t, x, t,y)
plt.semilogy(t, x, t,y)
plt.grid('on')
print("Execution time:",time.clock() - start_time, "seconds")
plt.show()