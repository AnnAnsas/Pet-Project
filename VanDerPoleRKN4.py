import matplotlib.pyplot as plt
import numpy as np
import json


def rk4N(r, t, h, N):
    rN = r
    hN = h / N
    for i in range (0, N):
        rN = rk4(rN, t, hN)
        t += hN
    return rN

def rk4(r, t, h):
    h2 = h / 2
    k1 = h * f(r, t)
    k2 = h * f(r + 0.5 * k1, t + h2)
    k3 = h * f(r + 0.5 * k2, t + h2)
    k4 = h * f(r + k3, t + h)
    return r + (k1 + 2 * k2 + 2 * k3 + k4) / 6

def Check (r, t, h, N, eps):
    rk = rk4 (r, t, h)
    rkN = rk4N (r, t, h, N)
    return bool (abs (rk[0] - rkN[0]) > eps or abs (rk[1] - rkN[1]) > eps)

def CheckF (r, t, h, N, eps):
    rk = rk4 (r, t,2 * h)
    rkN = rk4N (r, t, 2 * h, N)
    return bool (abs (rk[0] - rkN[0]) < eps and abs (rk[1] - rkN[1]) < eps)

def f(r, t):
    x, y = r[0], r[1]
    return np.array([fx(x, y, t), gy(x, y, t)], float)

def fx (x, y, t):
    return y

def gy (x, y, t):
    return -omega*omega*x + delta*( 1-alpha*x*x)*y


fl = open('conf0.json')
p = json.load(fl)
fl.close()

startx = 1e-2
starty = p['y0']
t = p['t0']
tmax = p['t1']
h = tmax / 10
N = 2
delta = p['delta']
alpha = p['alpha']
omega = p['omega']
eps1 = p['eps1']


xpoints, ypoints, dxpoints, dypoints, tpoints, hpoints = [], [], [], [], [], []
r = np.array([startx, starty], float)
hmax = 0
hmin = tmax
g = 0


while t <= tmax:
    xpoints.append(r[0])
    ypoints.append(r[1])
    tpoints.append(t)
    k = 0
    while (Check(r, t, h, N, eps1)):
        h = h / N
        k += 1
    if k == 0 and CheckF(r, t, h, N, eps1):
            h = 2 * h
    r = rk4 (r, t, h)
    hpoints.append(h)
    t += h
    if h > hmax:
        hmax = h
    if h < hmin:
        hmin = h
    g += 1
    
print ("Max step: ", hmax)
print ("Min step: ", hmin)
print ("Points: ", g)

l1 = ( delta + (delta ** 2 - 4 * omega ** 2) ** (1/2) ) / 2
l2 = ( delta - (delta ** 2 - 4 * omega ** 2) ** (1/2) ) / 2
print ("x_start = ", startx,",", "y_start = ", starty)
print ("l1 = ", l1,",", "l2 = ", l2)
if type(l1)==complex and l1.imag != 0:
    if l1.real > 0:
        print ("Неустойчивый фокус")
    if l1.real < 0:
        print ("Устойчивый фокус")
elif type(l1)==complex and l1.real == 0:
    print ("Центр")
elif l1 < 0 and l2 < 0:
    print("Устойчивый узел")
elif l1*l2 < 0:
    print("Седло")
elif l1 > 0 and l2 > 0:
    print ("Неустойчивый узел")

plt.figure()
fig, ax = plt.subplots()
ax.plot (tpoints, xpoints, label = 'x')
ax.plot (tpoints, ypoints, label = 'V')
ax.set_title('VDP')
ax.set_xlabel ('Time')
ax.set_ylabel ('Function')
ax.legend(loc = 2)
plt.grid(True)
plt.show()


plt.figure()
fig, ax = plt.subplots()
ax.plot (xpoints, ypoints)#, label = 'phase')
ax.set_title('V(x)')
ax.set_xlabel ('x')
ax.set_ylabel ('V')
#ax.legend(loc = 2)
plt.grid(True)
plt.savefig('xV.png')
plt.show()





