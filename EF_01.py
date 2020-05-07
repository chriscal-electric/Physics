import numpy as np
import sympy as sp

def v2m(v):
    m = sp.sqrt(v[0]**2 + v[1]**2)
    return m

def dist(v1, v2):
    dv = np.array([0.0, 0.0])
    dv[0] = v2[0] - v1[0]
    dv[1] = v2[1] - v1[1]
    return dv

def angle(p1, p2):
    ang = sp.atan2(p2[1] - p1[1], p2[0] - p1[0])
    return ang

def m2v(m, theta):
    v = np.array([0.0, 0.0])
    v[0] = m * sp.cos(theta)
    v[1] = m * sp.sin(theta)
    return v


# Three point loads, q1, q2, and qr are situated in a cartesian plane at points
# A, B and C, respectively. Calculate the total force that a q4 load receivese when
# situated at P

q1, q2, q3, q = -3e-6, 4e-6, 2e-6, 1e-6
A = [0,0]
B = [4,0]
C = [1,1]
P = [0,4]
k = 8.99e9


Coord = np.array([A, B, C])
Q = np.array([q1, q2, q3])
d = np.zeros(3)
theta = np.zeros(3)

for i in range(0, 3):
    d[i] = v2m(dist(Coord[i], P))
    theta[i] = angle(Coord[i], P)

Fm = np.zeros(3)

for i in range(0,3):
    Fm[i] = k * q * Q[i] / (d[i]**2)

Fv = np.zeros((3,2))

for i in range(0, 3):
    Fv[i] = m2v(Fm[i], theta[i])

Ftv = np.zeros(2)

for i in range(0,3):
    Ftv[0] = Ftv[0] + Fv[i,0]
    Ftv[1] = Ftv[1] + Fv[i,1]

print(Ftv)

Ftm = v2m(Ftv)
print(Ftm)

