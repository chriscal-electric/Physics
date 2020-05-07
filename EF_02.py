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

def NewtonRaphson(f, v, x):
    vac = sp.symbols('vac')
    nofe = len(f)
    it_max = 100
    X = np.zeros((it_max, nofe, 1))
    F = np.full((nofe, 1), vac)
    var = np.full((nofe), vac)

    for i in range(0,nofe):
        X[0,i,0] = x[i]
        F[i,0] = f[i]
        var[i] = v[i]

    DF = np.full((nofe,nofe), vac)
    DFX = np.zeros((it_max,nofe,nofe))
    DFX_inv = np.zeros((it_max,nofe,nofe))
    FX = np.zeros((it_max,nofe,1))
    H = np.zeros((it_max,nofe,1))


    tol = 0.00000001

    stop = False

    for i in range(0,nofe):
        for j in range(0,nofe):
            DF[i,j] = sp.diff(F[i,0], var[j])

    for it in range(0,it_max):

        if (stop == False):


            subby = '{'

            for i in range(0, nofe):
                subby = subby + 'var[' + str(i) + ']: X[it,' + str(i) + ',0]'
                if (i != nofe-1):
                    subby = subby + ',  '
                if (i == nofe-1):
                    subby = subby + '} '

            print(eval(subby))

            for i in range(0,nofe):
                for j in range(0,nofe):
                    DFX[it,i,j] = DF[i,j].subs(eval(subby)).evalf()

            for i in range(0,nofe):
                FX[it,i,0] = F[i,0].subs(eval(subby)).evalf()


            DFX_inv[it] = np.linalg.inv(DFX[it])


            H[it] = np.dot(DFX_inv[it], -FX[it])
            if(it+1 < it_max):

                for i in range(0,nofe):
                    X[it + 1, i, 0] = X[it, i, 0] + H[it, i, 0]

            else:
                print('NO CONVERGE')
                print('Prueba otros valores iniciales')
                stop = True

            if(abs(H[it,0,0]) <= tol):
                print('Hemos decidido parar en la iteración número: ' + str(it) + ' porque La H ya es menor que la tolerancia')
                return X[it]
                stop = True

k = 8.99e9
q1 = 2e-6
q2 = -3e-6

r1, r2 = sp.symbols('r1 r2')

E1 = k * abs(q1) / r1 ** 2
E2 = k * abs(q2) / r2 ** 2

V1 = k * q1 / r1
V2 = k * q2 / r2

print(E1)
print(E2)

f1 = E1 - E2
f2 = r1 - r2 -1

a, b = NewtonRaphson([f1, f2], [r1, r2], [0.5,-0.5])

f3 = V1 - V2

c, d = NewtonRaphson([f3, f2], [r1, r2], [0.5, -0.5])