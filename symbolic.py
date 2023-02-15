import sympy as sp
import numpy as np
from loguru import logger

logger.warning("Obliczanie symboliczne równania elementu...")

u_1, u_2, u_3, u_4, u_5, u_6, u_7, u_8 = sp.symbols('u_1 u_2 u_3 u_4 u_5 u_6 u_7 u_8')
a, s, t, b = sp.symbols('a s t b')
k_x, k_y, p, c, q, k1, h, Tot = sp.symbols('k_x k_y p c q k1 h Tot')
alpha, beta = sp.symbols('alpha beta')

u = sp.Matrix([[u_1], [u_2], [u_3], [u_4], [u_5], [u_6], [u_7], [u_8]])

r1 = np.ones([8, 1])
r2 = np.transpose(np.matrix([-a, 0, a, a, a, 0, -a, -a]))
r3 = np.transpose(np.matrix([-b, -b, -b, 0, b, b, b, 0]))
r4 = np.transpose(np.matrix(np.array(np.transpose(r2)) * np.array(np.transpose(r3))))
r5 = np.transpose(np.matrix(np.array(np.transpose(r2)) * np.array(np.transpose(r2))))
r6 = np.transpose(np.matrix(np.array(np.transpose(r3)) * np.array(np.transpose(r3))))
r7 = np.transpose(np.matrix(np.array(np.transpose(r4)) * np.array(np.transpose(r2))))
r8 = np.transpose(np.matrix(np.array(np.transpose(r4)) * np.array(np.transpose(r3))))

A = np.concatenate((r1, r2, r3, r4, r5, r6, r7, r8), axis=1)
A1 = np.matrix([[1, s, t, s * t, s ** 2, t ** 2, t * (s ** 2), s * (t ** 2)]])

A = sp.Matrix(A)
A = A.inv()
A = np.matrix(A)

N = A1.dot(A)
B = sp.zeros(2, 8)

for k in range(0, 8):
    B[0, k] = sp.diff(N[0, k], s)
    B[1, k] = sp.diff(N[0, k], t)

C = sp.Matrix([[k_x, 0], [0, k_y]])
C = np.matrix(C)
B = np.matrix(B)
Kk = B.transpose().dot(C).dot(B)
for k in range(0, 8):
    for i in range(0, 8):
        Kk[k, i] = sp.integrate(Kk[k, i], (s, -a, a))
        Kk[k, i] = sp.integrate(Kk[k, i], (t, -b, b))
Kp = N.transpose().dot(N)

for k in range(0, 8):
    for i in range(0, 8):
        Kp[k, i] = sp.integrate(Kp[k, i], (s, -a, a))
        Kp[k, i] = sp.integrate(Kp[k, i], (t, -b, b))

Kp = p * Kp

Nc1_2 = sp.zeros(1, 8)
Nc2_3 = sp.zeros(1, 8)
Nc3_4 = sp.zeros(1, 8)
Nc4_1 = sp.zeros(1, 8)

N = sp.Matrix(N)

for k in range(0, 8):
    Nc1_2[0, k] = N[0, k].subs(s, c).subs(t, -b)
    Nc2_3[0, k] = N[0, k].subs(s, a).subs(t, c)
    Nc3_4[0, k] = N[0, k].subs(s, -c).subs(t, b)
    Nc4_1[0, k] = N[0, k].subs(s, -a).subs(t, -c)

rq = sp.zeros(1, 8)

for k in range(0, 8):
    rq[0, k] = sp.integrate(N[0, k], (s, -a, a))
    rq[0, k] = sp.integrate(rq[0, k], (t, -b, b))
rq = q * rq

Nc1_2 = np.matrix(Nc1_2)
Nc2_3 = np.matrix(Nc2_3)
Nc3_4 = np.matrix(Nc3_4)
Nc4_1 = np.matrix(Nc4_1)

k_alpha_1_2 = Nc1_2.transpose().dot(Nc1_2)
k_alpha_2_3 = Nc2_3.transpose().dot(Nc2_3)
k_alpha_3_4 = Nc3_4.transpose().dot(Nc3_4)
k_alpha_4_1 = Nc4_1.transpose().dot(Nc4_1)

Nc1_2 = sp.Matrix(Nc1_2)
Nc2_3 = sp.Matrix(Nc2_3)
Nc3_4 = sp.Matrix(Nc3_4)
Nc4_1 = sp.Matrix(Nc4_1)

for k in range(0, 8):
    for i in range(0, 8):
        k_alpha_1_2[k, i] = sp.integrate(k_alpha_1_2[k, i], (c, -a, a))
        k_alpha_2_3[k, i] = sp.integrate(k_alpha_2_3[k, i], (c, -b, b))
        k_alpha_3_4[k, i] = sp.integrate(k_alpha_3_4[k, i], (c, -a, a))
        k_alpha_4_1[k, i] = sp.integrate(k_alpha_4_1[k, i], (c, -b, b))

k_alpha_1_2 = -alpha * k_alpha_1_2
k_alpha_2_3 = -alpha * k_alpha_2_3
k_alpha_3_4 = -alpha * k_alpha_3_4
k_alpha_4_1 = -alpha * k_alpha_4_1

r_beta_1_2 = sp.zeros(1, 8)
r_beta_2_3 = sp.zeros(1, 8)
r_beta_3_4 = sp.zeros(1, 8)
r_beta_4_1 = sp.zeros(1, 8)

for k in range(0, 8):
    r_beta_1_2[0, k] = sp.integrate(Nc1_2[0, k], (c, -a, a))
    r_beta_2_3[0, k] = sp.integrate(Nc2_3[0, k], (c, -b, b))
    r_beta_3_4[0, k] = sp.integrate(Nc3_4[0, k], (c, -a, a))
    r_beta_4_1[0, k] = sp.integrate(Nc4_1[0, k], (c, -b, b))

r_beta_1_2 = beta * r_beta_1_2
r_beta_2_3 = beta * r_beta_2_3
r_beta_3_4 = beta * r_beta_3_4
r_beta_4_1 = beta * r_beta_4_1

Kk = sp.Matrix(Kk)
Kk = Kk.subs(k_x, k1).subs(k_y, k1)


k_alpha_1_2 = sp.Matrix(k_alpha_1_2)
k_alpha_1_2 = k_alpha_1_2.subs(alpha, -h)
k_alpha_2_3 = sp.Matrix(k_alpha_2_3)
k_alpha_2_3 = k_alpha_2_3.subs(alpha, -h)
k_alpha_3_4 = sp.Matrix(k_alpha_3_4)
k_alpha_3_4 = k_alpha_3_4.subs(alpha, -h)
k_alpha_4_1 = sp.Matrix(k_alpha_4_1)
k_alpha_4_1 = k_alpha_4_1.subs(alpha, -h)

r_beta_1_2 = r_beta_1_2.subs(beta, Tot*h)
r_beta_2_3 = r_beta_2_3.subs(beta, Tot*h)
r_beta_3_4 = r_beta_3_4.subs(beta, Tot*h)
r_beta_4_1 = r_beta_4_1.subs(beta, Tot*h)

logger.info("Obliczenie równania elementu zakończone")
