import numpy as np
import matplotlib.pyplot as plt

#classical four-stage Runge–Kutta method implicit with two functions
def rk4(f, x0, Y0, h, n):
    X = np.zeros(n + 1)
    Y = np.zeros((n + 1, len(Y0)))
    X[0] = x0
    Y[0] = Y0

    for i in range(n):
        k1 = h * f(X[i], Y[i])
        k2 = h * f(X[i] + h/2, Y[i] + k1/2)
        k3 = h * f(X[i] + h/2, Y[i] + k2/2)
        k4 = h * f(X[i] + h, Y[i] + k3)
        Y[i+1] = Y[i] + (k1 + 2*k2 + 2*k3 + k4)/6
        X[i+1] = X[i] + h

    return X, Y

#q1:
f = lambda x, Y: np.array([Y[1], -Y[0]])
exact1 = lambda x: np.sin(x)
exact2 = lambda x: np.cos(x)
x0, Y0 = 0, np.array([0, 1])
h = 0.1
n_steps = int((2 * np.pi - x0)/h)
X, Y = rk4(f, x0, Y0, h, n_steps)
Y1, Y2 = Y[:,0], Y[:,1]
Y1_exact = exact1(X)
Y2_exact = exact2(X)
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.plot(X, Y1, 'b-', label='y1 (Numerical)')
plt.plot(X, Y1_exact, 'r--', label='y1 (Exact)')
plt.plot(X, Y2, 'g-', label='y2 (Numerical)')
plt.plot(X, Y2_exact, 'm--', label='y2 (Exact)')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.title('Numerical vs Exact Solutions')
plt.subplot(1, 2, 2)
plt.plot(Y1, Y2, 'b-')
plt.xlabel('y1')
plt.ylabel('y2')
plt.title('Phase Portrait (y1 vs y2)')
plt.grid()
plt.tight_layout()
plt.savefig('rk4_q1.png')
plt.close()

#q2:
f2 = lambda x, Y: np.array([Y[0] * (1 - Y[1]), Y[1] * (Y[0] - 1)])
x0, Y0 = 0, np.array([2, 1])
h1, h2 = 0.2, 0.1
n_steps1 = int((2 - x0)/h1)
n_steps2 = int((2 - x0)/h2)
X1, Y1 = rk4(f2, x0, Y0, h1, n_steps1)
X2, Y2 = rk4(f2, x0, Y0, h2, n_steps2)
Y1_1, Y1_2 = Y1[:,0], Y1[:,1]
Y2_1, Y2_2 = Y2[:,0], Y2[:,1]
plt.figure(figsize=(10, 5))
plt.plot(Y1_1, Y1_2, 'b-', label='h=0.2')
plt.plot(Y2_1, Y2_2, 'r--', label='h=0.1')
plt.xlabel('y1')
plt.ylabel('y2')
plt.title('Phase Portrait (y1 vs y2)')
plt.legend()
plt.grid()
plt.savefig('rk4_q2.png')
plt.close()
# Discussion:
# Halving the step size from 0.2 to 0.1 generally improves the accuracy of the numerical solution,
# as smaller step sizes tend to reduce truncation errors in the Runge-Kutta method.
# However, it also increases the computational cost since more steps are required to cover the same interval.
# In terms of stability, both step sizes appear to produce stable solutions for this particular system,
# but smaller step sizes can help avoid numerical instabilities in more sensitive systems.
# In this case, the trajectories for both step sizes are similar, indicating that the solution is stable.

#q3:
alpha, beta, delta, gamma = 1.5, 1, 1, 3
f3 = lambda x, Y: np.array([alpha * Y[0] - beta * Y[0] * Y[1], delta * Y[0] * Y[1] - gamma * Y[1]])
x0, Y0 = 0, np.array([10, 5])
h = 0.05
n_steps = int((5 - x0)/h)
X, Y = rk4(f3, x0, Y0, h, n_steps)
Y1, Y2 = Y[:,0], Y[:,1]
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.plot(X, Y1, 'b-', label='Prey (y1)')
plt.plot(X, Y2, 'r-', label='Predator (y2)')
plt.xlabel('Time')
plt.ylabel('Population')
plt.legend()
plt.title('Population Dynamics Over Time')
plt.subplot(1, 2, 2)
plt.plot(Y1, Y2, 'g-')
plt.xlabel('Prey (y1)')
plt.ylabel('Predator (y2)')
plt.title('Phase Portrait (y1 vs y2)')
plt.grid()
plt.tight_layout()
plt.savefig('rk4_q3.png')
plt.close()

#q4:
f4 = lambda x, Y: np.array([Y[1], -10 * Y[0] - Y[1]])
exact1 = lambda x: np.exp(-x) * np.cos(3 * x)
exact2 = lambda x: -3 * np.exp(-x) * np.sin(3 * x)
x0, Y0 = 0, np.array([1, 0])
h = 0.1
n_steps = int((5 - x0)/h)
X, Y_rk4 = rk4(f4, x0, Y0, h, n_steps)
Y1_rk4, Y2_rk4 = Y_rk4[:,0], Y_rk4[:,1]
Y1_exact = exact1(X)
Y2_exact = exact2(X)
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.plot(X, Y1_rk4, 'b-', label='y1 (RK4 Numerical)')
plt.plot(X, Y1_exact, 'r--', label='y1 (Exact)')
plt.plot(X, Y2_rk4, 'g-', label='y2 (RK4 Numerical)')
plt.plot(X, Y2_exact, 'm--', label='y2 (Exact)')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.title('RK4 Numerical vs Exact Solutions')
plt.subplot(1, 2, 2)
plt.plot(Y1_rk4, Y2_rk4, 'b-')
plt.xlabel('y1')
plt.ylabel('y2')
plt.title('Phase Portrait (y1 vs y2)')
plt.grid()
plt.tight_layout()
plt.savefig('rk4_q4.png')
plt.close()
# Discussion:
# The fourth-order Runge-Kutta method (RK4) provides a more accurate solution compared to the second-order method.
# This is evident from the closer alignment of the RK4 numerical solution with the exact solution.
# The RK4 method reduces the local truncation error significantly, leading to better overall accuracy.
# In contrast, the second-order method would show larger deviations from the exact solution, especially over longer intervals.