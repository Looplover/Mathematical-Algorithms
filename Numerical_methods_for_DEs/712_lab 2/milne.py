import numpy as np
import matplotlib.pyplot as plt

def milne_pc(f, x0, y0, h, n_steps):
    X = np.zeros(n_steps + 1)
    Y = np.zeros(n_steps + 1)
    X[0] = x0
    Y[0] = y0

    # Start-up steps using RK4
    for i in range(3):
        k1 = h * f(X[i], Y[i])
        k2 = h * f(X[i] + h/2, Y[i] + k1/2)
        k3 = h * f(X[i] + h/2, Y[i] + k2/2)
        k4 = h * f(X[i] + h, Y[i] + k3)
        Y[i+1] = Y[i] + (k1 + 2*k2 + 2*k3 + k4)/6
        X[i+1] = X[i] + h

    # Milne's Predictor-Corrector loop
    for i in range(3, n_steps):
        fn   = f(X[i],   Y[i])
        fn1  = f(X[i-1], Y[i-1])
        fn2  = f(X[i-2], Y[i-2])
        
        # Predictor
        Yp = Y[i-3] + (4*h/3) * (2*fn - fn1 + 2*fn2)
        
        fp = f(X[i] + h, Yp)
        
        # Corrector
        Yc = Y[i-1] + (h/3) * (fp + 4*fn + fn1)
        
        Y[i+1] = Yc
        X[i+1] = X[i] + h

    return X, Y

#q1:
f = lambda x, y: y - x**2 + 1
exact = lambda x: (x+1)**2 - 0.5*np.exp(x)
x0, y0 = 0, 0.5
h = 0.1
n_steps = int((1 - x0)/h)
X, Y = milne_pc(f, x0, y0, h, n_steps)
Y_exact1 = exact(X)

print("Question 1 :")
print("x      y (Milne)  y (Exact)   Error")
for i in range(len(X)):
    print(f"{X[i]:.1f}   {Y[i]:.6f}   {Y_exact1[i]:.6f}   {np.abs(Y[i] - Y_exact1[i]):.2e}")

# Plotting solution for q1
plt.figure()
plt.plot(X, Y, 'b-', label='Milne (Numerical)')
plt.plot(X, Y_exact1, 'r--', label='Exact')
plt.xlabel('x')
plt.ylabel('y')
plt.title("Milne's Method Solution (q1)")
plt.legend()
plt.grid(True)
plot_filename_q1_solution = 'milne_method_q1_solution.png'
plt.savefig(plot_filename_q1_solution)
plt.close()
print(f"Solution plot saved to {plot_filename_q1_solution}")

# Stability Analysis for q1
epsilons = np.logspace(-12, -2, 10)
max_relative_perturbations = []
for eps in epsilons:
    y0_perturbed_1 = y0 + eps
    X1_pert, Y1_pert = milne_pc(f, x0, y0_perturbed_1, h, n_steps)
    # Find the maximum perturbation over the whole interval and divide by epsilon
    max_relative_perturbations.append(np.max(np.abs(Y[1:] - Y1_pert[1:])) / eps)

plt.figure()
plt.loglog(epsilons, max_relative_perturbations, 'm-o', label='Max Perturbation / Epsilon')
plt.axhline(y=1, color='k', linestyle='--', label='Ratio = 1')
plt.xlabel('Initial Perturbation (epsilon)')
plt.ylabel('Max Perturbation / Initial Perturbation')
plt.title("Stability Analysis (q1)")
plt.legend()
plt.grid(True)
plot_filename_q1_stability = 'milne_method_q1_stability.png'
plt.savefig(plot_filename_q1_stability)
plt.close()
print(f"Stability plot saved to {plot_filename_q1_stability}\n")


# q2:
f2 = lambda x, y: np.cos(x) - y
exact2 = lambda x: 0.5 * (np.sin(x) + np.cos(x)) + 0.5 * np.exp(-x)
x0_2, y0_2 = 0, 1
h_2 = 0.2
n_steps_2 = int((2 - x0_2)/h_2)
X2, Y2 = milne_pc(f2, x0_2, y0_2, h_2, n_steps_2)

# Re-run loop to get predictor values for plotting
Y2_pred = np.zeros_like(Y2)
Y2_pred[:4] = Y2[:4]
Y2_corr = np.copy(Y2) 

for i in range(3, n_steps_2):
    fn   = f2(X2[i],   Y2_corr[i])
    fn1  = f2(X2[i-1], Y2_corr[i-1])
    fn2  = f2(X2[i-2], Y2_corr[i-2])
    
    Y2_pred[i+1] = Y2_corr[i-3] + (4*h_2/3) * (2*fn - fn1 + 2*fn2) # Predictor
    fp = f2(X2[i] + h_2, Y2_pred[i+1])
    Y2_corr[i+1] = Y2_corr[i-1] + (h_2/3) * (fp + 4*fn + fn1) # Corrector

Y2_exact = exact2(X2)

plt.figure()
plt.plot(X2, Y2_pred, 'r--', label='Predictor')
plt.plot(X2, Y2_corr, 'b-', label='Corrector')
plt.plot(X2, Y2_exact, 'g:', label='Exact')
plt.legend()
plt.xlabel('x')
plt.ylabel('y')
plt.title("Milne's Predictor-Corrector Method (q2)")

plot_filename_q2 = 'milne_method_q2.png'
plt.savefig(plot_filename_q2)
plt.close()

print(f"Question 2 :")
print(f"Plot saved to {plot_filename_q2}")

# Stability Analysis for q2
epsilons = np.logspace(-12, -2, 10)
max_relative_perturbations_2 = []
for eps in epsilons:
    y0_perturbed_2 = y0_2 + eps
    X2_pert, Y2_pert_raw = milne_pc(f2, x0_2, y0_perturbed_2, h_2, n_steps_2)
    Y2_pert_corr = np.copy(Y2_pert_raw)
    for i in range(3, n_steps_2):
        fn   = f2(X2_pert[i],   Y2_pert_corr[i])
        fn1  = f2(X2_pert[i-1], Y2_pert_corr[i-1])
        fn2  = f2(X2_pert[i-2], Y2_pert_corr[i-2])
        Yp = Y2_pert_corr[i-3] + (4*h_2/3) * (2*fn - fn1 + 2*fn2)
        fp = f2(X2_pert[i] + h_2, Yp)
        Y2_pert_corr[i+1] = Y2_pert_corr[i-1] + (h_2/3) * (fp + 4*fn + fn1)
    
    # Find the maximum perturbation over the whole interval and divide by epsilon
    max_relative_perturbations_2.append(np.max(np.abs(Y2_corr[1:] - Y2_pert_corr[1:])) / eps)

plt.figure()
plt.loglog(epsilons, max_relative_perturbations_2, 'm-o', label='Max Perturbation / Epsilon')
plt.axhline(y=1, color='k', linestyle='--', label='Ratio = 1')
plt.xlabel('Initial Perturbation (epsilon)')
plt.ylabel('Max Perturbation / Initial Perturbation')
plt.title("Stability Analysis (q2)")
plt.legend()
plt.grid(True)
plot_filename_q2_stability = 'milne_method_q2_stability.png'
plt.savefig(plot_filename_q2_stability)
plt.close()
print(f"Stability plot saved to {plot_filename_q2_stability}\n")


# q3: 
f3 = lambda x, y: -15*y + 15*np.sin(x)
x0_3, y0_3 = 0, 0
exact3 = lambda x: (225/226)*(np.sin(x) - (1/15)*np.cos(x) + (1/15)*np.exp(-15*x))
h_3 = 0.05
n_steps_3 = int((1 - x0_3)/h_3)
X3, Y3 = milne_pc(f3, x0_3, y0_3, h_3, n_steps_3)

# Re-run loop to get predictor values for plotting
Y3_pred = np.zeros_like(Y3)
Y3_pred[:4] = Y3[:4] 
Y3_corr = np.copy(Y3)

for i in range(3, n_steps_3):
    fn   = f3(X3[i],   Y3_corr[i])
    fn1  = f3(X3[i-1], Y3_corr[i-1])
    fn2  = f3(X3[i-2], Y3_corr[i-2])
    
    Y3_pred[i+1] = Y3_corr[i-3] + (4*h_3/3) * (2*fn - fn1 + 2*fn2) # Predictor
    fp = f3(X3[i] + h_3, Y3_pred[i+1])
    Y3_corr[i+1] = Y3_corr[i-1] + (h_3/3) * (fp + 4*fn + fn1) # Corrector
    
Y3_exact = exact3(X3)
Error3 = np.abs(Y3_corr - Y3_exact)

plt.figure()
plt.plot(X3, Y3_pred, 'r--', label='Predictor')
plt.plot(X3, Y3_corr, 'b-', label='Corrector')
plt.plot(X3, Y3_exact, 'g:', label='Exact')
plt.legend()
plt.xlabel('x')
plt.ylabel('y')
plt.title("Milne's Method on 3rd equation (q3)")

plot_filename_q3 = 'milne_method_q3.png'
plt.savefig(plot_filename_q3)
plt.close()

print(f"Question 3 :")
print(f"Plot saved to {plot_filename_q3}")

print("\nx      y (Milne)  y (Exact)   Error")
for i in range(len(X3)):
    print(f"{X3[i]:.2f}   {Y3_corr[i]:.6f}   {Y3_exact[i]:.6f}   {Error3[i]:.2e}")

# Stability Analysis for q3
epsilons = np.logspace(-12, -2, 10)
max_relative_perturbations_3 = []
for eps in epsilons:
    y0_perturbed_3 = y0_3 + eps
    X3_pert, Y3_pert_raw = milne_pc(f3, x0_3, y0_perturbed_3, h_3, n_steps_3)
    
    Y3_pert_corr = np.copy(Y3_pert_raw)
    for i in range(3, n_steps_3):
        fn   = f3(X3_pert[i],   Y3_pert_corr[i])
        fn1  = f3(X3_pert[i-1], Y3_pert_corr[i-1])
        fn2  = f3(X3_pert[i-2], Y3_pert_corr[i-2])
        
        Yp = Y3_pert_corr[i-3] + (4*h_3/3) * (2*fn - fn1 + 2*fn2) # Predictor
        fp = f3(X3_pert[i] + h_3, Yp)
        Y3_pert_corr[i+1] = Y3_pert_corr[i-1] + (h_3/3) * (fp + 4*fn + fn1) # Corrector
        
    # Find the maximum perturbation over the whole interval and divide by epsilon
    max_relative_perturbations_3.append(np.max(np.abs(Y3_corr[1:] - Y3_pert_corr[1:])) / eps)

plt.figure()
plt.loglog(epsilons, max_relative_perturbations_3, 'm-o', label='Max Perturbation / Epsilon')
plt.axhline(y=1, color='k', linestyle='--', label='Ratio = 1')
plt.xlabel('Initial Perturbation (epsilon)')
plt.ylabel('Max Perturbation / Initial Perturbation')
plt.title("Stability Analysis (q3)")
plt.legend()
plt.grid(True)
plot_filename_q3_stability = 'milne_method_q3_stability.png'
plt.savefig(plot_filename_q3_stability)
plt.close()
print(f"Stability plot saved to {plot_filename_q3_stability}")


# The instability observed is due to the stiffness of the equation.
# To improve stability, one could use a smaller step size or switch to a more stable method like implicit methods.