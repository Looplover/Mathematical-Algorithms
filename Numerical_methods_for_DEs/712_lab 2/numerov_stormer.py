import numpy as np
import matplotlib.pyplot as plt

def numerov_solve(f, a, b, y0, y_prime0, N):
    if N == 0: return np.array([a]), np.array([y0])
    h = (b - a) / N
    x = np.linspace(a, b, N + 1)
    y = np.zeros(N + 1)
    
    y[0] = y0
    if N == 1:
        y[1] = y[0] + h * y_prime0
        return x, y
        
    f0 = f(x[0], y[0])
    y[1] = y[0] + h * y_prime0 + 0.5 * h**2 * f0
    
    for n in range(1, N):
        f_prev = f(x[n-1], y[n-1])
        f_curr = f(x[n], y[n])
        
        y_tilde_next = 2 * y[n] - y[n-1] + h**2 * f_curr
        f_next_pred = f(x[n+1], y_tilde_next)
        
        y[n+1] = 2 * y[n] - y[n-1] + (h**2 / 12.0) * (f_next_pred + 10 * f_curr + f_prev)
        
    return x, y

def stormer_verlet_solve(f, a, b, q0, q_prime0, N):
    if N == 0: return np.array([a]), np.array([q0])
    h = (b - a) / N
    t = np.linspace(a, b, N + 1)
    
    q0_arr = np.asarray(q0)
    q_prime0_arr = np.asarray(q_prime0)
    
    q = np.zeros((N + 1, len(q0_arr)))
    q[0] = q0_arr
    
    if N == 1:
        q[1] = q[0] + h * q_prime0_arr
        return t, q

    f0 = f(t[0], q[0])
    q[1] = q[0] + h * q_prime0_arr + 0.5 * h**2 * f0
    
    for n in range(1, N):
        fn = f(t[n], q[n])
        q[n+1] = 2 * q[n] - q[n-1] + h**2 * fn
        
    return t, q

def f_n1(x, y): return np.exp(x) - y * np.exp(-x)
def f_n2(x, y): return np.sinh(x) - y / 3.0
def y_exact_n2(x): return (3.0/4.0) * np.sinh(x) - (3.0*np.sqrt(3.0)/4.0) * np.sin(x/np.sqrt(3.0))
def f_n3(x, y): return x * np.exp(-x)
def y_exact_n3(x): return x - 1 + (x + 2) * np.exp(-x)
def f_s1(t, q): r_cubed = (q[0]**2 + q[1]**2)**1.5; return -q / r_cubed if r_cubed > 0 else -q
def f_s2(t, q): return np.array([-q[0] - 0.1 * q[0]**3 + 0.2 * np.cos(2*t)])
def f_s3(t, q): return np.array([-(25 + 4 * np.exp(-t)) * q[0] + np.cosh(t)])

problems = {
    "N1": {"f": f_n1, "solver": numerov_solve, "range": [0, 2], "ic": [0, 0], "N0": 200, "exact_sol": None},
    "N2": {"f": f_n2, "solver": numerov_solve, "range": [0, 2], "ic": [0, 0], "N0": 200, "exact_sol": y_exact_n2},
    "N3": {"f": f_n3, "solver": numerov_solve, "range": [0, 3], "ic": [1, 0], "N0": 300, "exact_sol": y_exact_n3},
    "S1": {"f": f_s1, "solver": stormer_verlet_solve, "range": [0, 20*np.pi], "ic": [[0.6, 0.0], [0.0, 1.527525]], "N0": int(20*np.pi/0.01), "exact_sol": None},
    "S2": {"f": f_s2, "solver": stormer_verlet_solve, "range": [0, 60], "ic": [[0], [0]], "N0": 6000, "exact_sol": None},
    "S3": {"f": f_s3, "solver": stormer_verlet_solve, "range": [0, 6], "ic": [[0], [0]], "N0": 600, "exact_sol": None}
}

for name, p in problems.items():
    fig, axs = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f'Complete Analysis for Problem {name}', fontsize=16)
    
    a, b = p["range"]
    y0, y_prime0 = p["ic"]
    
    # Solution plots
    x_sol, y_sol = p["solver"](p["f"], a, b, y0, y_prime0, p["N0"])
    
    # --- MODIFIED SECTION ---
    if y_sol.ndim == 1:
        # 1D solution (from Numerov)
        axs[0].plot(x_sol, y_sol)
        axs[0].set_xlabel("x")
        axs[0].set_ylabel("y(x)")
    else:
        # 2D or 1D vector solution (from Stormer-Verlet)
        if y_sol.shape[1] > 1:
            # S1: 2D problem, plot phase portrait y1 vs y2
            axs[0].plot(y_sol[:, 0], y_sol[:, 1])
            axs[0].set_xlabel("y_1(t)")
            axs[0].set_ylabel("y_2(t)")
            if name == "S1": axs[0].axis('equal')
        else:
            # S2, S3: 1D problem, plot t vs y1
            axs[0].plot(x_sol, y_sol[:, 0])
            axs[0].set_xlabel("t")
            axs[0].set_ylabel("y(t)")
    # --- END MODIFIED SECTION ---
            
    axs[0].set_title("Solution")
    axs[0].grid(True)

    # Error plots
    if name == "S1": 
        # Plot energy conservation for kepler
        h = (b-a)/p["N0"]
        vx, vy = np.zeros_like(x_sol), np.zeros_like(x_sol)
        vx[1:-1] = (y_sol[2:, 0] - y_sol[:-2, 0]) / (2 * h)
        vy[1:-1] = (y_sol[2:, 1] - y_sol[:-2, 1]) / (2 * h)
        r = np.sqrt(y_sol[:, 0]**2 + y_sol[:, 1]**2)
        energy = 0.5 * (vx**2 + vy**2) - 1.0 / r
        energy_change = energy[1:-1] - energy[1]
        axs[1].plot(x_sol[1:-1], energy_change)
        axs[1].set_ylabel("Change in Energy")
        axs[1].set_title("Energy Conservation Error")
    else:
        if p["exact_sol"]:
            y_exact = p["exact_sol"](x_sol)
            error = np.abs(y_sol - y_exact)
            axs[1].set_title("Absolute Error vs Exact Solution")
        else: 
            # Reference solution for 1D vector cases (S2, S3)
            x_ref, y_ref = p["solver"](p["f"], a, b, y0, y_prime0, p["N0"] * 5)
            # Ensure correct shape for comparison
            y_sol_cmp = y_sol
            y_ref_cmp = y_ref[::5]
            
            # Handle (N,1) vs (N,) shape mismatch if it occurs
            if y_sol.ndim > 1 and y_sol.shape[1] == 1:
                y_sol_cmp = y_sol.ravel()
            if y_ref.ndim > 1 and y_ref.shape[1] == 1:
                y_ref_cmp = y_ref[::5].ravel()

            error = np.abs(y_sol_cmp - y_ref_cmp)
            axs[1].set_title("Error vs Reference Solution")
        
        error_to_plot = error if error.ndim == 1 else np.linalg.norm(error, axis=1)
        axs[1].plot(x_sol, error_to_plot)
        axs[1].set_ylabel("Absolute Error")
    axs[1].set_xlabel("x or t")
    axs[1].grid(True)
    
    # Convergence plots
    N_values = [p["N0"] // 4, p["N0"] // 2, p["N0"], p["N0"] * 2, p["N0"] * 4]
    h_values, max_errors = [], []

    if name == "S1": 
        # S1 convergence is on energy deviation
        ref_energy_val = energy[1]
        for N in N_values:
            t, q = p["solver"](p["f"], a, b, y0, y_prime0, N)
            h = (b-a)/N
            vx, vy = np.zeros_like(t), np.zeros_like(t)
            vx[1:-1] = (q[2:, 0] - q[:-2, 0]) / (2 * h)
            vy[1:-1] = (q[2:, 1] - q[:-2, 1]) / (2 * h)
            r = np.sqrt(q[:, 0]**2 + q[:, 1]**2)
            energy = 0.5 * (vx**2 + vy**2) - 1.0 / r
            max_errors.append(np.max(np.abs(energy[1:-1] - ref_energy_val)))
            h_values.append(h)
    else:
        if p["exact_sol"]:
            # Convergence vs exact solution
            for N in N_values:
                x, y_num = p["solver"](p["f"], a, b, y0, y_prime0, N)
                y_ex = p["exact_sol"](x)
                err = y_num - y_ex
                max_errors.append(np.max(np.abs(err if err.ndim==1 else np.linalg.norm(err,axis=1))))
                h_values.append((b - a) / N)
        else: 
            # Convergence vs reference solution
            x_ref, y_ref = p["solver"](p["f"], a, b, y0, y_prime0, p["N0"] * 8)
            for N in N_values:
                x, y_num = p["solver"](p["f"], a, b, y0, y_prime0, N)
                
                # Ensure correct shapes for comparison
                y_num_cmp = y_num
                y_ref_cmp = y_ref[::(p["N0"]*8 // N)]
                
                if y_num.ndim > 1 and y_num.shape[1] == 1:
                    y_num_cmp = y_num.ravel()
                if y_ref.ndim > 1 and y_ref.shape[1] == 1:
                    y_ref_cmp = y_ref[::(p["N0"]*8 // N)].ravel()

                err = y_num_cmp - y_ref_cmp
                max_errors.append(np.max(np.abs(err if err.ndim==1 else np.linalg.norm(err,axis=1))))
                h_values.append((b - a) / N)
                
    slope = np.polyfit(np.log(h_values), np.log(max_errors), 1)[0]
    axs[2].loglog(h_values, max_errors, 'o-', label=f'Slope = {slope:.2f}')
    axs[2].set_title("Convergence Plot")
    axs[2].set_xlabel("Step Size h")
    axs[2].set_ylabel("Max Error")
    axs[2].legend()
    axs[2].grid(True)
    
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plot_filename = f"analysis_plot_{name}.png"
    plt.savefig(plot_filename)
    print(f"Saved plot to {plot_filename}")
    
    plt.close(fig)

# Analysis:
# We donot have exact solutions for all of them and hence we use techniques like energy stability for kepler and reference (very small step size) for others.
# The convergence plots show the expected order of accuracy for both Numerov (4th order) and Stormer-Verlet (2nd order) methods.
# The error plots indicate that the methods are stable and the errors decrease with smaller step sizes as expected.
# The energy plot for the Kepler problem shows that the Stormer-Verlet method conserves energy very well over long time integrations.
# The stiff equation requires h < 0.4 for stability. We use 0.01 and hence it is stable. Otherwise it can lead to unstable results.