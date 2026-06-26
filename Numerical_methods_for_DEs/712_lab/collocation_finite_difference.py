import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 12, 'figure.figsize': (10, 6)})

def solve_problem_1():
    print("Solving Problem 1 (Non-Uniform Grid FVM)")
    p = lambda x: 1 + x**2
    q = lambda x: np.exp(x)
    f = lambda x: (1 + x**2) * np.pi**2 * np.sin(np.pi * x) - \
                  2 * np.pi * x * np.cos(np.pi * x) - 2 * x + \
                  np.exp(x) * (np.sin(np.pi * x) + x)
    
    y_exact_func = lambda x: np.sin(np.pi * x) + x
    
    bc_0_val = np.pi + 1  # y' + y = bc_0 at x=0
    bc_1_val = -2 * np.pi + 1 # 2y' - y = bc_1 at x=1

    N_values = [64, 128, 256, 512]
    errors = []
    all_h = []

    for N in N_values:
        s = np.linspace(0, 1, N + 1)
        x = 0.5 * (1 - np.cos(np.pi * s))
        
        h = np.diff(x) # h[i] = x[i+1] - x[i]. (length N)
                       # h[0]=h_1, h[N-1]=h_N
        
        x_mid = 0.5 * (x[1:] + x[:-1]) # x_mid[i] = x_{i+1/2}. (length N)
                                       # x_mid[0]=x_{1/2}, x_mid[N-1]=x_{N-1/2}
        
        H = 0.5 * (x[2:] - x[:-2]) # H[i] = H_{i+1}. (length N-1)
                                  # H[0] = (x[2]-x[0])/2 = H_1

        A = np.zeros((N + 1, N + 1))
        b = np.zeros(N + 1)
        
        # p_mid[i] = p(x_{i+1/2})
        p_mid = p(x_mid) 

        for i in range(1, N):
            p_plus  = p_mid[i]   # p(x_{i+1/2})
            p_minus = p_mid[i-1] # p(x_{i-1/2})
            h_plus  = h[i]       # h_{i+1}
            h_minus = h[i-1]     # h_i
            H_i     = H[i-1]     # H_i
            
            A[i, i-1] = -p_minus / h_minus
            A[i, i]   = p_plus / h_plus + p_minus / h_minus + q(x[i]) * H_i
            A[i, i+1] = -p_plus / h_plus
            b[i]      = f(x[i]) * H_i
        
        A[0, 0] = p_mid[0] / h[0] + p(x[0]) + q(x[0]) * h[0] / 2.0
        A[0, 1] = -p_mid[0] / h[0]
        b[0]    = f(x[0]) * h[0] / 2.0 + p(x[0]) * bc_0_val

        A[N, N-1] = -p_mid[N-1] / h[N-1]
        A[N, N]   = p(x[N]) / 2.0 + p_mid[N-1] / h[N-1] + q(x[N]) * h[N-1] / 2.0
        b[N]      = f(x[N]) * h[N-1] / 2.0 + p(x[N]) * bc_1_val / 2.0
        
        y_num = np.linalg.solve(A, b)
        y_exact = y_exact_func(x)
        
        error = np.linalg.norm(y_num - y_exact, np.inf)
        errors.append(error)
        all_h.append(np.max(h))
        print(f"N = {N:4d}, Max h = {np.max(h):.2e}, L_inf Error = {error: .2e}")

    h_values = np.array(all_h)
    errors = np.array(errors)
    order = np.log2(errors[:-1] / errors[1:])
    print(f"Observed order of convergence: {order}")

    plt.figure()
    plt.loglog(h_values, errors, 'bo-', label='L-inf Error (Non-Uniform)')
    plt.loglog(h_values, (h_values)**2, 'r--', label=r'$O(h^2)$')
    plt.title('Problem 1: FDM (Non-Uniform Grid) Error vs Max h')
    plt.xlabel('Max h')
    plt.ylabel('L-infinity Error')
    plt.legend()
    plt.grid(True, which="both", ls="--")
    
    filename = 'problem_1_error_convergence.png'
    plt.savefig(filename)
    plt.close()
    print(f"Saved plot to {filename}")


def bratu_newton(N, lambda_val, tol=1e-10, max_iter=100):
    h = 1.0 / N
    x = np.linspace(0, 1, N+1)
    
    y = np.zeros(N+1)
    cvg = False
    
    for k in range(max_iter):
        F = np.zeros(N-1, dtype=np.float64)
        J = np.zeros((N-1, N-1), dtype=np.float64)
        
        for i in range(N-1):
            # F(y) = -y'' - λe^y
            F[i] = -(y[i] - 2*y[i+1] + y[i+2])/h**2 - lambda_val*np.exp(y[i+1])
            
            # Jacobian entries
            J[i,i] = 2/h**2 - lambda_val*np.exp(y[i+1])
            if i < N-2:
                J[i,i+1] = -1/h**2  
            if i > 0:
                J[i,i-1] = -1/h**2

        # Solve J*dy = F
        try:
            dy = np.linalg.solve(J, -F)
        except np.linalg.LinAlgError:
            dy = np.zeros(N-1, dtype=np.float64)
            cvg = False
            break

        # Update solution
        y[1:-1] += dy
        
        # Check convergence
        if np.max(np.abs(dy)) < tol:
            cvg = True
            break
    return x, y , cvg

def solve_problem_2():
    print("\nSolving Problem 2 :")
    # Range of lambda values
    lambda_values = np.linspace(1, 4, 50)
    N_values = [100, 200, 400, 800]
    max_norms = {N: [] for N in N_values}

    # Compute solutions for different N and lambda
    lambda_crit = None
    for N in N_values:
        for lambda_val in lambda_values:
            x, y, cvg = bratu_newton(N, lambda_val)
            if(cvg):
                max_norms[N].append(np.max(np.abs(y)))
            else:
                lambda_crit = lambda_val
                break


    # Plot results
    plt.figure(figsize=(10, 6))
    for N in N_values:
        plt.plot(lambda_values[:len(max_norms[N])], max_norms[N], 
                label=f'N={N}', marker='o', markersize=3)
        
    # plot lambda critical line
    if lambda_crit is not None:
        plt.axvline(x=lambda_crit, color='k', linestyle='--', label='Critical λ')

    plt.xlabel('λ')
    plt.ylabel('||y||∞')
    plt.title('Bratu Problem: Maximum Norm vs λ')
    plt.legend()
    plt.grid(True)
    
    filename = 'problem_2_bratu_norm.png'
    plt.savefig(filename)
    plt.close()
    print(f"Saved plot to {filename}")

def solve_problem_3():
    print("\nSolving Problem 3 :")
    
    y_exact_func = lambda x, eps: (np.exp(x / eps) - 1) / (np.exp(1 / eps) - 1)
    
    N = 100
    h = 1.0 / N
    x = np.linspace(0, 1, N + 1)
    
    def solve_scheme(eps, scheme_type):
        A = np.zeros((N - 1, N - 1))
        b = np.zeros(N - 1)
        
        C_diff = eps / h**2
        C_conv = 1.0 / (2 * h)
        
        for i in range(N - 1):
            if scheme_type == 'central':
                A[i, i] = 2.0 * C_diff
                if i > 0:
                    A[i, i-1] = -C_diff - C_conv
                if i < N - 2:
                    A[i, i+1] = -C_diff + C_conv
            
            elif scheme_type == 'upwind1':
                A[i, i] = 2.0 * C_diff + 1.0 / h
                if i > 0:
                    A[i, i-1] = -C_diff - 1.0 / h
                if i < N - 2:
                    A[i, i+1] = -C_diff
            
            elif scheme_type == 'upwind2':
                if i == 0:
                    A[i, i] = 2.0 * C_diff + 1.0 / h
                    A[i, i+1] = -C_diff
                else: 
                    A[i, i] = 2.0 * C_diff + 3.0 / (2*h)
                    A[i, i-1] = -C_diff - 4.0 / (2*h)
                    if i > 1:
                        A[i, i-2] = 1.0 / (2*h)
            
        if scheme_type == 'central':
            b[0] = 0
            b[N-2] = (C_diff - C_conv) * 1.0
            
        elif scheme_type == 'upwind1':
            b[0] = 0
            b[N-2] = C_diff * 1.0
        
        elif scheme_type == 'upwind2':
            b[0] = 0
            b[1] = 0
            b[N-2] = C_diff * 1.0
            
        y_interior = np.linalg.solve(A, b)
        return np.concatenate(([0], y_interior, [1]))

    epsilons = [1e-2, 1e-3, 1e-4]
    
    for eps in epsilons:
        y_cd = solve_scheme(eps, 'central')
        y_u1 = solve_scheme(eps, 'upwind1')
        y_u2 = solve_scheme(eps, 'upwind2')
        y_ex = y_exact_func(x, eps)
        
        plt.figure()
        plt.plot(x, y_ex, 'k-', linewidth=3, label='Exact Solution')
        plt.plot(x, y_cd, 'r:o', markersize=5, label='Central Difference')
        plt.plot(x, y_u1, 'g--^', markersize=5, label='1st-Order Upwind')
        plt.plot(x, y_u2, 'b-.s', markersize=5, label='2nd-Order Upwind')
        
        peclet = h / (2 * eps)
        plt.title(f'Problem 3: Convection-Diffusion (N={N}, $\epsilon$={eps:.1e}, $P_h$={peclet:.1f})')
        plt.xlabel('x')
        plt.ylabel('y(x)')
        plt.legend()
        plt.grid(True)
        
        filename = f'problem_3_solution_eps_{eps:.0e}.png'
        plt.savefig(filename)
        plt.close()
        print(f"Saved plot to {filename}")

def poly_basis_der(degree, deriv_order, x):
    if deriv_order > degree:
        return 0.0
    
    if deriv_order == 0:
        return x**degree
    
    coeff = 1.0
    for i in range(deriv_order):
        coeff *= (degree - i)
        
    return coeff * (x**(degree - deriv_order))

def solve_problem_7():
    print("\nSolving Problem 7 :")
    
    N = 3
    A = np.zeros((N + 1, N + 1))
    b = np.zeros(N + 1)
    
    xc = [1.0 / 3.0, 2.0 / 3.0]
    
    def L(j, x):
        y_val   = poly_basis_der(j, 0, x)
        y_p_val = poly_basis_der(j, 1, x)
        y_pp_val= poly_basis_der(j, 2, x)
        return -(1 + x) * y_pp_val + 2 * x * y_p_val + y_val

    for j in range(N + 1):
        A[0, j] = poly_basis_der(j, 0, 0.0)
    b[0] = 0.0
    
    for j in range(N + 1):
        A[1, j] = poly_basis_der(j, 0, 1.0)
    b[1] = 1.0
    
    for j in range(N + 1):
        A[2, j] = L(j, xc[0])
    b[2] = np.sin(np.pi * xc[0])
    
    for j in range(N + 1):
        A[3, j] = L(j, xc[1])
    b[3] = np.sin(np.pi * xc[1])
    
    try:
        a = np.linalg.solve(A, b)
        print(f"Coefficients (a0, a1, a2, a3): {a}")
        
        x_plot = np.linspace(0, 1, 200)
        y_N = np.polyval(a[::-1], x_plot)
        
        P = np.poly1d(a[::-1])
        P_p = P.deriv(1)
        P_pp = P.deriv(2)
        
        y_N_p = P_p(x_plot)
        y_N_pp = P_pp(x_plot)
        
        L_yN = -(1 + x_plot) * y_N_pp + 2 * x_plot * y_N_p + y_N
        f_x = np.sin(np.pi * x_plot)
        R = L_yN - f_x
        
        plt.figure()
        plt.plot(x_plot, y_N, 'b-', label=r'$y_N(x)$ (Cubic Approx.)')
        plt.plot([0, 1], [0, 1], 'ro', label='BCs')
        plt.plot(xc, np.polyval(a[::-1], xc), 'gs', label='Collocation Pts')
        plt.title('Problem 7: Cubic Collocation Solution')
        plt.xlabel('x')
        plt.ylabel('y(x)')
        plt.legend()
        plt.grid(True)
        
        filename_sol = 'problem_7_solution.png'
        plt.savefig(filename_sol)
        plt.close()
        print(f"Saved plot to {filename_sol}")
        
        plt.figure()
        plt.plot(x_plot, R, 'r-', label=r'$R(x) = L(y_N) - f(x)$')
        plt.plot(xc, [0, 0], 'gs', label='Collocation Pts (R=0)')
        plt.title('Problem 7: Residual')
        plt.xlabel('x')
        plt.ylabel('Residual R(x)')
        plt.legend()
        plt.grid(True)
        
        filename_res = 'problem_7_residual.png'
        plt.savefig(filename_res)
        plt.close()
        print(f"Saved plot to {filename_res}")
        
    except np.linalg.LinAlgError:
        print("Error: The matrix A is singular. Cannot solve.")

def solve_problem_8():
    print("\nSolving Problem 8 :")
    
    N = 4
    A = np.zeros((N + 1, N + 1))
    b = np.zeros(N + 1)
    
    f = lambda x: -( 2*x*(np.pi*np.cos(np.pi*x) + 1 - 2*x) + \
                     (1+x**2)*(-np.pi**2*np.sin(np.pi*x) - 2) ) + \
                  (1+x)*(np.sin(np.pi*x) + x - x**2)
    
    xc = [0.5 * (1 + np.cos(i * np.pi / 4.0)) for i in range(1, 4)]
    
    def L(j, x):
        y_val   = poly_basis_der(j, 0, x)
        y_p_val = poly_basis_der(j, 1, x)
        y_pp_val= poly_basis_der(j, 2, x)
        return -(2 * x * y_p_val + (1 + x**2) * y_pp_val) + (1 + x) * y_val
        
    for j in range(N + 1):
        A[0, j] = poly_basis_der(j, 1, 0.0) + poly_basis_der(j, 0, 0.0)
    b[0] = np.pi + 1
    
    for j in range(N + 1):
        A[1, j] = 2 * poly_basis_der(j, 1, 1.0) - poly_basis_der(j, 0, 1.0)
    b[1] = -2 * np.pi - 2
    
    for j in range(N + 1):
        A[2, j] = L(j, xc[0])
    b[2] = f(xc[0])
    
    for j in range(N + 1):
        A[3, j] = L(j, xc[1])
    b[3] = f(xc[1])
    
    for j in range(N + 1):
        A[4, j] = L(j, xc[2])
    b[4] = f(xc[2])
    
    try:
        a = np.linalg.solve(A, b)
        print("Problem 8 Results :")
        print(f"Coefficients (a0...a4): {a}")
        
        x_plot = np.linspace(0, 1, 200)
        
        P = np.poly1d(a[::-1])
        y_N = P(x_plot)
        
        P_p = P.deriv(1)
        P_pp = P.deriv(2)
        
        y_N_p = P_p(x_plot)
        y_N_pp = P_pp(x_plot)
        
        L_yN = -(2 * x_plot * y_N_p + (1 + x_plot**2) * y_N_pp) + (1 + x_plot) * y_N
        f_plot = f(x_plot)
        R = L_yN - f_plot
        
        max_abs_residual = np.max(np.abs(R))
        print(f"Maximum absolute residual: {max_abs_residual:.2e}")

        plt.figure()
        plt.plot(x_plot, y_N, 'b-', label=r'$y_N(x)$ (Quartic Approx.)')
        plt.plot(xc, P(xc), 'gs', label='Collocation Pts')
        
        y_ex = np.sin(np.pi * x_plot) + x_plot - x_plot**2
        plt.plot(x_plot, y_ex, 'k--', label='Exact Solution (for comparison)')
        
        plt.title('Problem 8: Quartic Collocation Solution')
        plt.xlabel('x')
        plt.ylabel('y(x)')
        plt.legend()
        plt.grid(True)
        
        filename_sol = 'problem_8_solution.png'
        plt.savefig(filename_sol)
        plt.close()
        print(f"Saved plot to {filename_sol}")
        
        plt.figure()
        plt.plot(x_plot, R, 'r-', label=r'$R(x) = L(y_N) - f(x)$')
        plt.plot(xc, [0, 0, 0], 'gs', label='Collocation Pts (R=0)')
        plt.title('Problem 8: Residual')
        plt.xlabel('x')
        plt.ylabel('Residual R(x)')
        plt.legend()
        plt.grid(True)
        
        filename_res = 'problem_8_residual.png'
        plt.savefig(filename_res)
        plt.close()
        print(f"Saved plot to {filename_res}")

    except np.linalg.LinAlgError:
        print("Error: The matrix A is singular. Cannot solve.")

def solve_problem_9():
    print("\nSolving Problem 9 :")
    
    N = 7
    A = np.zeros((N + 1, N + 1))
    b = np.zeros(N + 1)
    
    EI = lambda x: 1 + x
    q = lambda x: 2 + x
    w = lambda x: 12 - 120*x - 240*x**2 + 2*x**3 + x**4 - 2*x**5 - x**6
    
    xc = [i / 5.0 for i in range(1, 5)]
    
    def L(j, x):
        y_val   = poly_basis_der(j, 0, x)
        y_ppp_val = poly_basis_der(j, 3, x)
        y_pppp_val = poly_basis_der(j, 4, x)
        
        return 2 * y_ppp_val + (1 + x) * y_pppp_val + (2 + x) * y_val

    for j in range(N + 1):
        A[0, j] = poly_basis_der(j, 0, 0.0)
    b[0] = 0.0
    
    for j in range(N + 1):
        A[1, j] = poly_basis_der(j, 1, 0.0)
    b[1] = 0.0
    
    for j in range(N + 1):
        A[2, j] = poly_basis_der(j, 2, 1.0)
    b[2] = 2.0
    
    for j in range(N + 1):
        A[3, j] = poly_basis_der(j, 3, 1.0)
    b[3] = 18.0
    
    for j in range(N + 1):
        A[4, j] = L(j, xc[0])
    b[4] = w(xc[0])
    
    for j in range(N + 1):
        A[5, j] = L(j, xc[1])
    b[5] = w(xc[1])
    
    for j in range(N + 1):
        A[6, j] = L(j, xc[2])
    b[6] = w(xc[2])
    
    for j in range(N + 1):
        A[7, j] = L(j, xc[3])
    b[7] = w(xc[3])

    try:
        a = np.linalg.solve(A, b)
        print("Problem 9 Results :")
        print("Computed Coefficients (a0...a7):")
        for i, coeff in enumerate(a):
            if np.abs(coeff) < 1e-12:
                print(f"  a{i} = 0.0")
            else:
                print(f"  a{i} = {coeff:.6f}")
        
        x_plot = np.linspace(0, 1, 201)
        
        P = np.poly1d(a[::-1])
        y_N = P(x_plot)
        
        P_ppp = P.deriv(3)
        P_pppp = P.deriv(4)
        
        L_yN = 2 * P_ppp(x_plot) + (1 + x_plot) * P_pppp(x_plot) + (2 + x_plot) * y_N
        w_plot = w(x_plot)
        R = L_yN - w_plot
        
        max_R_idx = np.argmax(np.abs(R))
        print(f"Residual is largest near x = {x_plot[max_R_idx]:.3f} (Value: {R[max_R_idx]:.2e})")

        plt.figure()
        plt.plot(x_plot, y_N, 'b-', label=r'$y_N(x)$ (Degree-7 Approx.)')
        plt.plot(xc, P(xc), 'gs', label='Collocation Pts')
        plt.title('Problem 9: 4th-Order Collocation Solution')
        plt.xlabel('x')
        plt.ylabel('y(x)')
        plt.legend()
        plt.grid(True)
        
        filename_sol = 'problem_9_solution.png'
        plt.savefig(filename_sol)
        plt.close()
        print(f"Saved plot to {filename_sol}")
        
        plt.figure()
        plt.plot(x_plot, R, 'r-', label=r'$R(x) = L(y_N) - w(x)$')
        plt.plot(xc, [0, 0, 0, 0], 'gs', label='Collocation Pts (R=0)')
        plt.title('Problem 9: Residual')
        plt.xlabel('x')
        plt.ylabel('Residual R(x)')
        plt.legend()
        plt.grid(True)
        
        filename_res = 'problem_9_residual.png'
        plt.savefig(filename_res)
        plt.close()
        print(f"Saved plot to {filename_res}")

    except np.linalg.LinAlgError:
        print("Error: The matrix A is singular. Cannot solve.")


if __name__ == "__main__":
    
    solve_problem_1()
    solve_problem_2()
    solve_problem_3()
    
    solve_problem_7()
    solve_problem_8()
    solve_problem_9()