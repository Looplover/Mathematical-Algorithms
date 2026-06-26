import numpy as np
import math
import matplotlib.pyplot as plt

np.set_printoptions(precision=6, suppress=True, linewidth=100)

# Crank-Nicolson
def solve_q1():

    print("Solving Q1: Crank-Nicolson (Single Mode):")
    
    # Parameters
    dx = 0.1
    dt = 0.005
    T_final = 0.10
    L = 1.0
    
    Nx = int(L / dx)
    Nt = int(T_final / dt)
    
    # Grid
    x = np.linspace(0, L, Nx + 1)
    
    # M = number of interior points
    M = Nx - 1 
    
    # Initial Condition
    u_n = np.sin(math.pi * x)
    
    # Interior points at n=0
    u_interior = u_n[1:-1]
    
    # Crank-Nicolson matrices
    r = dt / (2 * dx**2) # This is r/2 in some notations
    
    # A = (I + r*B), B = tridiag(1, -2, 1) -> A = tridiag(-r, 1+2r, -r)
    A = np.diag((1 + 2 * r) * np.ones(M)) + \
        np.diag(-r * np.ones(M - 1), k=1) + \
        np.diag(-r * np.ones(M - 1), k=-1)
    
    # B = (I - r*B) -> B = tridiag(r, 1-2r, r)
    B = np.diag((1 - 2 * r) * np.ones(M)) + \
        np.diag(r * np.ones(M - 1), k=1) + \
        np.diag(r * np.ones(M - 1), k=-1)

    # Time-stepping loop
    for n in range(Nt):
        # A * u^{n+1} = B * u^n
        rhs = B @ u_interior
        u_interior = np.linalg.solve(A, rhs)
        
    # Re-attach boundaries
    u_final = np.concatenate(([0], u_interior, [0]))
    
    # Exact solution
    u_exact = np.exp(-math.pi**2 * T_final) * np.sin(math.pi * x)
    
    # Error
    max_error = np.max(np.abs(u_final - u_exact))
    
    print(f"Parameters: dx={dx}, dt={dt}, T_final={T_final}")
    print(f"Values u(x_i, 0.10):")
    print(u_final)
    print(f"Max error against exact solution: {max_error:.2e}")
    
    plt.figure(figsize=(10, 6))
    plt.plot(x, u_final, 'bo-', label=f'Numerical (T={T_final})', markersize=4)
    plt.plot(x, u_exact, 'r--', label=f'Exact (T={T_final})')
    plt.title('Q1: Heat Eq. $u_t = u_{xx}$, $IC = \sin(\pi x)$')
    plt.xlabel('$x$')
    plt.ylabel(f'$u(x, {T_final})$')
    plt.legend()
    plt.grid(True)
    
    plot_filename = 'q1_heat_eq_sin_pi_x.png'
    plt.savefig(plot_filename)
    plt.close()
    print(f"Saved plot to {plot_filename}")
    
    print("-" * 50)


def solve_q2():
    
    print("Solving Q2: Crank-Nicolson (Higher Mode):")
    
    # Parameters
    dx = 0.05
    dt = 0.002
    T_final = 0.05
    L = 1.0
    
    Nx = int(L / dx)
    Nt = int(T_final / dt)
    
    # Grid
    x = np.linspace(0, L, Nx + 1)
    
    # M = number of interior points
    M = Nx - 1 
    
    # Initial Condition
    u_n = np.sin(2 * math.pi * x)
    u_interior = u_n[1:-1]
    
    # Crank-Nicolson matrices
    r = dt / (2 * dx**2)
    
    A = np.diag((1 + 2 * r) * np.ones(M)) + \
        np.diag(-r * np.ones(M - 1), k=1) + \
        np.diag(-r * np.ones(M - 1), k=-1)
    
    B = np.diag((1 - 2 * r) * np.ones(M)) + \
        np.diag(r * np.ones(M - 1), k=1) + \
        np.diag(r * np.ones(M - 1), k=-1)

    # Time-stepping loop
    for n in range(Nt):
        rhs = B @ u_interior
        u_interior = np.linalg.solve(A, rhs)
        
    # Re-attach boundaries
    u_final = np.concatenate(([0], u_interior, [0]))
    
    # Exact solution
    u_exact = np.exp(-(2 * math.pi)**2 * T_final) * np.sin(2 * math.pi * x)
    
    # Error
    max_error = np.max(np.abs(u_final - u_exact))
    
    print(f"Parameters: dx={dx}, dt={dt}, T_final={T_final}")
    print(f"Values u(x_i, 0.05):")
    print(u_final)
    print(f"Max error against exact solution: {max_error:.2e}")
    
    plt.figure(figsize=(10, 6))
    plt.plot(x, u_final, 'bo-', label=f'Numerical (T={T_final})', markersize=4)
    plt.plot(x, u_exact, 'r--', label=f'Exact (T={T_final})')
    plt.title('Q2: Heat Eq. $u_t = u_{xx}$, $IC = \sin(2\pi x)$')
    plt.xlabel('$x$')
    plt.ylabel(f'$u(x, {T_final})$')
    plt.legend()
    plt.grid(True)
    
    plot_filename = 'q2_heat_eq_sin_2pi_x.png'
    plt.savefig(plot_filename)
    plt.close()
    print(f"Saved plot to {plot_filename}")
    
    print("-" * 50)


def solve_q3():

    print("--- Solving Q3: Crank-Nicolson (Source Term) ---")
    
    # Parameters
    dx = 0.1
    dt = 0.005
    T_final = 0.10
    L = 1.0
    
    Nx = int(L / dx)
    Nt = int(T_final / dt)
    
    # Grid
    x = np.linspace(0, L, Nx + 1)
    x_interior = x[1:-1]
    
    # M = number of interior points
    M = Nx - 1 
    
    # Initial Condition
    u_n = np.zeros_like(x)
    u_interior = u_n[1:-1]
    
    # Crank-Nicolson matrices
    r = dt / (2 * dx**2)
    
    A = np.diag((1 + 2 * r) * np.ones(M)) + \
        np.diag(-r * np.ones(M - 1), k=1) + \
        np.diag(-r * np.ones(M - 1), k=-1)
    
    B = np.diag((1 - 2 * r) * np.ones(M)) + \
        np.diag(r * np.ones(M - 1), k=1) + \
        np.diag(r * np.ones(M - 1), k=-1)
    
    # Source term function
    def source(t_val):
        return 2 * np.exp(-t_val) * np.sin(math.pi * x_interior)

    # Time-stepping loop
    for n in range(Nt):
        t_n = n * dt
        t_n_plus_1 = (n + 1) * dt
        
        F_n = source(t_n)
        F_n_plus_1 = source(t_n_plus_1)
        
        # A * u^{n+1} = B * u^n + 0.5*dt*(F^{n+1} + F^n)
        rhs = B @ u_interior + 0.5 * dt * (F_n + F_n_plus_1)
        u_interior = np.linalg.solve(A, rhs)
        
    # Re-attach boundaries
    u_final = np.concatenate(([0], u_interior, [0]))
    
    print(f"Parameters: dx={dx}, dt={dt}, T_final={T_final}")
    print(f"Values u(x_i, 0.10):")
    print(u_final)
    
    plt.figure(figsize=(10, 6))
    plt.plot(x, u_final, 'bo-', label=f'Numerical (T={T_final})', markersize=4)
    plt.title('Q3: Heat Eq. $u_t = u_{xx} + 2e^{-t}\sin(\pi x)$')
    plt.xlabel('$x$')
    plt.ylabel(f'$u(x, {T_final})$')
    plt.legend()
    plt.grid(True)
    
    plot_filename = 'q3_heat_eq_source_term.png'
    plt.savefig(plot_filename)
    plt.close()
    print(f"Saved plot to {plot_filename}")
    
    print("-" * 50)

# ADI Methods for 2D Problems
def solve_q4():
    
    print("Solving Q4: ADI (Anisotropic Diffusion):")
    
    # Parameters
    dx = 0.1
    dy = 0.1
    dt = 0.002
    T_final = 0.05
    L = 1.0
    alpha = 2.0
    beta = 0.5
    
    Nx = int(L / dx)
    Ny = int(L / dy)
    Nt = int(T_final / dt)
    
    M = Nx - 1 # Interior points in x
    N = Ny - 1 # Interior points in y
    
    # Grids
    x = np.linspace(0, L, Nx + 1)
    y = np.linspace(0, L, Ny + 1)
    X_int, Y_int = np.meshgrid(x[1:-1], y[1:-1])
    X_full, Y_full = np.meshgrid(x, y)
    
    # Initial Condition
    u_n = np.sin(math.pi * X_full) * np.sin(math.pi * Y_full)
    u_interior = u_n[1:-1, 1:-1]
    
    # Peaceman-Rachford ADI matrices
    rx = alpha * dt / (2 * dx**2)
    ry = beta * dt / (2 * dy**2)
    
    # (I - rx*delta_x^2) u* = (I + ry*delta_y^2) u^n
    LHS_x = np.diag((1 + 2 * rx) * np.ones(M)) + \
            np.diag(-rx * np.ones(M - 1), k=1) + \
            np.diag(-rx * np.ones(M - 1), k=-1)
            
    RHS_y = np.diag((1 - 2 * ry) * np.ones(N)) + \
            np.diag(ry * np.ones(N - 1), k=1) + \
            np.diag(ry * np.ones(N - 1), k=-1)
    
    # (I - ry*delta_y^2) u^{n+1} = (I + rx*delta_x^2) u*
    RHS_x = np.diag((1 - 2 * rx) * np.ones(M)) + \
            np.diag(rx * np.ones(M - 1), k=1) + \
            np.diag(rx * np.ones(M - 1), k=-1)
            
    LHS_y = np.diag((1 + 2 * ry) * np.ones(N)) + \
            np.diag(-ry * np.ones(N - 1), k=1) + \
            np.diag(-ry * np.ones(N - 1), k=-1)

    # Time-stepping loop
    for n in range(Nt):
        # (Implicit in x)
        # We solve LHS_x @ U_star = (U_n @ RHS_y.T)
        rhs_1 = u_interior @ RHS_y.T
        u_star = np.linalg.solve(LHS_x, rhs_1)
        
        # (Implicit in y)
        # We solve (LHS_y @ U_n+1.T).T = RHS_x @ U_star
        # -> U_n+1 @ LHS_y.T = RHS_x @ U_star
        # -> LHS_y @ U_n+1.T = (RHS_x @ U_star).T
        rhs_2 = (RHS_x @ u_star).T
        u_interior_T = np.linalg.solve(LHS_y, rhs_2)
        u_interior = u_interior_T.T
        
    # Re-attach boundaries
    u_final = np.zeros((N + 2, M + 2))
    u_final[1:-1, 1:-1] = u_interior
    
    # Exact solution
    u_exact = np.exp(-(alpha * math.pi**2 + beta * math.pi**2) * T_final) * \
              np.sin(math.pi * X_full) * np.sin(math.pi * Y_full)
              
    max_error = np.max(np.abs(u_final - u_exact))
    
    print(f"Parameters: dx={dx}, dy={dy}, dt={dt}, T_final={T_final}, alpha={alpha}, beta={beta}")
    print(f"Values u(x_i, y_j, 0.05) (center value at 0.5, 0.5): {u_final[Ny//2, Nx//2]:.6f}")
    print("Full u(x_i, y_j, 0.05):\n", u_final) 
    print(f"Max error against exact solution: {max_error:.2e}")
    
    fig, axes = plt.subplots(1, 3, figsize=(20, 6), sharey=True)
    
    # Numerical Solution
    im1 = axes[0].imshow(u_final, extent=[0, L, 0, L], origin='lower', cmap='viridis', aspect='auto')
    axes[0].set_title(f'Numerical Solution $u_h(T={T_final})$')
    axes[0].set_xlabel('$x$')
    axes[0].set_ylabel('$y$')
    fig.colorbar(im1, ax=axes[0], fraction=0.046, pad=0.04)
    
    # Exact Solution
    im2 = axes[1].imshow(u_exact, extent=[0, L, 0, L], origin='lower', cmap='viridis', aspect='auto')
    axes[1].set_title(f'Exact Solution $u(T={T_final})$')
    axes[1].set_xlabel('$x$')
    fig.colorbar(im2, ax=axes[1], fraction=0.046, pad=0.04)
    
    # Error
    error = u_final - u_exact
    vmax = np.max(np.abs(error))
    im3 = axes[2].imshow(error, extent=[0, L, 0, L], origin='lower', cmap='coolwarm', vmin=-vmax, vmax=vmax, aspect='auto')
    axes[2].set_title(f'Error ($u_h - u$), Max: {max_error:.2e}')
    axes[2].set_xlabel('$x$')
    fig.colorbar(im3, ax=axes[2], fraction=0.046, pad=0.04, label='Error')

    fig.suptitle('Q4: ADI Anisotropic Diffusion $u_t = 2u_{xx} + 0.5u_{yy}$', fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95]) 
    
    plot_filename = 'q4_adi_anisotropic.png'
    plt.savefig(plot_filename)
    plt.close()
    print(f"Saved plot to {plot_filename}")
    
    print("-" * 50)


def solve_q5():

    print("Solving Q5: ADI (Linear Reaction-Diffusion):")
    
    # Parameters
    dx = 0.1
    dy = 0.1
    dt = 0.002
    T_final = 0.05
    L = 1.0
    alpha = 1.0
    beta = 1.0
    lam = 1.0
    
    Nx = int(L / dx)
    Ny = int(L / dy)
    Nt = int(T_final / dt)
    
    M = Nx - 1 # Interior points in x
    N = Ny - 1 # Interior points in y
    
    # Grids
    x = np.linspace(0, L, Nx + 1)
    y = np.linspace(0, L, Ny + 1)
    X_full, Y_full = np.meshgrid(x, y)
    
    # Initial Condition
    u_n = np.sin(math.pi * X_full) * np.sin(math.pi * Y_full)
    u_interior = u_n[1:-1, 1:-1]
    
    # ADI matrices with reaction term (C-N on reaction)
    rx = alpha * dt / (2 * dx**2)
    ry = beta * dt / (2 * dy**2)
    r_lambda = lam * dt / 4  # (lambda*dt/2) / 2
    
    # (I - rx*delta_x^2 + r_lambda*I)
    LHS_x = np.diag((1 + 2 * rx + r_lambda) * np.ones(M)) + \
            np.diag(-rx * np.ones(M - 1), k=1) + \
            np.diag(-rx * np.ones(M - 1), k=-1)
            
    # (I + ry*delta_y^2 - r_lambda*I)
    RHS_y = np.diag((1 - 2 * ry - r_lambda) * np.ones(N)) + \
            np.diag(ry * np.ones(N - 1), k=1) + \
            np.diag(ry * np.ones(N - 1), k=-1)
    
    # (I + rx*delta_x^2 - r_lambda*I)
    RHS_x = np.diag((1 - 2 * rx - r_lambda) * np.ones(M)) + \
            np.diag(rx * np.ones(M - 1), k=1) + \
            np.diag(rx * np.ones(M - 1), k=-1)
            
    # (I - ry*delta_y^2 + r_lambda*I)
    LHS_y = np.diag((1 + 2 * ry + r_lambda) * np.ones(N)) + \
            np.diag(-ry * np.ones(N - 1), k=1) + \
            np.diag(-ry * np.ones(N - 1), k=-1)

    # Time-stepping loop
    for n in range(Nt):
        # (Implicit in x)
        rhs_1 = u_interior @ RHS_y.T
        u_star = np.linalg.solve(LHS_x, rhs_1)
        
        # (Implicit in y)
        rhs_2 = (RHS_x @ u_star).T
        u_interior_T = np.linalg.solve(LHS_y, rhs_2)
        u_interior = u_interior_T.T
        
    # Re-attach boundaries
    u_final = np.zeros((N + 2, M + 2))
    u_final[1:-1, 1:-1] = u_interior
    
    # Exact solution
    u_exact = np.exp(-(2 * math.pi**2 + lam) * T_final) * \
              np.sin(math.pi * X_full) * np.sin(math.pi * Y_full)
              
    max_error = np.max(np.abs(u_final - u_exact))
    
    print(f"Parameters: dx={dx}, dy={dy}, dt={dt}, T_final={T_final}, lambda={lam}")
    print(f"Values u(x_i, y_j, 0.05) (center value at 0.5, 0.5): {u_final[Ny//2, Nx//2]:.6f}")
    print("Full u(x_i, y_j, 0.05):\n", u_final) 
    print(f"Max error against exact solution: {max_error:.2e}")
    
    fig, axes = plt.subplots(1, 3, figsize=(20, 6), sharey=True)
    
    # Numerical Solution
    im1 = axes[0].imshow(u_final, extent=[0, L, 0, L], origin='lower', cmap='viridis', aspect='auto')
    axes[0].set_title(f'Numerical Solution $u_h(T={T_final})$')
    axes[0].set_xlabel('$x$')
    axes[0].set_ylabel('$y$')
    fig.colorbar(im1, ax=axes[0], fraction=0.046, pad=0.04)
    
    # Exact Solution
    im2 = axes[1].imshow(u_exact, extent=[0, L, 0, L], origin='lower', cmap='viridis', aspect='auto')
    axes[1].set_title(f'Exact Solution $u(T={T_final})$')
    axes[1].set_xlabel('$x$')
    fig.colorbar(im2, ax=axes[1], fraction=0.046, pad=0.04)
    
    # Error
    error = u_final - u_exact
    vmax = np.max(np.abs(error))
    im3 = axes[2].imshow(error, extent=[0, L, 0, L], origin='lower', cmap='coolwarm', vmin=-vmax, vmax=vmax, aspect='auto')
    axes[2].set_title(f'Error ($u_h - u$), Max: {max_error:.2e}')
    axes[2].set_xlabel('$x$')
    fig.colorbar(im3, ax=axes[2], fraction=0.046, pad=0.04, label='Error')

    fig.suptitle('Q5: ADI Reaction-Diffusion $u_t = u_{xx} + u_{yy} - u$', fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    plot_filename = 'q5_adi_reaction_diffusion.png'
    plt.savefig(plot_filename)
    plt.close()
    print(f"Saved plot to {plot_filename}")
    
    print("-" * 50)


def solve_q6():

    print("Solving Q6: ADI (Source Term):")
    
    # Parameters
    dx = 0.1
    dy = 0.1
    dt = 0.002
    T_final = 0.05
    L = 1.0
    alpha = 1.0
    beta = 1.0
    
    Nx = int(L / dx)
    Ny = int(L / dy)
    Nt = int(T_final / dt)
    
    M = Nx - 1 # Interior points in x
    N = Ny - 1 # Interior points in y
    
    # Grids
    x = np.linspace(0, L, Nx + 1)
    y = np.linspace(0, L, Ny + 1)
    X_int, Y_int = np.meshgrid(x[1:-1], y[1:-1])
    
    # Initial Condition
    u_interior = np.zeros((N, M))
    
    # ADI matrices (same as Q4 with alpha=beta=1)
    rx = alpha * dt / (2 * dx**2)
    ry = beta * dt / (2 * dy**2)
    
    LHS_x = np.diag((1 + 2 * rx) * np.ones(M)) + \
            np.diag(-rx * np.ones(M - 1), k=1) + \
            np.diag(-rx * np.ones(M - 1), k=-1)
            
    RHS_y = np.diag((1 - 2 * ry) * np.ones(N)) + \
            np.diag(ry * np.ones(N - 1), k=1) + \
            np.diag(ry * np.ones(N - 1), k=-1)
            
    RHS_x = np.diag((1 - 2 * rx) * np.ones(M)) + \
            np.diag(rx * np.ones(M - 1), k=1) + \
            np.diag(rx * np.ones(M - 1), k=-1)
            
    LHS_y = np.diag((1 + 2 * ry) * np.ones(N)) + \
            np.diag(-ry * np.ones(N - 1), k=1) + \
            np.diag(-ry * np.ones(N - 1), k=-1)


    # Source term
    def source(t_val):
        return np.exp(-t_val) * np.sin(math.pi * X_int) * np.sin(math.pi * Y_int)

    # Time-stepping loop
    for n in range(Nt):
        t_n = n * dt
        t_n_half = t_n + dt / 2
        
        # Evaluate source at half-step
        F_half = source(t_n_half)
        
        # (Implicit in x)
        # (I-rx*Dxx)u* = (I+ry*Dyy)u^n + dt/2 * F^{n+1/2}
        rhs_1 = (u_interior @ RHS_y.T) + 0.5 * dt * F_half
        u_star = np.linalg.solve(LHS_x, rhs_1)
        
        # (Implicit in y)
        # (I-ry*Dyy)u^{n+1} = (I+rx*Dxx)u* + dt/2 * F^{n+1/2}
        rhs_2 = (RHS_x @ u_star).T + (0.5 * dt * F_half).T
        u_interior_T = np.linalg.solve(LHS_y, rhs_2)
        u_interior = u_interior_T.T
        
    # Re-attach boundaries
    u_final = np.zeros((N + 2, M + 2))
    u_final[1:-1, 1:-1] = u_interior
    
    print(f"Parameters: dx={dx}, dy={dy}, dt={dt}, T_final={T_final}")
    print(f"Values u(x_i, y_j, 0.05) (center value at 0.5, 0.5): {u_final[Ny//2, Nx//2]:.6f}")
    print("Full u(x_i, y_j, 0.05):\n", u_final) 
    
    plt.figure(figsize=(8, 7))
    im = plt.imshow(u_final, extent=[0, L, 0, L], origin='lower', cmap='viridis', aspect='auto')
    plt.colorbar(im, fraction=0.046, pad=0.04, label=f'$u_h(x,y,{T_final})$')
    plt.title('Q6: ADI with Source $u_t = u_{xx} + u_{yy} + F$')
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.tight_layout()
    
    plot_filename = 'q6_adi_source_term.png'
    plt.savefig(plot_filename)
    plt.close()
    print(f"Saved plot to {plot_filename}")
    
    print("-" * 50)

if __name__ == "__main__":
    solve_q1()
    solve_q2()
    solve_q3()
    solve_q4()
    solve_q5()
    solve_q6()