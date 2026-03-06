import numpy as np

def percus_yevick_Sq(q, phi, d):
    """
    Compute Percus-Yevick hard sphere structure factor S(q)
    
    Parameters:
        q   : array of scattering vector magnitudes in 1/nm or 1/Å
        phi : volume fraction (0 < phi < 0.49)
        d   : hard sphere diameter (same units as 1/q)
    Returns:
        S(q)
    """
    x = q * d
    
    # PY coefficients
    alpha = (1 + 2*phi)**2 / (1 - phi)**4
    beta = -6*phi * (1 + phi/2)**2 / (1 - phi)**4
    gamma = phi*(1 + 2*phi)**2 / (2*(1 - phi)**4)

    # helper functions
    f1 = (np.sin(x) - x*np.cos(x)) / x**3
    f2 = (2*x*np.sin(x) - (x**2 - 2)*np.cos(x) - 2) / x**4
    f3 = (-x**4*np.cos(x) + 4*((3*x**2 - 6)*np.cos(x) + 
          (x**3 - 6*x)*np.sin(x) + 6)) / x**6

    c_q = alpha*f1 + beta*f2 + gamma*f3
    S_q = 1 / (1 - 24*phi*c_q)  # final structure factor

    return S_q


# Example usage
if __name__ == "__main__":
    q = np.linspace(0.01, 30, 500)  # q range
    phi = 0.2                      # volume fraction
    d = 10                         # diameter (same units as q^-1)
    S = percus_yevick_Sq(q, phi, d)

    # plot if needed
    import matplotlib.pyplot as plt
    plt.plot(q, S)
    plt.xlabel("q")
    plt.ylabel("S(q)")
    plt.title("Percus-Yevick Hard Sphere Structure Factor")
    plt.show()
