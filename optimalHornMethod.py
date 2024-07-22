import numpy as np
from scipy.optimize import fsolve

# Constants
c = 299_792_458  # Speed of light in m/s



# Function to solve quartic equation (9-159)
def solve_A(G, lamb, a, b):
    eps_ap = 0.5  # Aperture efficiency (assumed)
    coefficients = [
        1,
        -(a + b),
        3 * G * lamb**2 / (8 * np.pi * eps_ap),
        0,
        -3 * G**2 * lamb**4 / (32 * np.pi**2 * eps_ap**2)
    ]
    roots = np.roots(coefficients)
    # Filter the real roots
    real_roots = [root for root in roots if np.isreal(root)]
    if real_roots:
        return np.real(real_roots[0])
    else:
        raise ValueError("No real solution found for A.")

# Function to calculate remaining dimensions
def calculate_dimensions(A, G, lamb, a, b):
    eps_ap = 0.51  # Aperture efficiency (assumed)
    B = G/(0.51*A*(4*(np.pi/pow(lamb,2))))
    R1 = A**2 / (3 * lamb)
    RH = R1*(A - a) / (A)
    lH = np.sqrt(pow((A/2), 2) + R1**2)
    R2 = B**2 / (2 * lamb)
    RE = R2*(B - b) / (B)
    lE = np.sqrt(pow((B/2), 2) + R2**2)
    return B, R1, RH, lH, R2, RE, lE

def solve_quartic(a, b, c, d, e):
    coeffs = [a, b, c, d, e]
    roots = np.roots(coeffs)
    real_roots = [root.real for root in roots if np.isreal(root)]
    return real_roots

# Step 1: Specify the desired gain, operating wavelength, and waveguide dimensions
G_dB = 20  # Desired gain in dB
G = 10**(G_dB / 10)  # Convert gain to linear scale
# freq = 33.25e9  # Operating frequency in Hz
freq = 22.25e9  # Operating frequency in Hz
lamb = c / freq  # Wavelength in meters
# WR-28 standard
# a1 = 0.007112  # Waveguide dimension a in meters
# b1 = 0.003556  # Waveguide dimension b in meters
a1 = 0.010668 # Waveguide dimension a in meters
b1 = 0.004318  # Waveguide dimension b in meters


# Solve the equation 1x^4 - ax^3 + 0 + (3*b*G*(Lam)^2 / 8*pi*E_ap)x - 3G^2 Lam^4 / 32pi^2 E_ap^2 = 0
eps_ap = 0.51
lamb = c/freq
a = 1
b = -a1
c1 = 0
d = (3*b1*G*pow(lamb,2))/(8*np.pi*eps_ap)
e = -(3*pow(G,2)*pow(lamb,4))/(32*pow(np.pi, 2)*pow(eps_ap,2))

ARoots = solve_quartic(a, b, c1, d, e)
print("Roots:", ARoots)

# # Step 2: Solve for A
# A = solve_A(G, lamb, a, b)

for roots in ARoots:
    if roots > 0:
        # Step 3: Find the remaining horn dimensions
        B, R1, RH, lH, R2, RE, lE = calculate_dimensions(roots, G, lamb, a1, b1)   

        # Step 4: Verify the solution
        s = (1/8)*pow((B/lamb),2)*(1/(R2/lamb))
        t = (1/8)*pow((roots/lamb),2)*(1/(R1/lamb))

        if round(RE, 4) == round(RH, 4) and round(s, 2) == 0.25 and round(t, 3) == 0.375:
            print(f"Final results:")
            print(f"A = {roots*1000:.6f} mm")
            print(f"B = {B*1000:.6f} mm")
            print(f"R1 = {R1*1000:.6f} mm")
            print(f"RH = {RH*1000:.6f} mm")
            print(f"lH = {lH*1000:.6f} mm")
            print(f"R2 = {R2*1000:.6f} mm")
            print(f"RE = {RE*1000:.6f} mm")
            print(f"lE = {lE*1000:.6f} mm")
            # print(f"s = {s:.6f}")
            # print(f"t = {t:.6f}")

            print(f"Resulting Gain: {(eps_ap*(4*np.pi/pow(lamb,2)*roots*B)):.6f}")