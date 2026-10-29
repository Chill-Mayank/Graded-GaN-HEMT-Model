"""
fermi.py - compact model step 1: the Fermi-Dirac integral of order 1/2.
Source [P]: Nazir et al., IEEE TED 72(9), 2025, Eq. (1) and Fig. 1(b).
"""
import numpy as np
from scipy.integrate import quad
from scipy.special import expit, gamma


def fermi_half_numerical(eta_F):
    """Exact F_1/2(eta_F) by numerical integration [P, Eq. (1)]:
    F_1/2 = (2/sqrt(pi)) * integral_0^inf sqrt(eta) / (1 + exp(eta - eta_F)) d(eta)
    eta_F = (E_F - E_C) / (k_B*T), no units."""

    def integrand(eta):
        # sqrt(eta): density of states;  expit(eta_F - eta) = 1/(1 + exp(eta - eta_F)): occupation
        return np.sqrt(eta) * expit(eta_F - eta)

    # The occupation falls from 1 to 0 near eta = eta_F, so integrate the two sides separately.
    # 60 kT above that, the occupation is below exp(-60) ~ 1e-26: nothing left to add.
    split = max(eta_F, 0.0)
    below, _ = quad(integrand, 0.0, split, epsabs=0, epsrel=1e-10)
    above, _ = quad(integrand, split, split + 60.0, epsabs=0, epsrel=1e-10)
    return 2.0 / np.sqrt(np.pi) * (below + above)

def fermi_half_paper(eta_F, a=9.6, b=2.13):
    """Closed-form F_1/2(eta_F): the paper's modified Aymerich-Humet formula [P, Eq. (2)].
    j = 1/2, a = 9.6, b = 2.13 [P]. Eq. (2) as printed leaves out the 2/sqrt(pi) of Eq. (1),
    so it is multiplied back in here [CALC, checked against Fig. 1(b)]."""
    j = 0.5
    high = (j + 1) * 2 ** (j + 1) / (b + eta_F + np.sqrt((eta_F - b) ** 2 + a)) ** (j + 1)
    low = np.exp(-eta_F) / gamma(j + 1)
    return 2.0 / np.sqrt(np.pi) / (high + low)


if __name__ == "__main__":
    # Self-test [CALC]: the exact value against the two textbook limits
    #   Boltzmann  (eta_F << 0):  F_1/2 -> exp(eta_F)
    #   degenerate (eta_F >> 0):  F_1/2 -> 4/(3*sqrt(pi)) * eta_F**1.5
    print(" eta_F    numerical    Boltzmann   degenerate")
    for eta_F in [-10.0, -5.0, 0.0, 5.0, 10.0]:
        exact = fermi_half_numerical(eta_F)
        boltzmann = np.exp(eta_F)
        degenerate = 4.0 / (3.0 * np.sqrt(np.pi)) * max(eta_F, 0.0) ** 1.5
        print(f"{eta_F:+6.1f}   {exact:.5e}  {boltzmann:.5e}  {degenerate:.5e}")
        # Piece 1b: the paper's closed form against the exact integral [CALC]
    print()
    print(" eta_F    numerical     Eq. (2)     error")
    for eta_F in [-10.0, -5.0, 0.0, 2.0, 3.4, 5.0, 10.0]:
        exact = fermi_half_numerical(eta_F)
        paper = fermi_half_paper(eta_F)
        print(f"{eta_F:+6.1f}   {exact:.5e}  {paper:.5e}  {100 * (paper / exact - 1):+6.2f} %")