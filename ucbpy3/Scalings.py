

# Montagner and Anderson 1989 Upper-Mantle constants
MA1989 = {'dlnEta_dlnXi': -2.5,
          'dlnVp_dlnVs':   0.5,
          'dlnPhi_dlnXi': -1.5,
          'dlnRho_dlnVs':  0.33}

# Brocher's crustal Vp, Rho fit
def Brocher(Vs):
    Vp = 0.9409 + 2.0947 * Vs - 0.8206 * Vs ** 2 + 0.2683 * Vs ** 3 - 0.0251 * Vs ** 4 
    Rho = 1.6612 * Vp - 0.4721 * Vp ** 2 + 0.0671 * Vp ** 3 - 0.0043 * Vp ** 4 + 0.000106 * Vp ** 5 
    return Vp, Rho
