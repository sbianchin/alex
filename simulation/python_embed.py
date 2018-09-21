import numpy as np
from scipy import integrate, optimize, interpolate

vec_xx_lepton = np.array([-17.05, -17.05, -13.95, -13.95, -13.95, -13.95,
                 -13.95, -13.95, -10.85, -10.85, -7.75, -10.85,
                 -7.75, -10.85, -7.75, -4.65, -10.85, -4.65, -4.65,
                 -1.55, -10.85, -7.75, -1.55, -7.75, -1.55, 1.55])

vec_yy_lepton = np.array([13.95, 10.85, 10.85, 7.75, 4.65, 1.55, -1.55, -4.65,
                 -4.65, -7.75, -7.75, -10.85, -10.85, -13.95, -13.95,
                 -13.95, -17.05, -17.05, -20.15, -20.15, -23.25, -23.25,
                 -23.25, -26.35, -26.35, -26.35])

best_pars = np.array([ -1.82547076, -25.51895648,  -4.90425332, -65.6441317 ])

def loss(pars, vec_xx_lepton, vec_yy_lepton):
    m1, b1, m2, b2 = pars
    #Assign points to each track
    ssr_1 = (m1*vec_xx_lepton + b1 - vec_yy_lepton)**2 #first SSR
    ssr_2 = (m2*vec_xx_lepton + b2 - vec_yy_lepton)**2 #second SSR
    vec_l1 = vec_xx_lepton[ssr_1 <= ssr_2]
    vec_l2 = vec_xx_lepton[ssr_1 > ssr_2]
    total_ssr = np.sum(ssr_1[ssr_1 <= ssr_2]) + np.sum(ssr_2[ssr_2 < ssr_1])
    return total_ssr
