from hazma.field_theory_helper_functions.common_functions import minkowski_dot
from hazma.parameters import GF, qe, Vus
from hazma.parameters import charged_pion_mass as mpi
from hazma.parameters import neutral_kaon_mass as mk0


def mat_elt_sqrd_k0_to_pi_l_nu_g(klist, ml):
    """
    Parameters
    ----------
    klist
        List of 4-momenta for the k0, pi, l and  nu respectively.
    """
    P = klist[0]
    p = klist[1]
    q1 = klist[2]
    q2 = klist[3]

    q1DOTq2 = minkowski_dot(q1, q2)
    pDOTP = minkowski_dot(p, P)
    pDOTq1 = minkowski_dot(p, q1)
    pDOTq2 = minkowski_dot(p, q2)
    PDOTq1 = minkowski_dot(P, q1)
    PDOTq2 = minkowski_dot(P, q2)

    return (8.*GF**2*(2.*PDOTq1*pDOTq2 + 2.*PDOTq1*PDOTq2 +
                      2.*pDOTq1*(pDOTq2 + PDOTq2) - mk0**2*q1DOTq2 -
                      mpi**2*q1DOTq2 - 2.*pDOTP*q1DOTq2)*Vus**2)
