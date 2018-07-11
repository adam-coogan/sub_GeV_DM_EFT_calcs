from hazma.field_theory_helper_functions.common_functions import minkowski_dot
from hazma.parameters import GF, qe, Vus
from hazma.parameters import charged_pion_mass as mpi
from hazma.parameters import neutral_kaon_mass as mk0


def mat_elt_sqrd_k0_to_pi_l_nu_g(klist, ml):
    """
    Parameters
    ----------
    klist
        List of 4-momenta for the k0, pi, l, nu and g respectively.
    """
    P = klist[0]
    p = klist[1]
    q1 = klist[2]
    q2 = klist[3]
    k = klist[4]

    kDOTp = minkowski_dot(k, p)
    q1DOTq2 = minkowski_dot(q1, q2)
    pDOTP = minkowski_dot(p, P)
    kDOTP = minkowski_dot(k, P)
    pDOTq1 = minkowski_dot(p, q1)
    pDOTq2 = minkowski_dot(p, q2)
    kDOTq1 = minkowski_dot(k, q1)
    kDOTq2 = minkowski_dot(k, q2)
    PDOTq1 = minkowski_dot(P, q1)
    PDOTq2 = minkowski_dot(P, q2)

    return (4.*qe**2*GF**2 *
            (-2.*kDOTp**3*(-kDOTq1 + ml**2) *
             (pDOTq2 + PDOTq2) + kDOTq1**2*mpi**2 *
             (-2.*kDOTq2*(pDOTq1 + PDOTq1) - 2.*pDOTq1*pDOTq2 -
              2.*PDOTq1*pDOTq2 - 2.*pDOTq1*PDOTq2 - 2.*PDOTq1*PDOTq2 -
              2.*kDOTq1*(kDOTq2 + pDOTq2 + PDOTq2) + 2.*kDOTP*q1DOTq2 +
              mk0**2*q1DOTq2 + mpi**2*q1DOTq2 + 2.*pDOTP*q1DOTq2) +
             kDOTp**2*(kDOTq2*(-kDOTq1 + ml**2)*(mk0**2 + mpi**2 + 2.*pDOTP) +
                       4.*kDOTq1**2*pDOTq2 - 2.*kDOTq1*ml**2*pDOTq2 +
                       4.*kDOTq1*pDOTq1*pDOTq2 - 2.*ml**2*pDOTq1*pDOTq2 +
                       2.*kDOTq1*PDOTq1*pDOTq2 - 2.*ml**2*PDOTq1*pDOTq2 +
                       2.*kDOTq1**2*PDOTq2 - 2.*kDOTq1*ml**2*PDOTq2 +
                       4.*kDOTq1*pDOTq1*PDOTq2 - 2.*ml**2*pDOTq1*PDOTq2 +
                       2.*kDOTq1*PDOTq1*PDOTq2 - 2.*ml**2*PDOTq1*PDOTq2 -
                       2.*kDOTP*(-kDOTq1 + ml**2)*(pDOTq2 + PDOTq2) -
                       kDOTq1*mk0**2*q1DOTq2 + mk0**2*ml**2*q1DOTq2 -
                       kDOTq1*mpi**2*q1DOTq2 + ml**2*mpi**2*q1DOTq2 -
                       2.*kDOTq1*pDOTP*q1DOTq2 + 2.*ml**2*pDOTP*q1DOTq2 -
                       2.*kDOTq1*pDOTq1*q1DOTq2) +
             kDOTp*kDOTq1 *
             (2.*kDOTq1**2*pDOTq2 + pDOTq1 *
              (-(kDOTq2*(mk0**2 + mpi**2 + 2.*pDOTP - 2.*pDOTq1 - 2.*PDOTq1)) +
               2.*(2.*PDOTq1*pDOTq2 + 2.*PDOTq1*PDOTq2 +
                   2.*pDOTq1*(pDOTq2 + PDOTq2) +
                   kDOTP*(pDOTq2 + PDOTq2 - q1DOTq2) - mk0**2*q1DOTq2 -
                   mpi**2*q1DOTq2 - 2.*pDOTP*q1DOTq2)) +
              kDOTq1*(-2.*kDOTq2*(mpi**2 + pDOTP - pDOTq1) +
                      (2.*kDOTP + mk0**2 - mpi**2 + 6.*pDOTq1 +
                       2.*PDOTq1)*pDOTq2 -
                      2.*((mpi**2 + pDOTP - 2.*pDOTq1)*PDOTq2 +
                          pDOTP*q1DOTq2))))*Vus**2) / \
        (kDOTp**2*kDOTq1**2)
