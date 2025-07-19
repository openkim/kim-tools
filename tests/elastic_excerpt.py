import numpy as np
import numpy.typing as npt

from kim_tools.symmetry_util.elasticity import voigt_elast_compon_eqn


def algebraically_reconstruct_matrix(
    elastic_constants: npt.ArrayLike, sgnum: int
) -> npt.ArrayLike:
    """
    From an elasticity matrix in Voigt order and a space group number, reconstruct
    the elasticity matrix based on the algebraic symmetry rules

    Returns:
        Reconstructed matrix
    """
    eqn = voigt_elast_compon_eqn(sgnum)
    reconstructed_matrix = np.zeros((6, 6))
    # first, figure out which constants are unique and extract them
    for i in range(1, 7):
        for j in range(i, 7):
            if (i, j) not in eqn:
                reconstructed_matrix[i - 1, j - 1] = elastic_constants[i - 1, j - 1]

    for dep_comp in eqn:
        if eqn[dep_comp] is None:
            continue
        else:
            dep_comp_zero_based = (dep_comp[0] - 1, dep_comp[1] - 1)
            for coeff, indep_comp in zip(*eqn[dep_comp]):
                indep_comp_zero_based = (indep_comp[0] - 1, indep_comp[1] - 1)
                reconstructed_matrix[dep_comp_zero_based] += (
                    coeff * reconstructed_matrix[indep_comp_zero_based]
                )

    reconstructed_matrix = (
        reconstructed_matrix
        + reconstructed_matrix.T
        - np.diag(reconstructed_matrix.diagonal())
    )

    return reconstructed_matrix
