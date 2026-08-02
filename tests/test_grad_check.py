import numpy as np
import pytest
from fwat import FwatModel


# number of leading parameters carried in log space per kltype (velocities + rho);
# their finite-difference kernels are scaled by the model to match convert_kl.
_N_LOG = {1: 3, 2: 5, 3: 3}


def user_defined_misfit(weight,cmodel:np.ndarray):
    """
    Example of user-defined misfit function.
    This function computes a simple L2 norm of the model parameters
    as a demonstration. Users can modify this function to implement
    their own misfit calculations based on specific requirements.

    Args:
        cmodel: np.ndarray
            Array of model parameters.

    Returns:
        float
            Computed misfit value.
    """
    # Simple L2 norm as an example misfit
    misfit = np.sum(weight * cmodel**2) * 0.5
    kl = weight * cmodel
    return misfit,kl

def get_base_model(mdtype, kltype):
    if mdtype != 'dtti':
        raise NotImplementedError(f"model type mdtype='{mdtype}' is not implemented")

    G0 = 0.05
    phi = np.pi / 3
    gcp = G0 * np.cos(2 * phi)
    gsp = G0 * np.sin(2 * phi)
    rho = 2.8
    eta = 1.0

    if kltype == 1:
        vs = 3.5
        vp = 1.732 * vs
        user_model = np.array([[vp], [vs], [rho], [gcp], [gsp]], dtype='f4')
    elif kltype == 2:
        vsv = 3.5
        vsh = 3.7
        vph = 1.732 * vsv
        vpv = 1.732 * vsh
        user_model = np.array([[vph], [vpv], [vsh], [vsv], [rho], [eta], [gcp], [gsp]], dtype='f4')
    elif kltype == 3:
        vsv = 3.5
        vsh = 3.7
        vph = 1.732 * vsv
        vpv = 1.732 * vsh
        vp = np.sqrt((2 * vph**2 + vpv**2) / 3)
        vs = np.sqrt((2 * vsh**2 + vsv**2) / 3)
        kappaa = (vph - vpv) / vpv
        kappab = (vsh - vsv) / vsv
        user_model = np.array([[vp], [vs], [rho], [kappaa], [kappab], [eta], [gcp], [gsp]], dtype='f4')
    else:
        raise NotImplementedError(f"dtti kltype={kltype} is not implemented")

    M = FwatModel(filename=None, mdtype=mdtype, kltype=kltype)
    cmodel = M.convert_model(user_model, backward=True)
    return cmodel


def grad_check(mdtype, kltype, delta=1.0e-3):
    """
    Compare analytic kernels (from convert_kl) against finite-difference kernels
    for a given (mdtype, kltype). Returns (direc_names, kl_fd, kl_analytic).
    """
    np.random.seed(42)

    # generate base model + a reference misfit/kernel
    cmodel = get_base_model(mdtype, kltype)
    weight = np.random.rand(cmodel.shape[0], 1)
    _, kl_base = user_defined_misfit(weight, cmodel)

    M = FwatModel(filename=None, mdtype=mdtype, kltype=kltype)

    # user model + analytic kernels
    model = M.convert_model(cmodel)
    _, kl_analytic = M.convert_kl(cmodel, kl_base)

    # finite-difference kernels
    kl_fd = np.zeros_like(kl_analytic)
    for i in range(model.shape[0]):
        model_p = model.copy()
        model_m = model.copy()
        model_p[i, 0] *= (1. + delta)
        model_m[i, 0] *= (1. - delta)

        cmodel_p = M.convert_model(model_p, backward=True)
        cmodel_m = M.convert_model(model_m, backward=True)
        misfit_p, _ = user_defined_misfit(weight, cmodel_p)
        misfit_m, _ = user_defined_misfit(weight, cmodel_m)

        kl_fd[i, 0] = (misfit_p - misfit_m) / (2. * delta * model[i, 0])

    # analytic kernels are w.r.t. log parameters for the leading n_log entries;
    # scale the FD kernels the same way so they are comparable.
    n_log = _N_LOG[kltype]
    kl_fd[:n_log, ...] *= model[:n_log, ...]

    return M.direc_names(), kl_fd, kl_analytic


@pytest.mark.parametrize("kltype", [1, 2, 3])
def test_grad_check_dtti(kltype):
    _, kl_fd, kl_analytic = grad_check("dtti", kltype)

    # guarded relative error; ~0.05% is finite-difference noise, a real
    # derivative bug would be orders of magnitude larger.
    rel_err = np.abs(kl_fd - kl_analytic) / (np.abs(kl_analytic) + 1.e-20)
    assert np.max(rel_err) < 1.e-2, (
        f"dtti kltype={kltype}: max FD gradient relative error "
        f"{np.max(rel_err) * 100:.4f}% exceeds 1%"
    )


def main():
    for kltype in (1, 2, 3):
        direc_names, kl_fd, kl_analytic = grad_check("dtti", kltype)
        diff = np.abs(kl_fd - kl_analytic) / (np.abs(kl_analytic) + 1.e-20) * 100
        for i in range(kl_analytic.shape[0]):
            print(f"Model type: dtti, KL type: {kltype}, Direction: {direc_names[i]}, "
                  f"kl_fd: {kl_fd[i,0]:.6e}, kl_analytic: {kl_analytic[i,0]:.6e}, "
                  f"relative error: {diff[i,0]}%")
        print("")


if __name__ == "__main__":
    main()
