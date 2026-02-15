#!/usr/bin/env python3
import sys
from pathlib import Path


def emit(payload):
    for key, value in payload.items():
        print(f"{key}={value}")


def vec_to_csv(x):
    return ",".join(f"{float(v):.17g}" for v in x.reshape(-1))


def dims_to_csv(x):
    shape = tuple(int(v) for v in x.shape)
    return ",".join(str(v) for v in shape)


def make_fixture(seed=17, n=10, k=8, n_times=6, n_energies=7):
    import numpy as np

    rng = np.random.default_rng(seed)
    q1, _ = np.linalg.qr(rng.normal(size=(n, k)))

    shift = 2
    p12 = np.zeros((n, n), dtype=float)
    p12[np.arange(n), (np.arange(n) + shift) % n] = 1.0
    q2 = p12 @ q1

    evals = np.arange(k, dtype=float)
    landmarks = np.array([1, 4, 7], dtype=int)  # 1-based for R parity

    return {
        "n_times": int(n_times),
        "n_energies": int(n_energies),
        "evals": evals,
        "phi1": q1,
        "phi2": q2,
        "landmarks_1b": landmarks,
    }


def run_probe():
    import numpy as np
    from pyFM.refine.icp import icp_iteration
    from pyFM.refine.zoomout import zoomout_iteration
    from pyFM.signatures.HKS_functions import auto_HKS
    from pyFM.signatures.WKS_functions import auto_WKS

    fx = make_fixture()
    evals = fx["evals"]
    phi1 = fx["phi1"]
    phi2 = fx["phi2"]
    landmarks_1b = fx["landmarks_1b"]
    landmarks_0b = landmarks_1b - 1

    abs_ev = np.sort(np.abs(evals))
    hks_t = np.geomspace(4 * np.log(10) / abs_ev[-1], 4 * np.log(10) / abs_ev[1], fx["n_times"])

    e_min = np.log(abs_ev[1])
    e_max = np.log(abs_ev[-1])
    wks_sigma = 7 * (e_max - e_min) / fx["n_energies"]
    wks_e = np.linspace(e_min + 2 * wks_sigma, e_max - 2 * wks_sigma, fx["n_energies"])

    hks = auto_HKS(evals, phi1, fx["n_times"], landmarks=None, scaled=True)
    hks_lm = auto_HKS(evals, phi1, fx["n_times"], landmarks=landmarks_0b, scaled=True)
    wks = auto_WKS(evals, phi1, fx["n_energies"], landmarks=None, scaled=True)
    wks_lm = auto_WKS(evals, phi1, fx["n_energies"], landmarks=landmarks_0b, scaled=True)

    k_map = 4
    c0 = np.eye(k_map, dtype=float)
    c_icp = icp_iteration(c0, phi1, phi2, use_adj=False, n_jobs=1)
    c_zoom = zoomout_iteration(c0, phi1, phi2, step=(1, 2), A2=None, n_jobs=1)

    return {
        "status": "ok",
        "n_times": str(fx["n_times"]),
        "n_energies": str(fx["n_energies"]),
        "evals": vec_to_csv(evals),
        "phi1_dims": dims_to_csv(phi1),
        "phi1": vec_to_csv(phi1),
        "phi2_dims": dims_to_csv(phi2),
        "phi2": vec_to_csv(phi2),
        "landmarks_1b": ",".join(str(int(v)) for v in landmarks_1b),
        "hks_time_grid": vec_to_csv(hks_t),
        "wks_energy_grid": vec_to_csv(wks_e),
        "wks_sigma": f"{float(wks_sigma):.17g}",
        "hks_dims": dims_to_csv(hks),
        "hks": vec_to_csv(hks),
        "hks_lm_dims": dims_to_csv(hks_lm),
        "hks_lm": vec_to_csv(hks_lm),
        "wks_dims": dims_to_csv(wks),
        "wks": vec_to_csv(wks),
        "wks_lm_dims": dims_to_csv(wks_lm),
        "wks_lm": vec_to_csv(wks_lm),
        "c0_dims": dims_to_csv(c0),
        "c0": vec_to_csv(c0),
        "icp_once_dims": dims_to_csv(c_icp),
        "icp_once": vec_to_csv(c_icp),
        "zoom_once_dims": dims_to_csv(c_zoom),
        "zoom_once": vec_to_csv(c_zoom),
    }


def main():
    repo_root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(repo_root / "pyFM"))

    try:
        import numpy  # noqa: F401
        import scipy  # noqa: F401
        import pyFM  # noqa: F401
    except Exception as exc:
        emit({"status": "error", "error": f"import_failure:{exc}"})
        return 0

    try:
        emit(run_probe())
        return 0
    except Exception as exc:
        emit({"status": "error", "error": f"runtime_failure:{exc}"})
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
