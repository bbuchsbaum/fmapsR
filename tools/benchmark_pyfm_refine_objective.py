#!/usr/bin/env python3
import sys
from pathlib import Path


def emit(payload):
    for key, value in payload.items():
        print(f"{key}={value}")


def vec_to_csv(x):
    return ",".join(f"{float(v):.17g}" for v in x.reshape(-1))


def int_vec_to_csv(x):
    return ",".join(str(int(v)) for v in x)


def dims_to_csv(x):
    return ",".join(str(int(v)) for v in x.shape)


def mat_to_csv(x):
    return vec_to_csv(x)


def make_refine_fixture(seed=31, n=14, k_full=9, k_init=4, icp_nit=3, zoom_nit=2):
    import numpy as np
    from pyFM.refine.icp import icp_refine
    from pyFM.refine.zoomout import zoomout_refine
    from pyFM.spectral.convert import FM_to_p2p, p2p_to_FM

    rng = np.random.default_rng(seed)

    q1, _ = np.linalg.qr(rng.normal(size=(n, k_full)))
    shift = 3
    p12 = np.zeros((n, n), dtype=float)
    p12[np.arange(n), (np.arange(n) + shift) % n] = 1.0
    q2 = p12 @ q1

    c0 = np.eye(k_init, dtype=float)
    sub_source = np.array([0, 2, 4, 6, 8, 10, 12], dtype=int)
    sub_target = np.array([1, 3, 5, 7, 9, 11, 13], dtype=int)

    icp_sub = icp_refine(
        c0,
        q1[sub_source, :],
        q2[sub_target, :],
        nit=icp_nit,
        use_adj=False,
        return_p2p=False,
        n_jobs=1,
        verbose=False,
    )
    zoom_sub = zoomout_refine(
        c0,
        q1,
        q2,
        nit=zoom_nit,
        step=(1, 1),
        A2=None,
        subsample=(sub_source, sub_target),
        return_p2p=False,
        n_jobs=1,
        verbose=False,
    )

    n_w = 12
    k_w = 6
    w = np.linspace(0.5, 1.7, n_w)
    u1, _ = np.linalg.qr(rng.normal(size=(n_w, k_w)))
    u2, _ = np.linalg.qr(rng.normal(size=(n_w, k_w)))
    phi1_w = u1 / np.sqrt(w)[:, None]
    phi2_w = u2 / np.sqrt(w)[:, None]
    c0_w = np.eye(4, dtype=float)
    p2p_w = FM_to_p2p(c0_w, phi1_w, phi2_w, use_adj=False, n_jobs=1)
    fm_weighted = p2p_to_FM(p2p_w, phi1_w[:, :4], phi2_w[:, :4], A2=w)

    return {
        "icp_nit": int(icp_nit),
        "zoom_nit": int(zoom_nit),
        "step": (1, 1),
        "phi1": q1,
        "phi2": q2,
        "c0": c0,
        "sub_source_0b": sub_source,
        "sub_target_0b": sub_target,
        "icp_sub": icp_sub,
        "zoom_sub": zoom_sub,
        "weights": w,
        "phi1_w": phi1_w,
        "phi2_w": phi2_w,
        "c0_w": c0_w,
        "p2p_w_0b": p2p_w,
        "fm_weighted": fm_weighted,
    }


def make_objective_fixture(seed=47):
    import numpy as np
    from pyFM.optimize.base_functions import (
        LB_commutation,
        LB_commutation_grad,
        descr_preservation,
        descr_preservation_grad,
        energy_func_std,
        grad_energy_std,
        oplist_commutation,
        oplist_commutation_grad,
    )

    rng = np.random.default_rng(seed)

    k1 = 5
    k2 = 4
    p = 3
    n_ops = 3
    w_descr = 0.35
    w_lap = 0.07
    w_comm = 0.65

    c = rng.normal(size=(k2, k1))
    a = rng.normal(size=(k1, p))
    b = rng.normal(size=(k2, p))

    ev1 = np.linspace(1.0, 2.4, k1)
    ev2 = np.linspace(0.9, 2.2, k2)
    ev_sqdiff = (ev2[:, None] - ev1[None, :]) ** 2
    denom = float(ev_sqdiff.sum())
    if denom > 0:
        ev_sqdiff = ev_sqdiff / denom

    op_list = []
    for _ in range(n_ops):
        op1 = rng.normal(size=(k1, k1))
        op2 = rng.normal(size=(k2, k2))
        op_list.append((op1, op2))

    e_descr = w_descr * descr_preservation(c, a, b)
    e_lap = w_lap * LB_commutation(c, ev_sqdiff)
    e_comm = w_comm * oplist_commutation(c, op_list)
    e_total = energy_func_std(
        c.reshape(-1),
        w_descr,
        w_lap,
        w_comm,
        0.0,
        a,
        b,
        op_list,
        [],
        ev_sqdiff,
    )

    grad_raw = (
        w_descr * descr_preservation_grad(c, a, b)
        + w_lap * LB_commutation_grad(c, ev_sqdiff)
        + w_comm * oplist_commutation_grad(c, op_list)
    )
    grad_locked = grad_energy_std(
        c.reshape(-1),
        w_descr,
        w_lap,
        w_comm,
        0.0,
        a,
        b,
        op_list,
        [],
        ev_sqdiff,
    ).reshape((k2, k1))

    return {
        "k1": k1,
        "k2": k2,
        "p": p,
        "n_ops": n_ops,
        "w_descr": w_descr,
        "w_lap": w_lap,
        "w_comm": w_comm,
        "c": c,
        "a": a,
        "b": b,
        "ev_sqdiff": ev_sqdiff,
        "op_list": op_list,
        "e_descr": float(e_descr),
        "e_lap": float(e_lap),
        "e_comm": float(e_comm),
        "e_total": float(e_total),
        "grad_raw": grad_raw,
        "grad_locked": grad_locked,
    }


def run_probe():
    ref = make_refine_fixture()
    obj = make_objective_fixture()

    out = {
        "status": "ok",
        "icp_nit": str(ref["icp_nit"]),
        "zoom_nit": str(ref["zoom_nit"]),
        "step": int_vec_to_csv(ref["step"]),
        "phi1_dims": dims_to_csv(ref["phi1"]),
        "phi1": mat_to_csv(ref["phi1"]),
        "phi2_dims": dims_to_csv(ref["phi2"]),
        "phi2": mat_to_csv(ref["phi2"]),
        "c0_dims": dims_to_csv(ref["c0"]),
        "c0": mat_to_csv(ref["c0"]),
        "sub_source_1b": int_vec_to_csv(ref["sub_source_0b"] + 1),
        "sub_target_1b": int_vec_to_csv(ref["sub_target_0b"] + 1),
        "icp_sub_dims": dims_to_csv(ref["icp_sub"]),
        "icp_sub": mat_to_csv(ref["icp_sub"]),
        "zoom_sub_dims": dims_to_csv(ref["zoom_sub"]),
        "zoom_sub": mat_to_csv(ref["zoom_sub"]),
        "weights": vec_to_csv(ref["weights"]),
        "phi1_w_dims": dims_to_csv(ref["phi1_w"]),
        "phi1_w": mat_to_csv(ref["phi1_w"]),
        "phi2_w_dims": dims_to_csv(ref["phi2_w"]),
        "phi2_w": mat_to_csv(ref["phi2_w"]),
        "c0_w_dims": dims_to_csv(ref["c0_w"]),
        "c0_w": mat_to_csv(ref["c0_w"]),
        "p2p_w_1b": int_vec_to_csv(ref["p2p_w_0b"] + 1),
        "fm_weighted_dims": dims_to_csv(ref["fm_weighted"]),
        "fm_weighted": mat_to_csv(ref["fm_weighted"]),
        "k1": str(obj["k1"]),
        "k2": str(obj["k2"]),
        "p": str(obj["p"]),
        "n_ops": str(obj["n_ops"]),
        "w_descr": f"{obj['w_descr']:.17g}",
        "w_lap": f"{obj['w_lap']:.17g}",
        "w_comm": f"{obj['w_comm']:.17g}",
        "A_dims": dims_to_csv(obj["a"]),
        "A": mat_to_csv(obj["a"]),
        "B_dims": dims_to_csv(obj["b"]),
        "B": mat_to_csv(obj["b"]),
        "C_dims": dims_to_csv(obj["c"]),
        "C": mat_to_csv(obj["c"]),
        "ev_sqdiff_dims": dims_to_csv(obj["ev_sqdiff"]),
        "ev_sqdiff": mat_to_csv(obj["ev_sqdiff"]),
        "e_descr": f"{obj['e_descr']:.17g}",
        "e_lap": f"{obj['e_lap']:.17g}",
        "e_comm": f"{obj['e_comm']:.17g}",
        "e_total": f"{obj['e_total']:.17g}",
        "grad_raw_dims": dims_to_csv(obj["grad_raw"]),
        "grad_raw": mat_to_csv(obj["grad_raw"]),
        "grad_locked_dims": dims_to_csv(obj["grad_locked"]),
        "grad_locked": mat_to_csv(obj["grad_locked"]),
    }

    for i, (op1, op2) in enumerate(obj["op_list"], start=1):
        out[f"op1_{i}_dims"] = dims_to_csv(op1)
        out[f"op1_{i}"] = mat_to_csv(op1)
        out[f"op2_{i}_dims"] = dims_to_csv(op2)
        out[f"op2_{i}"] = mat_to_csv(op2)

    return out


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
