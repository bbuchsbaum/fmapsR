#!/usr/bin/env python3
import argparse
import sys
import time
from pathlib import Path


def emit(lines):
    for k, v in lines.items():
        print(f"{k}={v}")


def make_problem(
    n: int,
    k: int,
    p: int,
    noise: float,
    seed: int,
    basis_noise: float,
    descriptor_corruption: float,
    eval_fraction: float,
):
    import numpy as np

    rng = np.random.default_rng(seed)

    theta = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    coords = np.column_stack([np.cos(theta), np.sin(theta)])

    shift = max(1, int(np.floor(0.15 * n)))
    p21 = (np.arange(n) + shift) % n  # target -> source, 0-based

    p12 = np.zeros((n, n), dtype=float)
    p12[np.arange(n), p21] = 1.0

    q, _ = np.linalg.qr(rng.normal(size=(n, k)))
    phi1 = q[:, :k]
    phi2 = p12 @ phi1

    if basis_noise > 0:
        q2, _ = np.linalg.qr(phi2 + basis_noise * rng.normal(size=phi2.shape))
        phi2 = q2[:, :k]

    src_desc = rng.normal(size=(n, p))
    tgt_desc = p12 @ src_desc + rng.normal(scale=noise, size=(n, p))

    if descriptor_corruption > 0:
        n_bad = max(1, int(round(descriptor_corruption * n)))
        bad_idx = rng.choice(n, size=n_bad, replace=False)
        tgt_desc[bad_idx] = rng.normal(size=(n_bad, p))

    eval_n = max(1, int(round(eval_fraction * n)))
    eval_idx = np.sort(rng.choice(n, size=eval_n, replace=False))

    return {
        "coords": coords,
        "p21": p21,
        "phi1": phi1,
        "phi2": phi2,
        "src_desc": src_desc,
        "tgt_desc": tgt_desc,
        "eval_idx": eval_idx,
    }


def run_once(
    n: int,
    k: int,
    p: int,
    noise: float,
    seed: int,
    maxit: int,
    icp_nit: int,
    basis_noise: float,
    descriptor_corruption: float,
    eval_fraction: float,
):
    import numpy as np
    from scipy.optimize import fmin_l_bfgs_b
    from pyFM.optimize.base_functions import energy_func_std, grad_energy_std
    from pyFM.refine.icp import icp_refine
    from pyFM.spectral.convert import FM_to_p2p

    prob = make_problem(
        n=n,
        k=k,
        p=p,
        noise=noise,
        seed=seed,
        basis_noise=basis_noise,
        descriptor_corruption=descriptor_corruption,
        eval_fraction=eval_fraction,
    )

    phi1 = prob["phi1"]
    phi2 = prob["phi2"]
    src_desc = prob["src_desc"]
    tgt_desc = prob["tgt_desc"]

    descr1 = phi1.T @ src_desc
    descr2 = phi2.T @ tgt_desc

    list_descr = []
    for i in range(p):
        op1 = phi1.T @ (src_desc[:, i, None] * phi1)
        op2 = phi2.T @ (tgt_desc[:, i, None] * phi2)
        list_descr.append((op1, op2))

    eval1 = np.arange(1, k + 1, dtype=float)
    eval2 = eval1 + 0.01
    ev_sqdiff = (eval1[None, :] - eval2[:, None]) ** 2
    denom = float(ev_sqdiff.sum())
    if denom > 0:
        ev_sqdiff = ev_sqdiff / denom

    args_opt = (
        1e-1,
        1e-3,
        0.5,
        0.0,
        descr1,
        descr2,
        list_descr,
        [],
        ev_sqdiff,
    )

    x0 = np.zeros((k, k), dtype=float)
    x0[0, 0] = 1.0

    t0 = time.perf_counter()
    out = fmin_l_bfgs_b(
        energy_func_std,
        x0.ravel(),
        fprime=grad_energy_std,
        args=args_opt,
        maxiter=maxit,
    )

    fm_12 = out[0].reshape(k, k)
    fm_icp = icp_refine(
        fm_12,
        phi1,
        phi2,
        nit=icp_nit,
        verbose=False,
        use_adj=False,
        return_p2p=False,
        n_jobs=1,
    )

    p2p = FM_to_p2p(fm_icp, phi1, phi2, use_adj=False, n_jobs=1)
    p21 = prob["p21"]
    eval_idx = prob["eval_idx"]

    d = np.linalg.norm(prob["coords"][:, None, :] - prob["coords"][None, :, :], axis=2)

    elapsed = time.perf_counter() - t0
    accuracy = float(np.mean(p2p[eval_idx] == p21[eval_idx]))
    geodesic_mean = float(np.mean(d[p2p[eval_idx], p21[eval_idx]]))

    scale_vals = d[np.isfinite(d) & (d > 0)]
    scale = float(np.mean(scale_vals)) if scale_vals.size > 0 else 1.0
    geodesic_norm = geodesic_mean / max(scale, 1e-12)

    gram = fm_icp.T @ fm_icp
    ident = np.eye(gram.shape[0], dtype=gram.dtype)

    return {
        "runtime_sec": elapsed,
        "objective": float(out[1]),
        "map_fro_norm": float(np.linalg.norm(fm_icp, ord="fro")),
        "map_orth_resid": float(np.linalg.norm(gram - ident, ord="fro")),
        "accuracy": accuracy,
        "geodesic_mean": geodesic_mean,
        "geodesic_normalized_mean": geodesic_norm,
        "eval_n": int(eval_idx.size),
    }


def main():
    parser = argparse.ArgumentParser(description="Synthetic pyFM parity benchmark with quality metrics")
    parser.add_argument("--n", type=int, default=120)
    parser.add_argument("--k", type=int, default=24)
    parser.add_argument("--p", type=int, default=18)
    parser.add_argument("--noise", type=float, default=0.02)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--maxit", type=int, default=120)
    parser.add_argument("--icp-nit", type=int, default=3)
    parser.add_argument("--basis-noise", type=float, default=0.0)
    parser.add_argument("--descriptor-corruption", type=float, default=0.0)
    parser.add_argument("--eval-fraction", type=float, default=1.0)
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(repo_root / "pyFM"))

    try:
        import numpy as np  # noqa: F401
        import scipy  # noqa: F401
        import pyFM  # noqa: F401
    except Exception as exc:
        emit({"status": "error", "error": f"import_failure:{exc}"})
        return 0

    try:
        out = run_once(
            n=args.n,
            k=args.k,
            p=args.p,
            noise=args.noise,
            seed=args.seed,
            maxit=args.maxit,
            icp_nit=args.icp_nit,
            basis_noise=args.basis_noise,
            descriptor_corruption=args.descriptor_corruption,
            eval_fraction=args.eval_fraction,
        )

        emit(
            {
                "status": "ok",
                "runtime_sec": f"{out['runtime_sec']:.9f}",
                "objective": f"{out['objective']:.9f}",
                "map_fro_norm": f"{out['map_fro_norm']:.9f}",
                "map_orth_resid": f"{out['map_orth_resid']:.9f}",
                "accuracy": f"{out['accuracy']:.9f}",
                "geodesic_mean": f"{out['geodesic_mean']:.9f}",
                "geodesic_normalized_mean": f"{out['geodesic_normalized_mean']:.9f}",
                "eval_n": str(out["eval_n"]),
                "python": sys.version.split()[0],
            }
        )
        return 0
    except Exception as exc:
        emit({"status": "error", "error": f"runtime_failure:{exc}"})
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
