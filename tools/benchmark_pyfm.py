#!/usr/bin/env python3
import argparse
import sys
import time
from pathlib import Path


def emit(lines):
    for k, v in lines.items():
        print(f"{k}={v}")


def main():
    parser = argparse.ArgumentParser(description="Synthetic pyFM optimization benchmark")
    parser.add_argument("--n", type=int, default=160)
    parser.add_argument("--k1", type=int, default=32)
    parser.add_argument("--k2", type=int, default=32)
    parser.add_argument("--p", type=int, default=24)
    parser.add_argument("--maxit", type=int, default=200)
    parser.add_argument("--icp-nit", type=int, default=5)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(repo_root / "pyFM"))

    try:
        import numpy as np
        from scipy.optimize import fmin_l_bfgs_b
        from pyFM.optimize.base_functions import energy_func_std, grad_energy_std
        from pyFM.refine.icp import icp_refine
        from pyFM.spectral.convert import FM_to_p2p
    except Exception as exc:
        emit({"status": "error", "error": f"import_failure:{exc}"})
        return 0

    try:
        rng = np.random.default_rng(args.seed)
        t0 = time.perf_counter()

        theta = np.linspace(0.0, 2.0 * np.pi, args.n, endpoint=False)
        coords = np.column_stack([np.cos(theta), np.sin(theta)])
        shift = max(1, int(np.floor(0.15 * args.n)))
        p21 = (np.arange(args.n) + shift) % args.n

        p12 = np.zeros((args.n, args.n), dtype=float)
        p12[np.arange(args.n), p21] = 1.0

        q1, _ = np.linalg.qr(rng.normal(size=(args.n, args.k1)))
        evects1 = q1[:, : args.k1]
        evects2 = p12 @ evects1

        desc1_n = rng.normal(size=(args.n, args.p))
        desc2_n = p12 @ desc1_n + rng.normal(scale=0.01, size=(args.n, args.p))
        descr1 = evects1.T @ desc1_n
        descr2 = evects2.T @ desc2_n

        list_descr = []
        for i in range(args.p):
            # Descriptor multiplication operators in reduced bases.
            op1 = evects1.T @ (desc1_n[:, i, None] * evects1)
            op2 = evects2.T @ (desc2_n[:, i, None] * evects2)
            list_descr.append((op1, op2))

        ev1 = np.arange(1, args.k1 + 1, dtype=float)
        ev2 = np.arange(1, args.k2 + 1, dtype=float)
        ev_sqdiff = (ev1[None, :] - ev2[:, None]) ** 2
        ev_sqdiff = ev_sqdiff / ev_sqdiff.sum()

        w_descr = 1e-1
        w_lap = 1e-3
        w_dcomm = 1.0

        args_opt = (
            w_descr,
            w_lap,
            w_dcomm,
            0.0,
            descr1,
            descr2,
            list_descr,
            [],
            ev_sqdiff,
        )

        x0 = np.zeros((args.k2, args.k1), dtype=float)
        x0[0, 0] = 1.0

        t_opt0 = time.perf_counter()
        out = fmin_l_bfgs_b(
            energy_func_std,
            x0.ravel(),
            fprime=grad_energy_std,
            args=args_opt,
            maxiter=args.maxit,
        )
        t_opt = time.perf_counter() - t_opt0

        FM_12 = out[0].reshape(args.k2, args.k1)

        t_icp0 = time.perf_counter()
        FM_icp = icp_refine(
            FM_12,
            evects1,
            evects2,
            nit=args.icp_nit,
            verbose=False,
            use_adj=False,
            return_p2p=False,
            n_jobs=1,
        )
        t_icp = time.perf_counter() - t_icp0
        elapsed = time.perf_counter() - t0
        p2p = FM_to_p2p(FM_icp, evects1, evects2, use_adj=False, n_jobs=1)
        d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=2)
        accuracy = float(np.mean(p2p == p21))
        geodesic_mean = float(np.mean(d[p2p, p21]))
        scale_vals = d[np.isfinite(d) & (d > 0)]
        scale = float(np.mean(scale_vals)) if scale_vals.size > 0 else 1.0
        geodesic_norm = geodesic_mean / max(scale, 1e-12)

        emit(
            {
                "status": "ok",
                "mode": "full_pipeline_opt_plus_icp",
                "runtime_sec": f"{elapsed:.9f}",
                "runtime_opt_sec": f"{t_opt:.9f}",
                "runtime_icp_sec": f"{t_icp:.9f}",
                "objective": f"{float(out[1]):.9f}",
                "accuracy": f"{accuracy:.9f}",
                "geodesic_normalized_mean": f"{geodesic_norm:.9f}",
                "funcalls": str(out[2].get("funcalls", "")),
                "nit": str(out[2].get("nit", "")),
                "warnflag": str(out[2].get("warnflag", "")),
                "python": sys.version.split()[0],
            }
        )
        return 0
    except Exception as exc:
        emit({"status": "error", "error": f"runtime_failure:{exc}"})
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
