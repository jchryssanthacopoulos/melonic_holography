"""Heatmap of the reduced (centre-integrated) N=2 SYK ladder kernel.

Plots, from research/tex/n2_dirac_ladders/sec_dirac.tex, eq. (reduced),

    K^(0)(x, y) = rho(x + y) V(y)
                = c * sgn(x + y) |x + y|^(-2 delta) * sgn(y) |y|^(-2 mu),

with c = v rho_0 = 2 mu (q - 1) C J b_psi^q, where
    C = 2 B(nu, 2 delta - 1) - B(nu, nu)   (eq. chain, continued Beta),
    J b_psi^q = tan(pi delta / 2) / (2 pi) (eq. props),
    delta = 1/q, nu = 1 - delta, mu = 1/2 - delta.
b_psi cancels, so c depends on q only. The same K^(0) is the charged-channel
kernel and the neutral kernel on {(f, -f')} (paragraph "The neutral channels
see the same kernel").

The kernel is singular on y = 0 (rung, |y|^-2mu) and on x = -y (rail,
|x+y|^-2delta, the reflection in the rung), so the colour scale is symmetric-log.

Usage:
    .venv/bin/python scripts/plot_reduced_kernel.py --q 4 --L 2 --n 801
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
from matplotlib.colors import SymLogNorm


def params(q):
    delta = mp.mpf(1) / q
    nu = 1 - delta
    mu = mp.mpf(1) / 2 - delta
    C = 2 * mp.beta(nu, 2 * delta - 1) - mp.beta(nu, nu)
    Jbq = mp.tan(mp.pi * delta / 2) / (2 * mp.pi)
    c = 2 * mu * (q - 1) * C * Jbq
    return float(delta), float(mu), float(C), float(c)


def kernel(x, y, delta, mu, c):
    s = x + y
    return c * np.sign(s) * np.abs(s) ** (-2 * delta) * np.sign(y) * np.abs(y) ** (-2 * mu)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--q", type=int, default=3, help="interaction order (q >= 3)")
    ap.add_argument("--L", type=float, default=2.0, help="plot window [-L, L]^2")
    ap.add_argument("--n", type=int, default=800, help="grid points per axis (even avoids the singular lines)")
    ap.add_argument("--linthresh", type=float, default=None, help="SymLogNorm linear threshold (default: |c|)")
    ap.add_argument("--out", default="results/figures", help="output directory")
    args = ap.parse_args()

    delta, mu, C, c = params(args.q)
    print(f"q={args.q}: delta={delta:.6g}, mu={mu:.6g}, C={C:.10g}, c=v*rho0={c:.10g}")

    # Cell-centred grid. With even n no y sample is 0; shifting x by a quarter
    # cell keeps x + y off zero as well (the symmetric grid would hit x = -y).
    edges = np.linspace(-args.L, args.L, args.n + 1)
    h = edges[1] - edges[0]
    xedges = edges + h / 4
    xmid = 0.5 * (xedges[:-1] + xedges[1:])
    ymid = 0.5 * (edges[:-1] + edges[1:])
    X, Y = np.meshgrid(xmid, ymid, indexing="xy")
    assert np.all(Y != 0) and np.all(X + Y != 0)
    K = kernel(X, Y, delta, mu, c)

    lin = args.linthresh if args.linthresh is not None else abs(c)
    vmax = np.percentile(np.abs(K), 99.5)
    norm = SymLogNorm(linthresh=lin, vmin=-vmax, vmax=vmax, base=10)

    plt.rcParams.update({"font.size": 12})
    fig, ax = plt.subplots(figsize=(6.4, 5.4), constrained_layout=True)
    im = ax.pcolormesh(xedges, edges, K, cmap="RdBu_r", norm=norm, shading="flat", rasterized=True)
    ax.plot([-args.L, args.L], [0, 0], "k--", lw=0.8, label=r"$y=0$: rung $|y|^{-2\mu}$")
    ax.plot([-args.L, args.L], [args.L, -args.L], "k:", lw=1.0, label=r"$x=-y$: rail $|x+y|^{-2\delta}$")
    ax.set_xlabel(r"$x=\tau_{12}$  (units of $J^{-1}$)")
    ax.set_ylabel(r"$y=\tau_{34}$  (units of $J^{-1}$)")
    ax.set_xlim(-args.L, args.L)
    ax.set_ylim(-args.L, args.L)
    ax.set_aspect("equal")
    ax.set_title(rf"$K^{{(0)}}(x,y)=\rho(x+y)\,\mathcal{{V}}(y)$,  $q={args.q}$")
    ax.legend(loc="upper right", fontsize=9, framealpha=0.85)
    cb = fig.colorbar(im, ax=ax, shrink=0.9)
    cb.set_label(rf"$K^{{(0)}}$  (symlog, linthresh $={lin:.3g}$)")

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    stem = out / f"reduced_kernel_q{args.q}"
    fig.savefig(stem.with_suffix(".png"), dpi=200)
    fig.savefig(stem.with_suffix(".pdf"))
    print(f"wrote {stem}.png and {stem}.pdf")


if __name__ == "__main__":
    main()
