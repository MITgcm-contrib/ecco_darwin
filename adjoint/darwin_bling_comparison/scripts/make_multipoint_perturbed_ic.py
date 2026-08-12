#!/usr/bin/env python3
"""Write +/-eta (relative, multiplicative) perturbed copies of a Darwin
ptracer IC file at MANY surface grid points (k=0) simultaneously, for the
multi-point FD-sensitivity extension. Mirrors Kay's GCHP validation
methodology (see notes/from_kay/gchp_adjoint_validation.pdf): each target
cell's concentration is scaled by (1+eta) or (1-eta) independently, so the
absolute perturbation naturally scales with each point's own local value
(avoids negative concentrations that a fixed absolute eps would cause for
sparse tracers like NO3 across a wide dynamic range) and needs no
per-tracer eps tuning. Points are assumed well-separated (>= MIN_SEP grid
cells) so each point's perturbation stays local and doesn't contaminate
another point's cost-function neighborhood -- verified separately.
"""
import sys
import json
import numpy as np

NZ, NY, NX = 15, 64, 128


def make_pair(base_file, out_plus, out_minus, points, eta):
    d = np.fromfile(base_file, dtype='>f4').reshape(NZ, NY, NX)
    d_plus = d.copy()
    d_minus = d.copy()
    applied = []
    for j, i in points:
        base_val = float(d[0, j, i])
        d_plus[0, j, i] = base_val * (1.0 + eta)
        d_minus[0, j, i] = base_val * (1.0 - eta)
        applied.append({"j": j, "i": i, "base": base_val, "delta": eta * base_val})
    d_plus.astype('>f4').tofile(out_plus)
    d_minus.astype('>f4').tofile(out_minus)
    return applied


if __name__ == '__main__':
    base_file, out_plus, out_minus, points_json, eta_str, log_out = sys.argv[1:7]
    points = json.loads(points_json)
    eta = float(eta_str)
    applied = make_pair(base_file, out_plus, out_minus, points, eta)
    with open(log_out, 'w') as fh:
        json.dump({"eta": eta, "points": applied}, fh, indent=2)
    print(f"{base_file}: {len(points)} points perturbed by eta={eta} "
          f"-> log written to {log_out}")
