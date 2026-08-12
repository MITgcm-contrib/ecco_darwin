#!/usr/bin/env python3
"""Write +eps/-eps perturbed copies of a Darwin ptracer IC file at one
surface grid point (k=0), for the DIC/ALK/O2/NO3/PO4/FeT FD-sensitivity
pilot. Base fields are (15, 64, 128) big-endian float32, MDS convention."""
import sys
import numpy as np

NZ, NY, NX = 15, 64, 128


def make_pair(base_file, out_plus, out_minus, j, i, eps):
    d = np.fromfile(base_file, dtype='>f4').reshape(NZ, NY, NX)
    base_val = float(d[0, j, i])

    d_plus = d.copy()
    d_plus[0, j, i] = base_val + eps
    d_plus.astype('>f4').tofile(out_plus)

    d_minus = d.copy()
    d_minus[0, j, i] = base_val - eps
    d_minus.astype('>f4').tofile(out_minus)

    print(f"{base_file}: base={base_val:.6g} eps={eps:.6g} "
          f"-> plus={base_val+eps:.6g} minus={base_val-eps:.6g}")


if __name__ == '__main__':
    base_file, out_plus, out_minus, j, i, eps = sys.argv[1:7]
    make_pair(base_file, out_plus, out_minus, int(j), int(i), float(eps))
