# file_module.py
import os
import numpy as np
from config import M
from variables import DEPTH, B, U

class OutputManager:
    """
    Fast writers for hydro + tracers + rates.
    Starts fresh each run for hydro/width; tracers/rates depend on overwrite_all.
    pH is always truncated on first write.
    """
    def __init__(self, prefix: str = "",
                 depth_fmt="%.6g", u_fmt="%.6g", tracer_fmt="%.15e",
                 overwrite_all: bool = True, buffer_bytes: int = 1 << 20):
        self.prefix = prefix
        self.depth_fmt  = depth_fmt
        self.u_fmt      = u_fmt
        self.tracer_fmt = tracer_fmt
        self.buffer_bytes = buffer_bytes
        self.overwrite_all = overwrite_all
        self._open_fhs = {}  # path -> handle

        # any of these filenames (with or without .dat) will be truncated on first open
        self._force_truncate_basenames = {"pH", "pH.dat"}

        os.makedirs(prefix or ".", exist_ok=True)

        # --- Hydro files: truncate now if overwrite_all, else append
        depth_path = os.path.join(prefix, "depth.dat")
        u_path     = os.path.join(prefix, "U.dat")
        hyd_mode = "wb" if overwrite_all else "ab"
        self.depth = open(depth_path, hyd_mode, buffering=buffer_bytes)
        self.u     = open(u_path,     hyd_mode, buffering=buffer_bytes)

        # width is static → always write fresh once per run
        with open(os.path.join(prefix, "width.dat"), "wb") as f:
            f.write(b"t\t")
            np.arange(M + 1, dtype=np.int32).tofile(f, sep="\t")
            f.write(b"\n0\t")
            np.ascontiguousarray(B[:M+1], dtype=np.float64).tofile(f, sep="\t", format=self.depth_fmt)
            f.write(b"\n")

    def _canonical_name(self, filename: str) -> str:
        # ensure .dat extension
        return filename if filename.endswith(".dat") else f"{filename}.dat"

    def _get_fh(self, filename: str):
        """
        Open tracer/rate file.
        - If overwrite_all=True → truncate on first open.
        - Additionally, pH (pH / pH.dat) always truncates on first open.
        """
        fname = self._canonical_name(filename)
        path = os.path.join(self.prefix, fname)

        fh = self._open_fhs.get(path)
        if fh is None:
            os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
            base = os.path.basename(fname)
            force_trunc = base in self._force_truncate_basenames
            mode = "wb" if (self.overwrite_all or force_trunc) else "ab"
            fh = open(path, mode, buffering=self.buffer_bytes)
            self._open_fhs[path] = fh
        return fh

    def _write_row(self, fh, t: float, arr: np.ndarray, fmt: str):
        fh.write(f"{t}\t".encode())
        np.ascontiguousarray(arr, dtype=np.float64).tofile(fh, sep="\t", format=fmt)
        fh.write(b"\n")

    # ----- public API -----
    def write_hyd(self, t: float):
        self._write_row(self.depth, t, DEPTH[:M+1], self.depth_fmt)
        self._write_row(self.u,     t, U[:M+1],     self.u_fmt)

    def write_tracer(self, co: np.ndarray, filename: str, t: float):
        fh = self._get_fh(filename)
        self._write_row(fh, t, co[1:M+1], self.tracer_fmt)

    def write_rates(self, values: np.ndarray, filename: str, t: float):
        fh = self._get_fh(filename)
        self._write_row(fh, t, values[1:M+1], self.tracer_fmt)

    def close(self):
        for fh in (self.depth, self.u, *self._open_fhs.values()):
            try: fh.close()
            except: pass
        self._open_fhs.clear()
