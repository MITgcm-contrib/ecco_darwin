"""
File handling module (translated from file.c)
"""

from config import M
from variables import DEPTH, B, D, U

def hydwrite(t):
    """Write physical variables (depth and width) to output files."""
    with open("depth.dat", "a") as fptr1, open("width.dat", "a") as fptr2:
        fptr1.write(f"{t}\t")
        fptr2.write(f"{t}\t")

        for i in range(M + 1):
            fptr1.write(f"{DEPTH[i]:f}\t")
            fptr2.write(f"{B[i]:f}\t")

        fptr1.write("\n")
        fptr2.write("\n")

def transwrite(co, s, t):
    """Write concentration values to output files."""
    with open(s, "a") as fptr1:
        fptr1.write(f"{t}\t")

        for i in range(1, M + 1):
            fptr1.write(f"{co[i]:.15e}\t")

        fptr1.write("\n")

def Rates(co, s, t):
    """Write biogeochemical reaction rates to output files."""
    with open(s, "a") as fptr1:
        fptr1.write(f"{t}\t")
        for i in range(1, M + 1):
            fptr1.write(f"{co[i]:.15e}\t")
        fptr1.write("\n")
