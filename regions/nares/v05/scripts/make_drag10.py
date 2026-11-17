"""Make code_drag10 overlay copies (Nares 2026-11-15, fix (a)): replace uVel/vVel(.,.,kSrf,bi,bj) by uOcnSI/vOcnSI(.,.,bi,bj)
in the ice-ocean drag routines, add the header include, and call SEAICE_OCNVEL_AVG at the start of each SEAICE_DYNSOLVER call."""
import re, sys
SRC, OUT = sys.argv[1], sys.argv[2]
pat = re.compile(r'\b([uv])Vel\(([^,()]*),([^,()]*),kSrf,bi,bj\)')
for f in ['seaice_lsr.F', 'seaice_oceandrag_coeffs.F', 'seaice_ocean_stress.F', 'seaice_dynsolver.F']:
    s = open(f'{SRC}/{f}').read()
    s, n = pat.subn(lambda m: f'{m.group(1)}OcnSI({m.group(2)},{m.group(3)},bi,bj)', s)
    assert '#include "DYNVARS.h"' in s, f
    s = s.replace('#include "DYNVARS.h"', '#include "DYNVARS.h"\n#include "SEAICE_OCNVEL.h"', 1)
    if f == 'seaice_dynsolver.F':
        key = 'C--   interface of dynamics with atmopheric forcing fields'
        assert s.count(key) == 1
        s = s.replace(key, 'C--   Nares 2026-11-15: ocean velocity for the ice-ocean drag = top 10 m mean\n      CALL SEAICE_OCNVEL_AVG( myThid )\n\n' + key)
    for ln in s.splitlines():
        if not ln[:1] in 'cC!#' and len(ln) > 72: sys.exit(f'line > 72 in {f}: {ln}')
    open(f'{OUT}/{f}', 'w').write(s)
    print(f, 'replacements', n, '| kSrf Vel left:', len(re.findall(r'Vel\([^)]*kSrf', s)))
