"""Make regional_run/code_swfix (2026-11-08, Dustin "Yes, I stage it"): use CDOMdecay only when useGCHEM is true.
Physics-only runs (useGCHEM = F) never set CDOMdecay, so the CDOM branch removed all ocean SW heating.
New files only: code_swfix/apply_forcing.F (from code_darwin) and code_swfix/seaice_growth.F (from code_capload)."""
import os, re, sys
R = '/nobackup/dcarrol2/downscaling/nares/regional_run'; O = R + '/code_swfix'
if os.path.exists(O): sys.exit('FATAL: code_swfix exists')
os.makedirs(O)
def sub_of(src, pos):
    i = src.rfind('SUBROUTINE', 0, pos); j = src.find('\n', i); return src[i:j].strip(), src[i:pos]
# apply_forcing.F
s = open(R + '/code_darwin/apply_forcing.F').read()
a = s.index('#ifdef DARWIN_ALLOW_CDOM\n       DO j=0,sNy+1\n        DO i=0,sNx+1\n         gT_arr(i,j) = gT_arr(i,j)\n     &        - Qsw(i,j,bi,bj)*( swfracb(1)*CDOMdecay')
b = s.index('#endif /* DARWIN_ALLOW_CDOM */', a)
blk = s[a:b]
assert s.count('CDOMdecay(i,j,k,bi,bj)') == 1
cdom, plain = blk.split('#else\n')
cdom = cdom[len('#ifdef DARWIN_ALLOW_CDOM\n'):]
name, head = sub_of(s, a); assert '#include "PARAMS.h"' in head, name
new = ('#ifdef DARWIN_ALLOW_CDOM\nC     Nares 2026-11-08: CDOMdecay is set only by DARWIN_FORCING (useGCHEM). Physics-only\n'
       'C     runs must use the plain SW heating, else CDOMdecay = 0 removes all ocean SW heat.\n'
       '      IF ( useGCHEM ) THEN\n' + cdom + '      ELSE\n' + plain + '      ENDIF\n#else\n' + plain)
s2 = s[:a] + new + s[b:]
open(O + '/apply_forcing.F', 'w').write(s2)
print('apply_forcing.F: patched in', name)
# seaice_growth.F
s = open(R + '/code_capload/seaice_growth.F').read()
old = ('#ifdef DARWIN_ALLOW_CDOM\n            tmpscal2=SEAICE_SWFrac * a_QSWbyATM_open(i,j) *\n'
       '     &               CDOMdecay(i,j,2,bi,bj)\n#else\n')
assert s.count(old) == 1
name, head = sub_of(s, s.index(old)); assert '#include "PARAMS.h"' in head, name
new = ('#ifdef DARWIN_ALLOW_CDOM\nC     Nares 2026-11-08: CDOMdecay only when useGCHEM (it is never set otherwise).\n'
       '            IF ( useGCHEM ) THEN\n             tmpscal2=SEAICE_SWFrac * a_QSWbyATM_open(i,j) *\n'
       '     &                CDOMdecay(i,j,2,bi,bj)\n            ELSE\n'
       '             tmpscal2=SEAICE_SWFrac * a_QSWbyATM_open(i,j)\n            ENDIF\n#else\n')
open(O + '/seaice_growth.F', 'w').write(s.replace(old, new))
print('seaice_growth.F: patched in', name, '(cap kept:', 'heffTooHeavy = 6. _d 0' in s, ')')
