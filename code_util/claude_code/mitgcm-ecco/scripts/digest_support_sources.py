#!/usr/bin/env python3
"""Turn the raw cache from fetch_support_sources.sh into topic-bucketed digests for distillation.

Usage:
  digest_support_sources.py [cache_dir] [--since YEAR]

Reads <cache>/mail/*.txt.gz (pipermail mbox) and <cache>/github/*.jsonl, and writes
<cache>/digest/<topic>.md: answered threads only (a reply from someone other than the asker),
quoted text and signatures stripped, messages truncated, threads with a core-developer reply
first, newest first within that. <cache>/digest/INDEX.md lists counts and sizes.

The digests are raw material: references/troubleshooting.md is written from them by hand/agents
(symptom -> cause -> fix, with links). Topic matching is keyword-based and approximate; a thread
can land in several buckets.
"""
import argparse, collections, email.utils, gzip, glob, json, os, re

ap = argparse.ArgumentParser()
ap.add_argument('cache', nargs='?', default=os.path.expanduser('~/.cache/mitgcm-ecco-sources'))
ap.add_argument('--since', type=int, default=2003)
A = ap.parse_args()
OUT = os.path.join(A.cache, 'digest'); os.makedirs(OUT, exist_ok=True)
LIST = 'http://mailman.mitgcm.org/pipermail/mitgcm-support'

# people whose answers carry weight (mail display names / GitHub logins, lower case substrings)
CORE = ['campin', 'jmc', 'losch', 'mjlosch', 'jahn', 'heimbach', 'forget', 'gaelforget', 'doddridge',
        'edoddridge', 'chris hill', 'christophernhill', 'adcroft', 'menemenlis', 'dimitris', 'molod',
        'jeff scott', 'jscott', 'timothy smith', 'timothyas', 'tim smith', 'fenty', 'ian fenty', 'ifenty',
        'an nguyen', 'atnguyen', 'ou wang', 'ouwang', 'dutkiewicz', 'stephanie', 'jm campin',
        'jean-michel', 'dfer', 'ryan abernathey', 'rabernat', 'mpatrick', 'mattpatrick', 'oliver jahn']

TOPICS = collections.OrderedDict([
  ('build',       r'genmake|optfile|compil|gfortran|ifort|intel fortran|linker|undefined reference|makefile|make depend|cpp_eeoptions|build_options|\.f:\d|segmentation fault at compile|mpif90|mpi_?include|netcdf.*lib|lnetcdf'),
  ('parallel',    r'\bmpi\b|mpirun|nPx|nPy|tile|exch2|exch |blank tile|sNx|OLx|overlap|ranks?|processor|openmp|multithread|scaling|hang'),
  ('instability', r'blow ?up|blew up|nan\b|\bnans\b|unstable|instabilit|cfl|courant|deltaT|time ?step|crash|diverg|grows?|explod|stop.*abnormal|abnormal end|monitor.*huge|exceeds|too large'),
  ('grid',        r'cubed.?sphere|\bllc\b|lat.?lon.?cap|curvilinear|horizgridfile|delX|delY|delR|bathymetry|topography|hfac|partial cell|vertical grid|z\*|rstar|sigma|r_low|depth file|grid file|tile00'),
  ('obcs',        r'\bobcs\b|open boundar|obeast|obwest|obnorth|obsouth|ob_[iyjx]|orlanski|sponge|boundary condition'),
  ('seaice',      r'sea.?ice|seaice|\bheff\b|\barea\b.*ice|thsice|ice dynamic|lsr|jfnk|evp\b|snow'),
  ('forcing',     r'\bexf\b|forcing|atmospheric|bulk formula|runoff|surface flux|qnet|empmr|fu\b|fv\b|wind stress|interpolat|cal_|calendar|startdate|period|cheapaml|relax|restoring|\bsss\b|\bsst\b'),
  ('diagnostics', r'diagnostic|data\.diagnostics|frequency|levels|timephase|fields\(|available_diagnostics|statistic|\bmnc\b|netcdf output|output file|\.meta|\.data\b|rdmds|budget'),
  ('restart',     r'pickup|restart|checkpoint|niter0|ntimesteps|pchkptfreq|chkptfreq|continue run'),
  ('adjoint',     r'adjoint|\btaf\b|tapenade|autodiff|\bad\b|tangent linear|\bctrl\b|\bcost\b|gradient|tamc|tape|store directive|optim|ecco|smooth|\bprofiles\b|grdchk|divided adjoint|openad'),
  ('tracers',     r'ptracer|gchem|\bdic\b|darwin|biogeochem|carbon|nutrient|tracer|gmredi|redi|bolus|kpp|ggl90|mixing|diffusiv|viscos|leith|smagorinsky|vertical mixing|convect'),
  ('physics',     r'nonhydrostatic|non-hydrostatic|free surface|implicit|eos|equation of state|teos|jmd95|mdjwf|linear eos|salt|freshwater|shelfice|ice shelf|rbcs|bottom drag|no.?slip|coriolis|beta plane|tides?|momentum|vorticity|mom_vecinv|advection scheme|tempadvscheme|cg2d|cg3d|solver|converge'),
])
TOPIC_RE = {k: re.compile(v, re.I) for k, v in TOPICS.items()}

def is_core(who): w = who.lower(); return any(c in w for c in CORE)

def clean(body, limit):
    keep, lines = [], body.splitlines()
    for ln in lines:
        s = ln.strip()
        if s.startswith('>'): continue
        if re.match(r'^(-------------- next part|An HTML attachment was scrubbed|_{10,}|-- ?$|Sent from my)', s): break
        if re.match(r'^On .{0,120}(wrote|écrit|schrieb):?\s*$', s) or re.match(r'^-{2,} ?Original Message', s): break
        if re.match(r'^(From|Sent|To|Cc|Subject|Date):', s) and keep and not keep[-1].strip(): break
        if re.match(r'^(URL|Name|Type|Size|Desc): ', s): continue
        keep.append(ln.rstrip())
    t = re.sub(r'\n{3,}', '\n\n', '\n'.join(keep)).strip()
    return t if len(t) <= limit else t[:limit] + ' [...]'

def norm_subj(s):
    s = re.sub(r'\[MITgcm-support\]', '', s, flags=re.I)
    s = re.sub(r'^(\s*(re|aw|fwd?|sv|antw|r)\s*(\[\d+\])?\s*:\s*)+', '', s, flags=re.I)
    return re.sub(r'\s+', ' ', s).strip().lower()

# ---------- mailing list ----------
msgs = []
for path in sorted(glob.glob(os.path.join(A.cache, 'mail', '*.txt.gz'))):
    month = os.path.basename(path)[:-7]
    if int(month[:4]) < A.since: continue
    txt = gzip.open(path, 'rt', errors='replace').read()
    for raw in re.split(r'\n(?=From \S+ at \S+\s+\w{3} \w{3}\s+\d+ [\d:]+ \d{4}\n)', '\n' + txt):
        if not raw.strip(): continue
        head, _, body = raw.strip('\n').partition('\n\n')
        h = {}
        for m in re.finditer(r'^(From|Date|Subject|Message-ID|In-Reply-To):[ \t]*(.*(?:\n[ \t]+.*)*)', head, re.M):
            h[m.group(1)] = re.sub(r'\s+', ' ', m.group(2)).strip()
        if 'Subject' not in h: continue
        who = re.sub(r'^.*\((.*)\)\s*$', r'\1', h.get('From', '?'))
        try: dt = email.utils.parsedate_to_datetime(h.get('Date', '')).strftime('%Y-%m-%d')
        except Exception: dt = month
        msgs.append(dict(who=who, date=dt, month=month, subj=h['Subject'], key=norm_subj(h['Subject']),
                         mid=h.get('Message-ID', ''), irt=h.get('In-Reply-To', ''), body=body))

# thread by In-Reply-To where possible, else by normalised subject within +-3 months of the root
byid, threads = {}, collections.OrderedDict()
for m in msgs:
    root = byid.get(m['irt'])
    if root is None:
        cand = threads.get(m['key'])
        root = cand if cand and m['key'] and abs(int(m['month'][:4]) - int(cand[0]['month'][:4])) <= 1 else None
    if root is None:
        root = []; threads[m['key'] if m['key'] not in threads else m['key'] + '#' + m['mid']] = root
    root.append(m); byid[m['mid']] = root

# ---------- GitHub ----------
gh = []
gdir = os.path.join(A.cache, 'github')
if os.path.exists(os.path.join(gdir, 'issues.jsonl')):
    com = collections.defaultdict(list)
    for fn in ('comments.jsonl', 'review_comments.jsonl'):
        p = os.path.join(gdir, fn)
        if os.path.exists(p):
            for ln in open(p):
                c = json.loads(ln); com[c['issue']].append(c)
    for ln in open(os.path.join(gdir, 'issues.jsonl')):
        i = json.loads(ln)
        if int(i['created_at'][:4]) < A.since: continue
        cs = sorted(com.get(i['number'], []), key=lambda c: c['created_at'])
        # PRs: keep fixes / bugs only (they describe a defect and its cure); issues: keep all with discussion
        if i['is_pr'] and not re.search(r'fix|bug|wrong|error|incorrect|crash|nan|broken|issue|problem', (i['title'] or '') + ' ' + ' '.join(i['labels']), re.I):
            continue
        if not i['is_pr'] and not cs: continue
        kind = 'PR' if i['is_pr'] else 'issue'
        posts = [dict(who=i['user'], date=i['created_at'][:10], body=i['body'] or '')] + \
                [dict(who=c['user'], date=c['created_at'][:10], body=(f"[{c['path']}] " if c.get('path') else '') + (c['body'] or '')) for c in cs]
        gh.append(dict(src='github', title=f"{kind} #{i['number']}: {i['title']} [{i['state']}]",
                       url=f"https://github.com/MITgcm/MITgcm/{'pull' if i['is_pr'] else 'issues'}/{i['number']}",
                       posts=posts, labels=i['labels']))

# ---------- assemble ----------
items = []
for t in threads.values():
    askers = {t[0]['who'].lower()}
    replies = [m for m in t[1:] if m['who'].lower() not in askers]
    if not replies: continue
    m0 = t[0]
    url = f"{LIST}/{m0['month']}/thread.html"
    posts = [dict(who=m['who'], date=m['date'], body=m['body']) for m in t]
    items.append(dict(src='mail', title=re.sub(r'\[MITgcm-support\]\s*', '', m0['subj']), url=url, posts=posts, labels=[]))
items += gh

def render(it):
    core = any(is_core(p['who']) for p in it['posts'][1:])
    out = [f"### {it['title']}", f"{it['url']} | {it['posts'][0]['date']} | {len(it['posts'])} posts{' | CORE' if core else ''}"]
    for k, p in enumerate(it['posts'][:12]):
        b = clean(p['body'], 2500 if k == 0 else 1800)
        if b: out.append(f"-- {p['who']}{' (core)' if is_core(p['who']) else ''}, {p['date']}:\n{b}")
    return core, '\n'.join(out) + '\n'

buckets = collections.defaultdict(list)
for it in items:
    core, txt = render(it)
    text = it['title'] + ' ' + ' '.join(p['body'][:3000] for p in it['posts'][:3])
    hits = [k for k, r in TOPIC_RE.items() if r.search(text)]
    # bucket by title first (precise), else by the strongest body match (most keyword hits)
    th = [k for k, r in TOPIC_RE.items() if r.search(it['title'])]
    if th: hits = th[:2]
    elif hits: hits = sorted(hits, key=lambda k: -len(TOPIC_RE[k].findall(text)))[:2]
    else: hits = ['other']
    for k in hits: buckets[k].append((core, it['posts'][0]['date'], txt))

idx = ['# Digest index', '', f'cache: {A.cache}', f'mail messages: {len(msgs)}, answered threads/issues: {len(items)}', '',
       '| topic | items | core-answered | chars |', '|---|---|---|---|']
for k in list(TOPICS) + ['other']:
    v = sorted(buckets.get(k, []), key=lambda x: (not x[0], x[1]), reverse=False)
    v = sorted(v, key=lambda x: x[1], reverse=True); v = sorted(v, key=lambda x: not x[0])
    body = '\n'.join(x[2] for x in v)
    with open(os.path.join(OUT, f'{k}.md'), 'w') as f: f.write(f'# {k}\n\n' + body)
    idx.append(f'| {k} | {len(v)} | {sum(x[0] for x in v)} | {len(body)} |')
open(os.path.join(OUT, 'INDEX.md'), 'w').write('\n'.join(idx) + '\n')
print('\n'.join(idx))
