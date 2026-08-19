import os, re, json, time, argparse, warnings
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.simbad import Simbad
from astroquery.ipac.ned import Ned
from astroquery.ipac.irsa.irsa_dust import IrsaDust
warnings.filterwarnings('ignore')

# clean the raw uv sn catalog before the full run:
# dedup near-duplicate rows, coord cross-match against simbad+ned (canonical name, type, z, host),
# grab galactic e(b-v) per target, drop non-uv rows and the two extended snrs.
# runs on native windows python (astroquery), NOT wsl. resumable via a per-target json cache.

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
from paths import CATALOG_RAW as RAW, CATALOG as OUT, XMATCH_CACHE as CACHE
HOST_EBV = os.path.join(ROOT, 'reference', 'host_ebv.csv')   # curated host reddening (see file header)

DROP = {'N132D-KNOT', 'E0102-HOTSPOT'}   # extended snrs, point-source box does not fit them

# simbad: ask for object type + redshift up front
sb = Simbad()
sb.TIMEOUT = 60
try:
    sb.add_votable_fields('otype', 'rvz_redshift')   # newer astroquery field names
except Exception:
    sb.add_votable_fields('otype', 'z_value')


def _val(tab, col):
    if tab is None or col not in tab.colnames or len(tab) == 0:
        return None
    v = tab[col][0]
    try:
        if hasattr(v, 'mask') and v.mask:
            return None
    except Exception:
        pass
    s = str(v).strip()
    return None if s in ('', '--', 'nan', 'None') else v


def name_variants(name):
    # proposer names are messy (SN-2005IP, PTF11KLY, M82-SN). try a few normalizations.
    n = name.strip()
    vs = [n, n.replace('-', ''), n.replace('-', ' ')]
    m = re.match(r'(SN)[-_ ]?(\d{4})\s*([A-Za-z]*)', n, re.I)
    if m:
        yr, tail = m.group(2), m.group(3).lower()
        vs += [f'SN {yr}{tail}', f'SN{yr}{tail}']
    return list(dict.fromkeys(v for v in vs if v))


def simbad_by_name(name):
    for v in name_variants(name):
        try:
            t = sb.query_object(v)
        except Exception:
            t = None
        if t is not None and len(t):
            return (str(_val(t, 'main_id') or '').strip() or None,
                    str(_val(t, 'otype') or '').strip() or None,
                    _val(t, 'rvz_redshift'))
    return None, None, None


def ned_z_by_name(name):
    for v in name_variants(name):
        try:
            t = Ned.query_object(v)
        except Exception:
            continue
        z = _val(t, 'Redshift')
        if z is not None and np.isfinite(float(z)):
            return float(z)
    return None


def coord_fallback(ra, dec):
    # name did not resolve. cone search, prefer an SN row for the name and any galaxy z for rest-frame.
    coord = SkyCoord(ra, dec, unit='deg')
    canonical = otype = z = host = None
    try:
        t = sb.query_region(coord, radius=8 * u.arcsec)
        if t is not None and len(t):
            sep = coord.separation(SkyCoord(t['ra'], t['dec'], unit='deg'))
            order = np.argsort(sep)
            sn = [i for i in order if str(t['otype'][i]).startswith('SN')]
            pick = sn[0] if sn else order[0]
            canonical = str(t['main_id'][pick]).strip()
            otype = str(t['otype'][pick]).strip()
            zc = _val(t[[pick]], 'rvz_redshift')
            if zc is not None:
                z = float(zc)
    except Exception:
        pass
    if z is None:
        try:
            n = Ned.query_region(coord, radius=8 * u.arcsec)
            if n is not None and len(n):
                gal = [r for r in n if 'G' in str(r['Type']) and r['Redshift'] is not None
                       and np.isfinite(float(r['Redshift']))]
                if gal:
                    z = float(gal[0]['Redshift']); host = str(gal[0]['Object Name']).strip()
        except Exception:
            pass
    return canonical, otype, z, host


def xmatch_one(name, ra, dec):
    out = {'canonical': None, 'otype': None, 'z': None, 'z_source': None, 'host': None, 'flags': []}
    canonical, otype, z_sb = simbad_by_name(name)
    if canonical is None:
        out['flags'].append('name_unresolved')
        canonical, otype, z_c, host = coord_fallback(ra, dec)
        out['host'] = host
        z_sb = z_c
    out['canonical'] = canonical
    out['otype'] = otype
    # z priority: ned-by-name (well populated for SNe) -> simbad/coord -> 0
    z_ned = ned_z_by_name(name)
    if z_ned is not None:
        out['z'] = z_ned; out['z_source'] = 'ned'
    elif z_sb is not None:
        out['z'] = float(z_sb); out['z_source'] = 'simbad'
    else:
        out['z'] = 0.0; out['z_source'] = 'none'; out['flags'].append('no_z')
    return out


def ebv_one(ra, dec):
    coord = SkyCoord(ra, dec, unit='deg')
    try:
        t = IrsaDust.get_query_table(coord, section='ebv')
        # sfd mean e(b-v)
        for c in t.colnames:
            if 'SFD' in c and 'mean' in c.lower():
                return float(t[c][0])
        return float(t['ext SFD mean'][0])
    except Exception:
        return np.nan


def tns_one(name, canonical, ra, dec):
    # tns cross-match: canonical name + type + z + DISCOVERY DATE (for the day-phase labels, the
    # gap for every non-2023ixf SN). query by tns name first (light), cone-search fallback (heavier)
    # only when a name won't resolve. tns_query rate-limits + caches internally.
    import tns_query as tq
    from astropy.time import Time
    out = {'tns_name': None, 'tns_type': None, 'tns_z': None,
           'tns_disc_date': None, 'tns_disc_mjd': None}
    obj = None
    for nm in (name, canonical):
        tn = tq.tns_name(nm) if nm else None
        if tn:
            try:
                obj = tq.get_object(tn)
            except Exception:
                obj = None
            if obj:
                break
    if obj is None:
        try:
            res = tq.search_coord(ra, dec, radius=5.0)
            if res:
                obj = tq.get_object(res[0]['objname'])
        except Exception:
            obj = None
    if isinstance(obj, dict):
        out['tns_name'] = ((obj.get('name_prefix') or '') + ' ' + (obj.get('objname') or '')).strip() or None
        out['tns_type'] = (obj.get('object_type') or {}).get('name')
        z = obj.get('redshift')
        out['tns_z'] = float(z) if z not in (None, '', 'null') else None
        dd = obj.get('discoverydate')
        out['tns_disc_date'] = dd
        if dd:
            try:
                out['tns_disc_mjd'] = round(float(Time(str(dd).replace(' ', 'T')).mjd), 3)
            except Exception:
                pass
    return out


def cached(key, fn):
    os.makedirs(CACHE, exist_ok=True)
    p = os.path.join(CACHE, key + '.json')
    if os.path.exists(p):
        with open(p) as f:
            return json.load(f)
    v = fn()
    with open(p, 'w') as f:
        json.dump(v, f)
    return v


def join_host_ebv(final):
    # join the curated host E(B-V) (catalog ebv is MW-only) from the authoritative reference file. shared by
    # the full build and the network-free --host-sync path so the catalog mirror never drifts from the source.
    final['host_ebv'] = 0.0; final['host_ebv_err'] = np.nan; final['host_ebv_src'] = 'none (MW-only)'
    if os.path.exists(HOST_EBV):
        he = pd.read_csv(HOST_EBV, comment='#')
        hmap = {str(n).upper(): row for n, row in zip(he['name'], he.to_dict('records'))}
        for i, nm in final['name'].items():
            h = hmap.get(str(nm).upper())
            if h is not None:
                final.at[i, 'host_ebv'] = float(h['host_ebv'])
                final.at[i, 'host_ebv_err'] = float(h['host_ebv_err']) if pd.notna(h['host_ebv_err']) else np.nan
                final.at[i, 'host_ebv_src'] = h['host_ebv_src']
        print(f"host E(B-V): {int((final['host_ebv_src'] != 'none (MW-only)').sum())} SNe from {HOST_EBV}")
    return final


def host_sync():
    # re-join host_ebv onto the EXISTING catalog (no network). run this after any reference/host_ebv.csv edit;
    # the network-heavy full rebuild is not needed just to propagate a host reddening.
    final = pd.read_csv(OUT)
    join_host_ebv(final)
    final.to_csv(OUT, index=False)
    print(f'host-sync: {OUT} re-joined from {HOST_EBV}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--limit', type=int, default=None)   # smoke test on the first N rows
    ap.add_argument('--host-sync', action='store_true',  # network-free: just re-join host_ebv onto the catalog
                    help='re-join reference/host_ebv.csv onto the existing catalog and exit (no network queries)')
    a = ap.parse_args()
    if a.host_sync:
        host_sync(); return

    df = pd.read_csv(RAW)
    if a.limit:
        df = df.head(a.limit)
    df = df[~df['name'].isin(DROP)].copy()          # drop the two snrs
    df = df[df['has_uv'] == True].copy()            # drop rows with no uv grating
    print(f'{len(df)} rows after dropping snrs + non-uv')

    rows = []
    for i, r in df.iterrows():
        ra, dec = float(r['ra']), float(r['dec'])
        key = f"{r['name']}_{ra:.5f}_{dec:.5f}".replace('/', '_')
        xm = cached('xm_' + key, lambda: xmatch_one(str(r['name']), ra, dec))
        ebv = cached('ebv_' + key, lambda: {'ebv': ebv_one(ra, dec)})['ebv']
        tns = cached('tns_' + key, lambda: tns_one(str(r['name']), xm.get('canonical'), ra, dec))
        # tns redshift as the last resort before defaulting to 0: recent + survey-named SNe (ZTF, AT,
        # PSN, internal names) mostly don't resolve in ned/simbad but DO have a tns z we already fetched
        if xm.get('z_source') == 'none' and tns.get('tns_z') not in (None, '', 0, 0.0):
            try:
                ztns = float(tns['tns_z'])
                if np.isfinite(ztns) and ztns > 0:
                    xm['z'] = ztns
                    xm['z_source'] = 'tns'
                    xm['flags'] = [fl for fl in xm.get('flags', []) if fl != 'no_z']
            except Exception:
                pass
        rows.append({**r.to_dict(), **xm, **tns, 'ebv': ebv,
                     'flags': ','.join(xm['flags']) or 'none'})
        print(f"{r['name']:20s} -> {str(xm['canonical']):22s} z={xm['z']:.5f}({xm['z_source']}) "
              f"ebv={ebv:.3f} {xm['otype']}  tns={tns['tns_name']} disc={tns['tns_disc_mjd']}")
        time.sleep(0.1)

    clean = pd.DataFrame(rows)

    # dedup: group rows that resolved to the same canonical simbad id (merge n_spec + gratings).
    # rows with no canonical fall back to a tight coord grouping.
    clean['group'] = clean['canonical'].fillna('')
    unresolved = clean['group'] == ''
    coords = SkyCoord(clean.loc[unresolved, 'ra'].values, clean.loc[unresolved, 'dec'].values, unit='deg')
    tags = []
    used = []
    for c in coords:
        m = [j for j, u0 in enumerate(used) if c.separation(u0) < 3 * u.arcsec]
        if m:
            tags.append(f'coord{m[0]}')
        else:
            used.append(c); tags.append(f'coord{len(used)-1}')
    clean.loc[unresolved, 'group'] = tags

    merged = []
    for g, sub in clean.groupby('group'):
        best = sub.loc[sub['n_spec'].idxmax()].to_dict()
        grat = sorted(set(';'.join(sub['gratings'].dropna().astype(str)).split(';')) - {''})
        instr = sorted(set(';'.join(sub['instr'].dropna().astype(str)).split(';')) - {''})
        best['gratings'] = ';'.join(grat)
        best['instr'] = ';'.join(instr)
        best['n_spec'] = int(sub['n_spec'].sum())
        best['n_merged'] = len(sub)
        merged.append(best)
    final = pd.DataFrame(merged).drop(columns=['group'])
    final = final.sort_values('n_spec', ascending=False).reset_index(drop=True)

    # join curated host E(B-V) from the authoritative reference (catalog ebv is MW-only).
    join_host_ebv(final)

    cols = ['name', 'canonical', 'ra', 'dec', 'z', 'z_source', 'ebv', 'host_ebv', 'host_ebv_err',
            'host_ebv_src', 'otype', 'host',
            'instr', 'gratings', 'n_spec', 'n_merged', 'is_remnant', 'has_uv', 'classification',
            'tns_name', 'tns_type', 'tns_z', 'tns_disc_date', 'tns_disc_mjd', 'flags']
    final = final[[c for c in cols if c in final.columns]]
    final.to_csv(OUT, index=False)
    print(f'\nwrote {OUT}: {len(final)} unique SNe (from {len(clean)} rows)')
    print('no_z:', int((final['z_source'] == 'none').sum()),
          ' | ned z:', int((final['z_source'] == 'ned').sum()),
          ' | simbad z:', int((final['z_source'] == 'simbad').sum()))


if __name__ == '__main__':
    main()
