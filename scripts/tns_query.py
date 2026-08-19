import os, json, time
import requests

# minimal tns client. bot creds come from .env (never printed). get_object(name) is the light,
# name-based call; search_coord is a cone (heavier - tns asks to use it sparingly), used only when
# a name won't resolve. >=1s between calls + backoff on 429. results cached per key so re-runs of
# the catalog cross-match don't re-hit the api.
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
from paths import TNS_CACHE as CACHE
BASE = 'https://www.wis-tns.org/api/get'


def _creds():
    env = {}
    with open(os.path.join(ROOT, '.env')) as f:
        for line in f:
            s = line.strip()
            if s and not s.startswith('#') and '=' in s:
                k, v = s.split('=', 1)
                env[k.strip()] = v.strip()
    return env['tns_id'], env['tns_bot_name'], env['tns_api_key']


_ID, _NAME, _KEY = _creds()
_HEADERS = {'user-agent': 'tns_marker{"tns_id":%s,"type":"bot","name":"%s"}' % (_ID, _NAME)}
_last = [0.0]


def _post(endpoint, data, tries=4):
    for i in range(tries):
        gap = time.time() - _last[0]
        if gap < 1.1:
            time.sleep(1.1 - gap)                 # honour the >=1s rate limit
        r = requests.post(f'{BASE}/{endpoint}', headers=_HEADERS,
                          data={'api_key': _KEY, 'data': json.dumps(data)}, timeout=60)
        _last[0] = time.time()
        if r.status_code == 429:                  # too many requests -> back off and retry
            time.sleep(2 ** i * 2)
            continue
        r.raise_for_status()
        return r.json().get('data')
    return None


def _cached(key, fn):
    os.makedirs(CACHE, exist_ok=True)
    p = os.path.join(CACHE, key.replace('/', '_') + '.json')
    if os.path.exists(p):
        return json.load(open(p))
    v = fn()
    json.dump(v, open(p, 'w'))
    return v


def tns_name(name):
    # strip the SN/AT prefix + separators -> the bare tns designation, e.g. 'SN2023IXF' -> '2023ixf'.
    import re
    m = re.match(r'\s*(?:SN|AT)[\s_-]?(\d{4}[A-Za-z]+)\s*$', str(name), re.I)
    return m.group(1)[:4] + m.group(1)[4:].lower() if m else None


def get_object(name):
    # object details by tns name (no prefix). data is the object dict directly in this api version;
    # older versions nest it under 'reply', so accept both. returns the object dict or None.
    d = _cached('obj_' + name, lambda: _post('object', {'objname': name}))
    if isinstance(d, dict) and isinstance(d.get('reply'), dict):
        return d['reply']
    return d if isinstance(d, dict) else None


def search_coord(ra, dec, radius=5.0):
    # cone search -> list of {objname, prefix, objid}. heavier; only when a name won't resolve.
    d = _cached(f'cone_{ra:.5f}_{dec:.5f}',
                lambda: _post('search', {'ra': float(ra), 'dec': float(dec),
                                         'radius': radius, 'units': 'arcsec'}))
    if isinstance(d, dict):
        return d.get('reply', [])
    return d or []
