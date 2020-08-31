from __future__ import (
    absolute_import,
    division,
    print_function,
)

import os
import re
from os.path import (
    basename,
    getsize,
    isfile,
    join,
)

import yaml

from planemo.conda_verify.const import (
    FIELDS,
    LICENSE_FAMILIES,
)
from planemo.conda_verify.utils import (
    all_ascii,
    get_bad_seq,
    memoized,
)

PEDANTIC = True
sel_pat = re.compile(r'(.+?)\s*\[(.+)\]$')
name_pat = re.compile(r'[a-z0-9_][a-z0-9_\-\.]*$')
version_pat = re.compile(r'[\w\.]+$')
url_pat = re.compile(r'(ftp|http(s)?)://')


class RecipeError(Exception):
    pass


def ns_cfg(cfg):
    plat = cfg['plat']
    py = cfg['PY']
    np = cfg['NPY']
    for x in py, np:
        assert isinstance(x, int), x
    return dict(
        nomkl=False,
        debug=False,
        linux=plat.startswith('linux-'),
        linux32=bool(plat == 'linux-32'),
        linux64=bool(plat == 'linux-64'),
        armv7l=False,
        arm=False,
        ppc64le=False,
        osx=plat.startswith('osx-'),
        unix=plat.startswith(('linux-', 'osx-')),
        win=plat.startswith('win-'),
        win32=bool(plat == 'win-32'),
        win64=bool(plat == 'win-64'),
        x86=plat.endswith(('-32', '-64')),
        x86_64=plat.endswith('-64'),
        py=py,
        py3k=bool(30 <= py < 40),
        py2k=bool(20 <= py < 30),
        py26=bool(py == 26),
        py27=bool(py == 27),
        py33=bool(py == 33),
        py34=bool(py == 34),
        py35=bool(py == 35),
        np=np,
    )


def select_lines(data, namespace):
    lines = []
    for line in data.splitlines():
        line = line.rstrip()
        m = sel_pat.match(line)
        if m:
            if PEDANTIC:
                x = m.group(1).strip()
                # error on comment, unless the whole line is a comment
                if '#' in x and not x.startswith('#'):
                    raise RecipeError("found commented selector: %s" % line)
            cond = m.group(2)
            if eval(cond, namespace, {}):
                lines.append(m.group(1))
            continue
        lines.append(line)
    return '\n'.join(lines) + '\n'


@memoized
def yamlize(data):
    res = yaml.safe_load(data)
    # ensure the result is a dict
    if res is None:
        res = {}
    return res


def parse(data, cfg):
    if cfg is not None:
        data = select_lines(data, ns_cfg(cfg))
    # ensure we create new object, because yamlize is memoized
    return dict(yamlize(data))


def get_field(meta, field, default=None):
    section, key = field.split('/')
    submeta = meta.get(section)
    if submeta is None:
        submeta = {}
    res = submeta.get(key)
    if res is None:
        res = default
    return res


def check_name(name):
    if name:
        name = str(name)
    else:
        raise RecipeError("package name missing")
    if not name_pat.match(name) or name.endswith(('.', '-', '_')):
        raise RecipeError("invalid package name '%s'" % name)
    seq = get_bad_seq(name)
    if seq:
        raise RecipeError("'%s' is not allowed in "
                          "package name: '%s'" % (seq, name))


def check_version(ver):
    if ver:
        ver = str(ver)
    else:
        raise RecipeError("package version missing")
    if not version_pat.match(ver):
        raise RecipeError("invalid version '%s'" % ver)
    if ver.startswith(('_', '.')) or ver.endswith(('_', '.')):
        raise RecipeError("version cannot start or end with '_' or '.': %s" %
                          ver)
    seq = get_bad_seq(ver)
    if seq:
        raise RecipeError("'%s' not allowed in version '%s'" % (seq, ver))


def check_build_number(bn):
    if not (isinstance(bn, int) and bn >= 0):
        raise RecipeError("build/number '%s' (not a positive interger)" % bn)


def check_requirements(meta):
    for req in get_field(meta, 'requirements/run', []):
        name = req.split()[0]
        if not name_pat.match(name):
            raise RecipeError("invalid run requirement name '%s'" % name)


def check_license_family(meta):
    if not PEDANTIC:
        return
    lf = get_field(meta, 'about/license_family',
                   get_field(meta, 'about/license'))
    if lf not in LICENSE_FAMILIES:
        print("""\
Error: license_family is invalid: %s
Note that about/license_family falls back to about/license.
Allowed license families are:""" % lf)
        for x in LICENSE_FAMILIES:
            print("  - %s" % x)
        raise RecipeError("wrong license family")


def check_url(url):
    if not url_pat.match(url):
        raise RecipeError("not a valid URL: %s" % url)


def check_about(meta):
    summary = get_field(meta, 'about/summary')
    if summary and len(summary) > 80:
        msg = "summary exceeds 80 characters"
        if PEDANTIC:
            raise RecipeError(msg)
        else:
            print("Warning: %s" % msg)

    for field in ('about/home', 'about/dev_url', 'about/doc_url',
                  'about/license_url'):
        url = get_field(meta, field)
        if url:
            check_url(url)

    check_license_family(meta)


hash_pat = {'md5': re.compile(r'[a-f0-9]{32}$'),
            'sha1': re.compile(r'[a-f0-9]{40}$'),
            'sha256': re.compile(r'[a-f0-9]{64}$')}


def check_source(meta):
    src = meta.get('source')
    if not src:
        return
    fn = src.get('fn')
    if fn:
        for ht in 'md5', 'sha1', 'sha256':
            hexgigest = src.get(ht)
            if hexgigest and not hash_pat[ht].match(hexgigest):
                raise RecipeError("invalid hash: %s" % hexgigest)
        url = src.get('url')
        if url:
            check_url(url)

    git_url = src.get('git_url')
    if git_url and (src.get('git_tag') and src.get('git_branch')):
        raise RecipeError("cannot specify both git_branch and git_tag")


def validate_meta(meta):
    for section in meta:
        if PEDANTIC and section not in FIELDS:
            raise RecipeError("Unknown section: %s" % section)
        submeta = meta.get(section)
        if submeta is None:
            submeta = {}
        for key in submeta:
            if PEDANTIC and key not in FIELDS[section]:
                raise RecipeError("in section %r: unknown key %r" %
                                  (section, key))

    check_name(get_field(meta, 'package/name'))
    check_version(get_field(meta, 'package/version'))
    check_build_number(get_field(meta, 'build/number', 0))
    check_requirements(meta)
    check_about(meta)
    check_source(meta)


def validate_files(recipe_dir, meta):
    for field in 'test/files', 'source/patches':
        flst = get_field(meta, field)
        if not flst:
            continue
        for fn in flst:
            if PEDANTIC and fn.startswith('..'):
                raise RecipeError("path outsite recipe: %s" % fn)
            path = join(recipe_dir, fn)
            if isfile(path):
                continue
            raise RecipeError("no such file '%s'" % path)


def iter_cfgs():
    for py in 27, 34, 35:
        for plat in 'linux-64', 'linux-32', 'osx-64', 'win-32', 'win-64':
            yield dict(plat=plat, PY=py, NPY=111)


def dir_size(dir_path):
    return sum(sum(getsize(join(root, fn)) for fn in files)
               for root, unused_dirs, files in os.walk(dir_path))


def check_dir_content(recipe_dir):
    disallowed_extensions = (
        '.tar', '.tar.gz', '.tar.bz2', '.tar.xz',
        '.so', '.dylib', '.la', '.a', '.dll', '.pyd',
    )
    for root, unused_dirs, files in os.walk(recipe_dir):
        for fn in files:
            fn_lower = fn.lower()
            if fn_lower.endswith(disallowed_extensions):
                if PEDANTIC:
                    raise RecipeError("found: %s" % fn)
                else:
                    print("Warning: found: %s" % fn)
            path = join(root, fn)
            # only allow small archives for testing
            if (PEDANTIC and fn_lower.endswith(('.bz2', '.gz')) and getsize(path) > 512):
                raise RecipeError("found: %s (too large)" % fn)

    if basename(recipe_dir) == 'icu':
        return

    # check total size od recipe directory (recursively)
    kb_size = dir_size(recipe_dir) / 1024
    kb_limit = 512
    if PEDANTIC and kb_size > kb_limit:
        raise RecipeError("recipe too large: %d KB (limit %d KB)" %
                          (kb_size, kb_limit))

    if PEDANTIC:
        try:
            with open(join(recipe_dir, 'build.sh'), 'rb') as fi:
                data = fi.read()
            if data and not data.decode('utf-8').startswith(('#!/bin/bash\n',
                                                             '#!/bin/sh\n')):
                raise RecipeError("not a bash script: build.sh")
        except IOError:
            pass


def render_jinja2(recipe_dir):
    import jinja2

    loaders = [jinja2.FileSystemLoader(recipe_dir)]
    env = jinja2.Environment(loader=jinja2.ChoiceLoader(loaders))
    template = env.get_or_select_template('meta.yaml')
    return template.render(environment=env)


def validate_recipe(recipe_dir, pedantic=True):
    global PEDANTIC
    PEDANTIC = bool(pedantic)

    meta_path = join(recipe_dir, 'meta.yaml')
    with open(meta_path, 'rb') as fi:
        data = fi.read()
    if PEDANTIC and not all_ascii(data):
        raise RecipeError("non-ASCII in: %s" % meta_path)
    if b'{{' in data:
        if PEDANTIC:
            raise RecipeError("found {{ in %s (Jinja templating not allowed)" %
                              meta_path)
        else:
            data = render_jinja2(recipe_dir)
    else:
        data = data.decode('utf-8')

    check_dir_content(recipe_dir)

    for cfg in iter_cfgs():
        meta = parse(data, cfg)
        validate_meta(meta)
        validate_files(recipe_dir, meta)
