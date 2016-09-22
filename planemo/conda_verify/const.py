LICENSE_FAMILIES = set("""
AGPL
GPL2
GPL3
LGPL
BSD
MIT
Apache
PSF
Public-Domain
Proprietary
Other
""".split())

FIELDS = {
    'package': {'name', 'version'},
    'source': {'fn', 'url', 'md5', 'sha1', 'sha256',
               'git_url', 'git_tag', 'git_branch',
               'patches', 'hg_url', 'hg_tag'},
    'build': {'features', 'track_features',
              'number', 'entry_points', 'osx_is_app', 'noarch',
              'preserve_egg_dir', 'win_has_prefix', 'no_link',
              'ignore_prefix_files', 'msvc_compiler',
              'detect_binary_files_with_prefix',
              'always_include_files'},
    'requirements': {'build', 'run'},
    'app': {'entry', 'icon', 'summary', 'type', 'cli_opts'},
    'test': {'requires', 'commands', 'files', 'imports'},
    'about': {'license', 'license_url', 'license_family', 'license_file',
              'summary', 'description', 'home', 'doc_url', 'dev_url'},
}

MAGIC_HEADERS = {
    '\xca\xfe\xba\xbe': 'MachO-universal',
    '\xce\xfa\xed\xfe': 'MachO-i386',
    '\xcf\xfa\xed\xfe': 'MachO-x86_64',
    '\xfe\xed\xfa\xce': 'MachO-ppc',
    '\xfe\xed\xfa\xcf': 'MachO-ppc64',
    'MZ\x90\x00': 'DLL',
    '\x7fELF': 'ELF',
}

DLL_TYPES = {
    0x0: 'UNKNOWN', 0x1d3: 'AM33', 0x8664: 'AMD64', 0x1c0: 'ARM',
    0xebc: 'EBC', 0x14c: 'I386', 0x200: 'IA64', 0x9041: 'M32R',
    0x266: 'MIPS16', 0x366: 'MIPSFPU', 0x466: 'MIPSFPU16', 0x1f0: 'POWERPC',
    0x1f1: 'POWERPCFP', 0x166: 'R4000', 0x1a2: 'SH3', 0x1a3: 'SH3DSP',
    0x1a6: 'SH4', 0x1a8: 'SH5', 0x1c2: 'THUMB', 0x169: 'WCEMIPSV2',
}
