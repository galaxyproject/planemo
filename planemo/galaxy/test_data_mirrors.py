"""Fetch remote test data from mirrors that address files by their content hash.

A mirror serves each file at ``<mirror>/<hash function>/<first 2 hex chars>/<hash value>``,
with the hash in lowercase hex, e.g.
``https://depot.galaxyproject.org/iwc-test-data/sha1/9c/9cda46eac7cbb29508975327e410f729a4f1681c``.
Job inputs with a remote ``location`` and declared ``hashes`` are pointed at the
first mirror that has a matching file, and keep their original URL if none has.
Galaxy verifies the declared hashes when fetching, so the mirror does not need to be
trusted. The mirror must be reachable from Galaxy, not just from planemo.
"""

import functools
import os
import re
from typing import (
    Any,
    Dict,
    List,
    Optional,
    Set,
    Tuple,
)

import requests

from planemo.io import info

# Galaxy hash function names, with the mirror directory and hex digest length for each,
# in order of preference
HASH_FUNCTIONS = {
    "SHA-256": ("sha256", 64),
    "SHA-512": ("sha512", 128),
    "SHA-1": ("sha1", 40),
    "MD5": ("md5", 32),
}
# (connect, read) timeouts in seconds
TIMEOUT = (3, 10)


class HashMirrors:
    """Look up remote files on mirrors by the hashes declared for them."""

    def __init__(self, mirrors: List[str], session: Optional[requests.Session] = None) -> None:
        self.mirrors = [mirror if mirror.endswith("/") else f"{mirror}/" for mirror in mirrors]
        self.session = session or requests.Session()
        self._available: Dict[str, bool] = {}
        # mirrors that could not be reached, not asked again
        self._unreachable: Set[str] = set()

    def mirror_url(self, hashes: List[Dict[str, Any]]) -> Optional[str]:
        """Return the URL of a mirror copy of the file with ``hashes``, or None if no mirror has one."""
        declared = {
            entry.get("hash_function"): str(entry.get("hash_value")).lower()
            for entry in hashes
            if isinstance(entry, dict) and entry.get("hash_value")
        }
        for function, (directory, length) in HASH_FUNCTIONS.items():
            value = declared.get(function)
            # the value becomes part of the URL, so only accept a well-formed digest
            if value is None or not re.fullmatch(f"[0-9a-f]{{{length}}}", value):
                continue
            for mirror in self.mirrors:
                url = f"{mirror}{directory}/{value[:2]}/{value}"
                if self._is_available(mirror, url):
                    return url
        return None

    def _is_available(self, mirror: str, url: str) -> bool:
        if mirror in self._unreachable:
            return False
        if url not in self._available:
            try:
                response = self.session.head(url, timeout=TIMEOUT, allow_redirects=True)
                self._available[url] = response.status_code == 200
            except requests.RequestException as e:
                info(f"Not using test data mirror {mirror}: {e}")
                self._unreachable.add(mirror)
                return False
        return self._available[url]

    def rewrite_job(self, job: Dict[str, Any]) -> int:
        """Point remote files in ``job`` at mirror copies, returning how many were rewritten.

        Only Files that galaxy-tool-util stages with their ``name`` and ``hashes`` are
        considered: job inputs, items of lists and elements of collections. Files in
        records and composite data are staged without their hashes and are left alone.
        """
        return sum(self._rewrite_item(value) for value in job.values())

    def _rewrite_item(self, item: Any) -> int:
        if isinstance(item, list):
            return sum(self._rewrite_item(value) for value in item)
        if not isinstance(item, dict):
            return 0
        if item.get("class") == "Collection":
            return sum(self._rewrite_item(element) for element in item.get("elements") or [])
        if item.get("class") == "File":
            return self._rewrite_file(item)
        return 0

    def _rewrite_file(self, file: Dict[str, Any]) -> int:
        if not file.get("hashes") or file.get("composite_data"):
            return 0
        # same choice as galaxy-tool-util's staging
        key = "location" if file.get("location") else "path"
        location = file.get(key)
        if not isinstance(location, str) or not location.startswith(("http://", "https://")):
            return 0
        url = self.mirror_url(file["hashes"])
        if not url:
            return 0
        info(f"Using mirror copy {url} of {location}")
        if not file.get("name"):
            # the dataset name galaxy-tool-util would derive from the original URL
            file["name"] = os.path.basename(location)
        file[key] = url
        return 1


@functools.lru_cache(maxsize=None)
def hash_mirrors(mirrors: Tuple[str, ...]) -> HashMirrors:
    """Return a :class:`HashMirrors` shared by all jobs staged by this process."""
    return HashMirrors(list(mirrors))
