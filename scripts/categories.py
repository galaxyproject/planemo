from __future__ import print_function

import requests


categories = requests.get("https://testtoolshed.g2.bx.psu.edu/api/categories").json()
print("CURRENT_CATEGORIES = [")
for c in map(lambda c: c["name"], categories):
    print('    "%s",' % c)
print("]")
