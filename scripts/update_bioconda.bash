#!/bin/bash

set -e

BIOCONDA_PROJECT="jmchilton/bioconda-recipes"
RECIPE="recipes/planemo"
VERSION=`python -c "import xmlrpclib; print xmlrpclib.ServerProxy('https://pypi.python.org/pypi').package_releases('planemo')[0]"`
URL=`python -c "import xmlrpclib; import re; print re.escape([s for s in xmlrpclib.ServerProxy('https://pypi.python.org/pypi').release_urls('planemo', '$VERSION') if s['filename'].endswith('.tar.gz')][0]['url'])"`
MD5SUM=`md5sum dist/planemo-$VERSION.tar.gz | cut -d' ' -f1`

if [ ! -d bioconda ];
then
    git clone git@github.com:$BIOCONDA_PROJECT.git bioconda-recipes
fi
cd bioconda-recipes

sed -E -i "s/^  version: .*$/  version: \"$VERION\"" $RECIPE/meta.yaml
sed -E -i "s/^  url: .*$/  url: $URL/" $RECIPE/meta.yaml
sed -E -i "s/^  md5: .*$/  md5: $MD5SUM/" $RECIPE/meta.yaml

conda build $RECIPE --channel bioconda --channel r
git add 
git add $RECIPE/meta.yaml
git commit -m "Rev planemo to version $VERSION"
git push origin master
