#!/bin/bash

set -e

HOMEBREW_TAP="galaxyproject/homebrew-tap"
VERSION=`python -c "import xmlrpclib; print xmlrpclib.ServerProxy('https://pypi.python.org/pypi').package_releases('planemo')[0]"`
URL=`python -c "import xmlrpclib; import re; print re.escape([s for s in xmlrpclib.ServerProxy('https://pypi.python.org/pypi').release_urls('planemo', '$VERSION') if s['filename'].endswith('.tar.gz')][0]['url'])"`
SHA256=`sha256sum dist/planemo-$VERSION.tar.gz | cut -d' ' -f1`

if [ ! -d homebrew-tap ];
then
    git clone git@github.com:$HOMEBREW_TAP.git homebrew-tap
fi
cd homebrew-tap

brew uninstall planemo || true
echo $URL
sed -E -i "s/^  url.*$/  url \"$URL\"/" planemo.rb
echo $SHA256
sed -i "s/^  sha256.*$/  sha256 \"$SHA256\"/" planemo.rb
brew install planemo.rb > output

git add planemo.rb
git commit -m "Rev planemo to version $VERSION"
git push origin master
