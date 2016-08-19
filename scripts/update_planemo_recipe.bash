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
#For your reference the SHA256 is: efc829aa1c579c5d8cace3a3da46284aa1a78fbec80f3a5a31a73e08f5d2bc6e
sha_line=`awk '/For your reference the SHA256 is:/' output`
rm output
echo "sha_line is $sha_line"
IFS=" " read -a sha_line_array <<< "$sha_line"
sha=${sha_line_array[6]}
echo "Updating SHA256 to $sha"
sed -i "s/^  sha256.*$/  sha256 \"$sha\"/" planemo.rb

git add planemo.rb
git commit -m "Rev planemo to version $VERSION"
git push origin master
