#!/bin/bash

HOMEBREW_TAP="galaxyproject/homebrew-tap"
VERSION=$1
shift

if[ -d homebrew-tap ];
then
    git clone git@github.com:$HOMEBREW_TAP.git homebrew-tap
fi
cd homebrew-tap

brew uninstall planemo
sed -E -i "s/planemo-([0-9]+)\.([0-9]+)\.([0-9]+)\.tar\.gz/planemo-$VERSION.tar.gz/" planemo.rb
sed -i "s/^  sha256.*$/  sha256 \"\"/" planemo.rb
brew install planemo.rb > output
#For your reference the SHA256 is: efc829aa1c579c5d8cace3a3da46284aa1a78fbec80f3a5a31a73e08f5d2bc6e
sha_line=`awk '/For your reference the SHA256 is:/' output`
echo "sha_line is $sha_line"
IFS=" " read -a sha_line_array <<< "$sha_line"
sha=${sha_line_array[6]}
echo "Updating SHA256 to $sha"
sed -i "s/^  sha256.*$/  sha256 \"$sha\"/" planemo.rb

git add planemo.rb
git commit -m "Rev planemo to version $VERSION"
git push
