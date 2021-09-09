#!/bin/sh -e

git checkout master
git merge develop
git push
git tag $CBVERSION
git push origin $CBVERSION

# this assumes that build is a child of the top level directory
cd ..
make pip
