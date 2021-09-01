#!/bin/sh -e

git checkout master
git merge develop
git push
git tag $CBRELEASE
git push origin $CBRELEASE
make pip
