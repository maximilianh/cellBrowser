#!/bin/sh -e

git checkout master
git merge develop
git push
git tag $CBVERSION
git push origin $CBVERSION
make pip
