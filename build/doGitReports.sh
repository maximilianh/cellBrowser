#!/bin/sh -e

directory=/hive/groups/qa/CB-git-reports-history/$CBVERSION

if test -d $directory
then
    echo "Reports have already been done for $CBVERSION.  If you want to re-do them then:"
    echo "rm -rf /hive/groups/qa/CB-git-reports-history/$CBVERSION"
    false
fi

mkdir  $directory
indexFile=$directory/index.html

echo "<html><head><title>git-reports</title></head><body>" > $indexFile
echo "<h1>GIT changes: cell browser</h1>" >> $indexFile
echo "<ul>" >> $indexFile
echo "<li><a href=review/index.html>review</a> ($CBLASTDATE to $CBDATE)" >> $indexFile
echo "<li><a href=/CB-git-reports-history/>Previous versions</a>" >> $indexFile
echo "</body></html>" >> $indexFile

#mv /hive/groups/qa/CB-git-reports/* /hive/groups/qa/CB-git-reports-history/$CBLASTVERSION
#git log --tags --simplify-by-decoration --pretty="format:%ai %d"

git-reports $CBLASTVERSION HEAD $CBLASTDATE $CBDATE review .. $directory review

rm -rf /hive/groups/qa/CB-git-reports/*
cp -r $directory/* /hive/groups/qa/CB-git-reports
find /hive/groups/qa/CB-git-reports -exec touch {} \;

#cp /hive/groups/qa/CB-git-reports-history/v0.7.9/index.html /hive/groups/qa/CB-git-reports
