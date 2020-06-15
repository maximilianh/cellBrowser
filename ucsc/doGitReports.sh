mkdir /hive/groups/qa/CB-git-reports-history/v0.7.9
mv /hive/groups/qa/CB-git-reports/* /hive/groups/qa/CB-git-reports-history/v0.7.9
git log --tags --simplify-by-decoration --pretty="format:%ai %d"


git-reports v0.7.9 HEAD 2020-04-20 2020-06-12 review `pwd` /hive/groups/qa/CB-git-reports review

cp /hive/groups/qa/CB-git-reports-history/v0.7.9/index.html /hive/groups/qa/CB-git-reports
