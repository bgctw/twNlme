# in zz_archive\09\bef\singleTree
# treeData.RData
Wutzler08BeechStem <- rs[is.finite(rs$stem),c("author","stand","alt","si", "age","stockden", "dbh",      "height" ,  "stem" )]
save(Wutzler08BeechStem, file="Wutzler08BeechStem.RData")

#move to data directory of this package
