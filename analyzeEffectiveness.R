# First distances
require(xtable)
print("Distances")
df = read.csv("tmp_analysis.csv", comment.char = "#")
results = aggregate(df$Value, by=list(df$Algorithm,df$Measure),FUN=mean)
result_rows = results[grepl("mean distance", results$Group.2),]
df = data.frame(Algorithm=result_rows$Group.1, Mean=result_rows$x)
result_rows = results[grepl("stddev distance", results$Group.2),]
df[["Std. Dev."]]=result_rows$x
result_rows = results[grepl("median distance", results$Group.2),]
df[["Median"]]=result_rows$x
result_rows = results[grepl("25% distance", results$Group.2),]
df[["25%"]]=result_rows$x
result_rows = results[grepl("75% distance", results$Group.2),]
df[["75%"]]=result_rows$x
df[["Effectiveness"]]=df$Mean/min(df$Mean)
result_rows = results[grepl("algorithm time", results$Group.2),]
df[["Time taken"]]=result_rows$x
df=df[order(df$Effectiveness),]
print(df)
xtable(df)
# Ranks
print("*******")
print("Ranks")
result_rows = results[grepl("mean rank", results$Group.2),]
df = data.frame(Algorithm=result_rows$Group.1, Mean=result_rows$x)
result_rows = results[grepl("stddev rank", results$Group.2),]
df[["Std. Dev."]]=result_rows$x
result_rows = results[grepl("median rank", results$Group.2),]
df[["Median"]]=result_rows$x
result_rows = results[grepl("25% rank", results$Group.2),]
df[["25%"]]=result_rows$x
result_rows = results[grepl("75% rank", results$Group.2),]
df[["75%"]]=result_rows$x
df[["Effectiveness"]]=df$Mean/min(df$Mean)
df=df[order(df$Effectiveness),]
print(df)
xtable(df)
print("****************")