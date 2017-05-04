df = read.csv("timing_tests.csv", comment.char = "#")
df = df[grepl("algorithm time", df$Measure),]
results = aggregate(df$Value, by=list(df$Algorithm, df$Measure),FUN=mean)
results = data.frame(Algorithm=results[,1],Speed=results[,3],Speedup=max(results[,3])/results[,3])
results=results[order(-results$Speedup),]
require(xtable)
xtable(results)

