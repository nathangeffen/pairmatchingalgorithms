df = read.csv("timing_agent_increase_tests.csv", comment.char = "#")
df = df[grepl("algorithm time", df$Measure),]
results = aggregate(df$Value, by=list(df$ID),FUN=mean)
results = data.frame(Agents=results[,1],Speed=round(results[,2]))
require(xtable)
xtable(results)

