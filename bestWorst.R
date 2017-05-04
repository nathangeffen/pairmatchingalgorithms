require(xtable)
print("Best and worst results")
df = read.csv("tmp_analysis.csv", comment.char = "#")
results = df[grepl("mean rank",df$Measure),]
min_results = aggregate(results$Value, by=list(results$Algorithm),FUN=min)
max_results = aggregate(results$Value, by=list(results$Algorithm),FUN=max)
sd_results = aggregate(results$Value, by=list(results$Algorithm),FUN=sd)
mean_results = aggregate(results$Value, by=list(results$Algorithm),FUN=mean)
df = data.frame(Algorithm=mean_results$Group.1, Mean=mean_results$x)
df["Best"] = min_results$x
df["Worst"] = max_results$x
df["sd"] = sd_results$x
df = df[order(df$Mean),]
print(df)
xtable(df)