require(xtable)
df000 = read.csv("output5k_ra_000.csv", comment.char = "#")
df025 = read.csv("output5k_ra_025.csv", comment.char = "#")
df050 = read.csv("output5k_ra_050.csv", comment.char = "#")
df075 = read.csv("output5k_ra_075.csv", comment.char = "#")
df100 = read.csv("output5k_ra_100.csv", comment.char = "#")

# Distance

results = df000
result_rows = results[grepl("mean distance", results$Measure),]
result_rows = aggregate(result_rows$Value,by=list(result_rows$Algorithm),FUN=mean)
df = data.frame(Algorithm=result_rows$Group.1) 
df["Mean 0.00"] = result_rows$x
df["Ratio to best 0.00"] = result_rows$x / min(result_rows$x)

results = df025
result_rows = results[grepl("mean distance", results$Measure),]
result_rows = aggregate(result_rows$Value,by=list(result_rows$Algorithm),FUN=mean)
df["Mean 0.25"] = result_rows$x
df["Ratio to best 0.25"] = result_rows$x / min(result_rows$x)

results = df050
result_rows = results[grepl("mean distance", results$Measure),]
result_rows = aggregate(result_rows$Value,by=list(result_rows$Algorithm),FUN=mean)
df["Mean 0.50"] = result_rows$x
df["Ratio to best 0.50"] = result_rows$x / min(result_rows$x)

results = df075
result_rows = results[grepl("mean distance", results$Measure),]
result_rows = aggregate(result_rows$Value,by=list(result_rows$Algorithm),FUN=mean)
df["Mean 0.75"] = result_rows$x
df["Ratio to best 0.75"] = result_rows$x / min(result_rows$x)

results = df100
result_rows = results[grepl("mean distance", results$Measure),]
result_rows = aggregate(result_rows$Value,by=list(result_rows$Algorithm),FUN=mean)
df["Mean 1.00"] = result_rows$x
df["Ratio to best 1.00"] = result_rows$x / min(result_rows$x)

df = df[order(df["Ratio to best 0.50"]),]
print(df)
xtable(df,digits=4)

# Rankings

results = df000
result_rows = results[grepl("mean rank", results$Measure),]
result_rows = aggregate(result_rows$Value,by=list(result_rows$Algorithm),FUN=mean)
df = data.frame(Algorithm=result_rows$Group.1) 
df["Mean 0.00"] = result_rows$x
df["Ratio to best 0.00"] = result_rows$x / min(result_rows$x)

results = df025
result_rows = results[grepl("mean rank", results$Measure),]
result_rows = aggregate(result_rows$Value,by=list(result_rows$Algorithm),FUN=mean)
df["Mean 0.25"] = result_rows$x
df["Ratio to best 0.25"] = result_rows$x / min(result_rows$x)

results = df050
result_rows = results[grepl("mean rank", results$Measure),]
result_rows = aggregate(result_rows$Value,by=list(result_rows$Algorithm),FUN=mean)
df["Mean 0.50"] = result_rows$x
df["Ratio to best 0.50"] = result_rows$x / min(result_rows$x)

results = df075
result_rows = results[grepl("mean rank", results$Measure),]
result_rows = aggregate(result_rows$Value,by=list(result_rows$Algorithm),FUN=mean)
df["Mean 0.75"] = result_rows$x
df["Ratio to best 0.75"] = result_rows$x / min(result_rows$x)

results = df100
result_rows = results[grepl("mean rank", results$Measure),]
result_rows = aggregate(result_rows$Value,by=list(result_rows$Algorithm),FUN=mean)
df["Mean 1.00"] = result_rows$x
df["Ratio to best 1.00"] = result_rows$x / min(result_rows$x)

df = df[order(df["Ratio to best 0.50"]),]

print(df)
xtable(df)