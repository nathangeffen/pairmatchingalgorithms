print("K analysis")
df = read.csv("output20k_k.csv", comment.char = "#")
df = df[grepl("mean rank",df$Measure),]
results = aggregate(df$Value,by=list(df$k),FUN=mean)
png(filename="EffectivenessK.png")
plot(results, type="b",xlab="k",ylab="Mean ranking")
dev.off()