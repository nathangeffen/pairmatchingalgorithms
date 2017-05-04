print("C analysis")
df = read.csv("output20k_c.csv", comment.char = "#")
df = df[grepl("mean rank",df$Measure),]
results = aggregate(df$Value,by=list(df$c),FUN=mean)
png(filename="EffectivenessC.png")
plot(results, type="b",xlab="c",ylab="Mean ranking")
dev.off()