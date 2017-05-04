df = read.csv("positions.csv", comment.char = "#")
x=seq(1,5000)

vals=df[grepl("RKPM", df$Algorithm),]
plot(x, vals$Rank, type="l",col="blue")

vals=df[grepl("DCPM", df$Algorithm),]
plot(x, vals$Rank, type="l",col="green")
 
vals=df[grepl("WSPM", df$Algorithm),]
plot(x, vals$Rank, type="l",col="yellow")
 
vals=df[grepl("CSPM", df$Algorithm),]
plot(x, vals$Rank, type="l",col="red")

vals=df[grepl("BFPM", df$Algorithm),]
plot(x, vals$Rank, type="l",col="violet")

vals=df[grepl("Blossom V", df$Algorithm),]
plot(x, vals$Rank, type="l",col="purple")

# results_mean = aggregate(df$Rank,by=list(df$Algorithm),FUN=mean)
# results_median = aggregate(df$Rank,by=list(df$Algorithm),FUN=median)

