d <- read.table("results.csv", header=TRUE, sep="\t");
pdf("runTimeVsN.pdf");
x = 
plot(log2(d$N), 
     log2(d$symmetric), 
     xlab=c("Log2(# tips)"), 
     ylab=c("Log2(Running time in sec)"),
     main=c(""), 
     type="l",
     lty=2);
legend(4, 9, c("Symmetric", "Pectinate"), lty=c(2, 1));

points(log2(d$N), log2(d$symmetric), pch=16)
lines(log2(d$N), log2(d$pectinate), lty=1)
points(log2(d$N), log2(d$pectinate), pch=16)
dev.off();

print (log2(286.362) - log2(0.017))/7
print (log2(3274.848) - log2(0.032))/7
