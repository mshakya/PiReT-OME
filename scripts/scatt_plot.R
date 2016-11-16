cmd_args = commandArgs(TRUE)

infile1 <- cmd_args[1]
filedir <- cmd_args[2]

infile=paste(filedir, '/', infile1, sep="")

fullset <- read.table(file = infile, sep = "\t", stringsAsFactors = FALSE,  header=TRUE)
sigset <- subset(fullset, fullset$significant=="yes")
insigset <- subset(fullset, fullset$significant=="no")
xname<-colnames(fullset)[2]
yname<-colnames(fullset)[3]
#xmax<-max(log10(fullset[,2]))
#xmin<-min(log10(fullset[,2]))
#ymax<-max(log10(fullset[,3]))
#ymin<-min(log10(fullset[,3]))

#xsigval<-log10(sigset[,2])
#ysigval<-log10(sigset[,3])
#xinsigval<-log10(insigset[,2])
#yinsigval<-log10(insigset[,3])

xmax<-max(fullset[,2])
xmin<-min(fullset[,2])
ymax<-max(fullset[,3])
ymin<-min(fullset[,3])

xsigval<-sigset[,2]
ysigval<-sigset[,3]
xinsigval<-insigset[,2]
yinsigval<-insigset[,3]

fileout=paste(filedir, '/figures/', infile1, ".pdf", sep="")
pdf(fileout)

plot(xinsigval, yinsigval, xlab=xname, log = "xy", ylab=yname, col=1, main="scattered RPKM plot", xlim=c(xmin, xmax), ylim=c(ymin, ymax))
#par(new=TRUE)
#plot(xsigval, ysigval, xlab=xname, ylab=yname, log = "xy", col="red", main="scattered RPKM plot", xlim=c(xmin, xmax), ylim=c(ymin, ymax))
#abline(a=0,b=1, lty=2, col="green")
points(xsigval, ysigval,col=2,pch=3)
legend("topleft",c("not significant","significant"),col=c(1,2),pch=c(1,3))
dev.off()
