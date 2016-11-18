library(optparse)

option_list <- list(
	make_option(
		c("--ifile"), type="character",
		help="Input table file"),
	make_option(
		c("--sfile"), type="character",
		help="Input table file"),	
	make_option(
		c("-o", "--ofile"), type="character",
		help="Output graph file", metavar="File")
)
opt = parse_args(OptionParser(option_list = option_list))

idata = read.delim(opt$ifile, header=FALSE)
sdata = read.delim(opt$sfile, header=FALSE)

pdf(opt$ofile, width = 14, height=14)

plot(
	jitter(idata$V6, amount=0.4), 
	jitter(sdata$V6, amount=0.4), 
	xlab='real alignment scores',
	ylab='shuffled alignment scores',
	pch=16,
	cex=0.5,
	col=rgb(0, 0, 0, 0.25),
	xlim=c(0,35),
	ylim=c(0,35)
)

min(idata$V6)
min(sdata$V6)

max(idata$V6)
max(sdata$V6)

hist(
	sdata$V6,
	xlim=c(0,35),
	col=rgb(1, 0, 0, 0.25),
	breaks=seq(-0.5, 50.5, 1),
	xaxp = c(0, 50, 50),
	freq=TRUE,
	xlab="Blue = real, Red = shuffled",
	main="Alignment Score Distribution"
)

hist(
	idata$V6,
	xlim=c(0,35),
	col=rgb(0, 0, 1, 0.25),
	breaks=seq(-0.5, 50.5, 1),
	add=TRUE,
	freq=TRUE
)

hist(
	idata$V7,
	xlim=c(0,200),
	col=rgb(0, 0, 1, 0.25),
	breaks=seq(-0.5, 200.5, 1),
	xaxp = c(0, 200, 200),
	freq=TRUE,
	xlab="Blue = real, Red = shuffled",
	main="Alignment Center Pos Distribution"
)

hist(
	sdata$V7,
	xlim=c(0,200),
	col=rgb(1, 0, 0, 0.25),
	breaks=seq(-0.5, 200.5, 1),
	add=TRUE,
	freq=TRUE
)

dev.off()