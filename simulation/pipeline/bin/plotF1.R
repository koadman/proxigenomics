plotF1<-function(ip,fw,fuw) {
	dw <- read.table(fw,col.names=c('fmeas','recall','prec'))
	duw <- read.table(fuw,col.names=c('fmeas','recall','prec'))
	plot(ip,dw$fmeas,type='b',col='red',ylim=c(0,1))
	lines(ip,duw$fmeas,type='b',col='blue',ylim=c(0,1))
	legend('bottomright',legend=c('weighted','unweighted'),col=c('red','blue'),pch=20)
}