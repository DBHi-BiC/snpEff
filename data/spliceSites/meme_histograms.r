
dens <- function( d, th, title ) {
	q <- quantile(d,th)
	sub=paste('Threshold:', q);
	plot(density(d), main=title, xlab=sub);
	abline( v=q, col='red')
}

d.entropy <- as.vector( read.table('donor.entropy.txt', header = FALSE)[[1]] )
a.entropy <- as.vector( read.table('acceptor.entropy.txt', header = FALSE)[[1]] )
d.prob <- as.vector( read.table('donor.prob.txt', header = FALSE)[[1]] )
a.prob <- as.vector( read.table('acceptor.prob.txt', header = FALSE)[[1]] )
branch.score <- as.vector( read.table('branchScores.txt', header = FALSE)[[1]] )

par( mfrow=c(2,2) );
dens(d.entropy, 0.05, 'Donor: Entropy');
dens(a.entropy, 0.05, 'Acceptor: Entropy');

dens(d.prob, 0.95, 'Donor: Probabilities');
dens(a.prob, 0.95, 'Acceptor: Probabilities');

par( mfrow=c(1,1) );
dens(branch.score, 0.99, 'Best U12 branch scores');
