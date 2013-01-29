

file <- 'testHg3763Chr1.branchDonorScore.txt';
file <- 'GRCh37.66.branchDonorScore.txt';
file <- 'hg19.branchDonorScore.txt';

d <- read.table(file, header = FALSE, col.names=c('score','rand','rand2') );

plot( density( d$rand ), main='PWM score histogram', sub='Donor: Red, Random: Blue and Black', xlab='Score', col='black');
lines( density( d$rand2 ), col='blue');
lines( density( d$score ), col='red');

