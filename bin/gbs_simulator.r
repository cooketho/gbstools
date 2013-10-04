# gbs_simulator.r
# simulate GBS SNPs and output in vcf format

# Draw genotypes
# Draw reads
# calculate likelihoods
# format output

draw_reads <- function(g, dp){
  ref <- rnbinom(n=1, size=psi, mu=mu * g[1]/2)
  nonref <- rnbinom(n=1, size=psi, mu=mu * g[2]/2)
  l_homref <- (1 - e)^ref * (e/3)^nonref
  l_het <- (((1 - e) + e/3)/2)^ref * (((1 - e) + e/3)/2)^nonref
  l_homnonref <- (e/3)^ref * (1 - e)^nonref
  lik <- c(l_homref, l_het, l_homnonref)
  lik <- lik/max(lik)
  lik <- round(-10 * log(lik, base=10))
  return(c(ref, nonref, lik))
}

genotype <- function (g){
  if (all(g == c(2,0,0))){
    genotype <- c(0, 0)
  } else if (all(g == c(1,1,0))){
    genotype <- c(0, 1)
  } else if (all(g == c(1,0,1))){
    genotype <- c(0, 2)
  } else if (all(g == c(0,2,0))){
    genotype <- c(1, 1)
  } else if (all(g == c(0,1,1))){
    genotype <- c(1, 2)
  } else if (all(g == c(0,0,2))){
    genotype <- c(2, 2)
  }
  return(genotype)
}

N <- 1000
n <- 6
mu <- 40
psi <- 25
e <- 0.01

sim <- data.frame(matrix(rep(NA, N*(n + 9)), nrow=N, ncol=(n + 9)))
for (i in 1:N){
  dropout <- runif(1, 0, 0.25)
  # dropout <- 0
  freq <- runif(1, 0, 1 - dropout)
  phi <- c(freq, 1 - dropout - freq, dropout)
  genotypes <- t(rmultinom(n=n, size=2, prob=phi))
  samples <- t(apply(genotypes, 1, function (g) {draw_reads(g=g, dp=mu)}))
  gt <- t(apply(genotypes, 1, function (g) {genotype(g)}))
  gt <- apply(gt, 1, function (x) {paste(x[1], x[2], sep='/')})
  gt[gt=='2/2'] <- './.'
  dp <- rowSums(samples[,1:2])  
  ad <- apply(samples[,1:2], 1, function (x) {paste(x[1], x[2], sep=',')})
  ad[dp==0] <- '.'
  pl <- apply(samples[,3:5], 1, function (x) {paste(x[1], x[2], x[3], sep=',')})
  pl[dp==0] <- '.'
  dp[dp==0] <- '.'
  samples <- paste(gt, ad, dp, pl, sep=':')
  
  nondropout <- genotypes <- genotypes[genotypes[,3]==0,]
  if (class(nondropout_genotypes)=='integer'){
    ac1 <- 
  ac1 <- colSums(genotypes[genotypes[,3]==0,])
  dropout_genotypes = genotypes[genotypes[,3]==1,]
  if (class(dropout_genotypes)=='integer'){
    ac2 <- 2*dropout_genotypes
  } else {
    ac2 <- 2*colSums(dropout_genotypes)}

    
  ac_apparent <- ac1 + ac2
  af_apparent <- sprintf('%.2f', ac_apparent[2]/sum(ac_apparent[1:2]))
  af_apparent <- paste('AF', af_apparent, sep='=')
  af_true <- sprintf('%.2f', colSums(genotypes, na.rm=T)/sum(genotypes))
  af_true <- paste('AFTRUE', paste(af_true[1], af_true[2], af_true[3], sep=','), sep='=')
  info <- paste(af_apparent, af_true, sep=';')
  site <- c(1, i, '.', 'A', 'T,N', '.', '.', info, 'GT:AD:DP:PL')
  site <- c(site, samples)
  sim[i,] <- site
}
names(sim) <- c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 1:n)
write.table(sim, 'sim.null.vcf', col.names=T, row.names=F, sep='\t', quote=F)
