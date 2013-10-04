import sys
import gbstools

script, vcf = sys.argv

reader = gbstools.Reader(vcf, disp=2.5, snpmode=True)
writer = gbstools.Writer(sys.stdout, template=vcf, disp=2.5, lineterminator='\n')

for snp in reader:
   if snp.param:
      while not snp.param[-1]['converged'] and not snp.param[-1]['fail']:
         snp.update_param(snp.param)
         if len(snp.param) > 40:
            snp.param[-1]['fail'] = True
      while not snp.null_param[-1]['converged'] and not snp.null_param[-1]['fail']:
         snp.update_param(snp.null_param)
         if len(snp.null_param) > 40:
            snp.null_param[-1]['fail'] = True
      snp.lik_ratio = snp.likelihood_ratio()
      snp.update_info()
   writer.write_record(snp)
   sys.stdout.flush()
