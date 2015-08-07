#!/usr/bin/env python
import vcf
import sys

header = ('CHROM', 'POS', 'DCount', 'AC', 'ACgbs', 'Hets', 'HetsGBS', 'HetsMissing', 'Missing', 'Background', 'AncestralRS', 'DLR', 'DFreq', 'AFH0', 'AFH1', 'LambdaH0', 'LambdaH1', 'DigestH0', 'DigestH1', 'EMFailH0', 'EMFailH1', 'Dispersion', 'DP')
print '#' + '\t'.join(header)

keys = ('DCount', 'AC', 'ACgbs', 'Hets', 'HetsGBS', 'HetsMissing', 'Missing', 'Background', 'AncestralRS', 'DLR', 'DFreq', 'AFH0', 'AFH1', 'LambdaH0', 'LambdaH1', 'DigestH0', 'DigestH1', 'EMFailH0', 'EMFailH1', 'Dispersion')
reader = vcf.Reader(sys.stdin)
for snp in reader:
    line = []
    line.append(snp.CHROM)
    line.append(snp.POS)
    for key in keys:
        try:
            val = snp.INFO[key]
            if type(val) == list:
                line.append(val[0])
            else:
                line.append(val)
        except:
            line.append('.')
    try:
        dp = [sample['DP'] for sample in snp if sample['DP']]
        line.append(sum(dp))
    except:
        line.append(0)
    print '\t'.join([str(i) for i in line])
