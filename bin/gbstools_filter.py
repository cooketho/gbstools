import vcf.filters

class GBStoolsFilter(vcf.filters.Base):
    'Filter sites by GBStools restriction site polymorphism likelihood ratio test'

    name = 'GBStools'

    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--critical', type=float, default=2.71, help='maximum passing value for LRT')
        parser.add_argument('--dfreq', type=float, default=0.05, help='maximum passing value for DFreq')
        parser.add_argument('--insertmin', type=float, help='minimum InsMed for test')
        parser.add_argument('--insertmax', type=float, help='maximum InsMed for test')
        parser.add_argument('--insertmad', type=float, default=60.0, help='maximum InsMAD for test')

    def __init__(self, args):
        self.threshold = args.critical
        self.insertmin = args.insertmin
        self.insertmax = args.insertmax
        self.insertmad = args.insertmad
        self.dfreq = args.dfreq

    def __call__(self, record):
        try:
            if 'EMFailH0' in record.INFO:
                pass
            elif 'EMFailH1' in record.INFO:
                pass
            elif record.INFO['InsMed'][0] < self.insertmin:
                pass
            elif record.INFO['InsMed'][0] > self.insertmax:
                pass
            elif record.INFO['InsMAD'][0] > self.insertmad:
                pass
            elif record.INFO['DLR'][0] > self.threshold or record.INFO['DFreq'][0] > self.dfreq:
                return record.INFO['DLR'][0]
        except:
            pass

    def filter_name(self):
        return self.name
