import argparse
import re
import math
from abc import ABC, abstractmethod

# Globals ----------------------------------------------------------------

MODES = ['sclust', 'sequenza', 'titan'] # add new seg filetypes here
LOH_TYPES = ['neutral', 'deletion', 'any']

# Classes ----------------------------------------------------------------

class Parser:
    ''' Extend this class and implement is_header(), parse_segment() methods.
        get_loh_flag() is optional. 
    '''
    def __init__(self, stream, sample, loh_type):
        self.stream   = open(stream, 'r')
        self.sample   = sample
        self.loh_type = loh_type

    @abstractmethod
    def is_header(self, line):
        ''' Return true if line is part of header '''
        pass

    @abstractmethod
    def parse_segment(self, line):
        ''' Return Segment object '''
        pass

    def get_loh_flag(self):
        ''' Return string indicating if segment is LOH or not. '''
        pass

    def calculate_logratio(self, cn):
        cn   = int(cn)
        logr = None
        if cn == 0:
            logr = '-Inf'
        else:
            logr = math.log(cn, 2) - 1
        return(str(logr))

class SequenzaParser(Parser):
    def __init__(self, stream, sample, loh_type):
        super().__init__(stream, sample, loh_type)

    def is_header(self, line):
        if line.split('\t', 1)[0] == "chromosome":
            return True
        if line.split('\t',1)[0] == "\"chromosome\"":
            return True
        else:
            return False

    def parse_segment(self, line):
        _line = line.split('\t')
        chrm, start, end = _line[0:3]
        if chrm.startswith('"'):
            chrm = chrm[1:]
        if chrm.endswith('"'):
            chrm = chrm[0:-1]
        cn, a, b = _line[9:12]
        loh_flag = self.get_loh_flag(cn, a, b)
        logr = self.calculate_logratio(cn)

        return(Segment(chrm, start, end, cn, logr, self.sample, loh_flag))

    def get_loh_flag(self, cn, a, b):
        cn = int(cn)
        loh_flag = '0'
        if self.loh_type == 'neutral':
            if a == 'NA':
                loh_flag = 'NA'
            elif cn == 2 and int(a) == 2:
                loh_flag = '1'
        elif self.loh_type == 'deletion':
            if a == 'NA' or b == 'NA':
                loh_flag = 'NA'
            elif cn == 1 and (int(a) + int(b)) == 1:
                loh_flag = '1'
        elif self.loh_type == 'any':
            if b == 'NA':
                loh_flag = 'NA'
            elif cn <= 2 and int(b) == 0:
                loh_flag = '1'
        return(loh_flag)

class SClustParser(Parser):
    def __init__(self, stream, sample, loh_type):
        super().__init__(stream, sample, loh_type)

    def is_header(self, line):
        return(False)

    def parse_segment(self, line):
        _line = line.split('\t')
        chrm, start, end, x, cn = _line[1:] # not sure what the x column represents
        loh_flag = self.get_loh_flag()
        logr = self.calculate_logratio(cn)
        return(Segment(chrm, start, end, cn, logr, self.sample, loh_flag))

    def get_loh_flag(self):
        return('NA')

class TitanParser(Parser):
    def __init__(self, stream, sample, loh_type):
        super().__init__(stream, sample, loh_type)

    def is_header(self, line):
        p = re.compile('start', re.I)
        m = re.search(p, line)
        return(m)

    def parse_segment(self, line):
        _line = line.split('\t')
        chrm, start, end = _line[1:4]
        logr = _line[6]
        call = _line[8]
        cn = _line[9]
        loh_flag = self.get_loh_flag(call)

        return(Segment(chrm, start, end, cn, logr, self.sample, loh_flag))

    def get_loh_flag(self, call):
        loh_flag = '0'
        if self.loh_type == 'neutral':
            if call == 'NLOH':
                loh_flag = '1'
        elif self.loh_type == 'deletion':
            if call == 'DLOH':
                loh_flag = '1'
        elif self.loh_type == 'any':
            if 'LOH' in call:
                loh_flag = '1'
        return(loh_flag)

class Segment:
    def __init__(self, chrm, start, end, cn, logr, sample, loh_flag = 'NA'):
        self.chrm     = chrm
        self.start    = start
        self.end      = end
        self.loh      = loh_flag
        self.cn       = cn
        self.cn_state = self.get_cnv_state()
        self.logr     = logr
        self.sample   = sample

    def get_cnv_state(self):
        cn_state = 'NEUT'
        cn = int(self.cn)

        if cn > 3:
            cn_state = 'AMP'
        elif cn == 3:
            cn_state = 'GAIN'
        elif cn == 1:
            cn_state = 'HETD'
        elif cn == 0:
            cn_state = 'HOMD'
        return(cn_state)

    def to_igv(self, prepend):
        if prepend:
            self.chrm = "chr"+self.chrm
        return('\t'.join([self.sample, self.chrm, self.start, \
                          self.end, self.loh, self.logr]))

    def to_oncocircos(self, prepend):
        if prepend:
            self.chrm = "chr"+self.chrm
        return('\t'.join([self.sample, self.chrm, self.start, \
                          self.end, self.loh, self.logr, self.cn_state]))

# Argument Parser --------------------------------------------------------

def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('--include_loh', '--loh_type', help = 'Type of LOH events to consider for LOH flag column. Default: any.', 
                        choices = LOH_TYPES, default = 'any', const = 'any', nargs = '?', dest = 'loh_type')
    parser.add_argument('--prepend_chr', help = 'Prepend "chr" to chromosome column.', action = 'store_true')
    parser.add_argument('--oncocircos', help = 'Output OncoCircos compatible output.', action = 'store_true')
    parser.add_argument('--gistic', help = 'Output GISTIC2.0 compatible output.', action = 'store_true')
    parser.add_argument('--mode', help = 'Input file type.', choices = MODES, required = True)
    parser.add_argument('--sample', '--sequenza_sample', help = 'Sample name for IGV ID column.', dest = 'sample', required = True)
    parser.add_argument('seg_file', help = 'Segmentation file to convert to IGV format.')
    return(parser.parse_args())
 
# Main -------------------------------------------------------------------

def main():
    args       = parse_arguments()
    mode       = args.mode
    sample     = args.sample
    seg_file   = args.seg_file
    prepend    = args.prepend_chr
    loh_type   = args.loh_type
    oncocircos = args.oncocircos
    gistic     = args.gistic
    parser     = None

    if gistic:
        prepend = True # backwards compatibility

    # - implement your own parser for any new segment file types by extending Parser
    # - add filetype to MODES list above
    # - add condition for new filetype here and instantiate your own Parser as parser
    if mode == 'sequenza': 
        parser = SequenzaParser(seg_file, sample, loh_type)
    elif mode == 'sclust':
        parser = SClustParser(seg_file, sample, loh_type)
    elif mode == 'titan':
        parser = TitanParser(seg_file, sample, loh_type)

    # header
    if not oncocircos:
        print('ID\tchrom\tstart\tend\tLOH_flag\tlog.ratio')

    for line in parser.stream:
        if parser.is_header(line):
            continue
        seg = parser.parse_segment(line)

        if oncocircos:
            print(seg.to_oncocircos(prepend))
        else:
            print(seg.to_igv(prepend))

if __name__ == '__main__':
    main()
