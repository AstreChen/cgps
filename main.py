
import subprocess
from optparse import OptionParser

def config_option():
    usage = 'Usage: main.py -e expfile -p phefile -d datatype -o outdir'
    p = OptionParser(usage)
    p.add_option(
        '-e','--expfile',dest='expfile',action='store_true',
        help='')
    p.add_option(
        '-p','--phefile',dest='phefile',action='store_true',
        help='')
    p.add_option(
        '-d','--datatype',dest='datatype',action='store_true',
        help='')
    p.add_option(
        '-o','--outdir',dest='outdir',action='store_true',
        help='')
    opt, args = p.parse_args()
    return p,opt,args 

opt_parser, opt, args = config_option()

if len(sys.argv) == 1:
    opt_parser.print_help()
    sys.exit(1)

cmdline = 'Rscript test1.R '+
' '.join(['-e', opt.expfile, '-p', opt.phefile, '-d', opt.datatype,'-o',opt.outdir])

subprocess.Popen(cmdline)
