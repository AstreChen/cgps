import sys, os
import subprocess
from optparse import OptionParser

def get_dir():
    srcf = os.environ['HOME']+'/.cgpsrc'
    if not os.path.isfile(srcf):
        sys.exit('please add .cgpsrc under HOME folder\n')
    with open(srcf,'r') as f:
        tmp = f.readline()
    cgps_home = tmp.split('=')[1].strip()
    return cgps_home

def config_option():
    usage = 'Usage: main.py -e expfile -p phefile -d datatype -s species -o outdir'
    p = OptionParser(usage)
    p.add_option(
        '-e','--expfile',dest='expfile',action='store',
        help='')
    p.add_option(
        '-p','--phefile',dest='phefile',action='store',
        help='')
    p.add_option(
        '-d','--datatype',dest='datatype',action='store',
        help='')
    p.add_option(
        '-s','--species',dest='spe',action='store',
        help='')
    p.add_option(
        '-o','--outdir',dest='outdir',action='store',
        help='')
    opt, args = p.parse_args()
    return p,opt,args 

cgps_home = get_dir()
opt_parser, opt, args = config_option()
print opt.expfile,opt.phefile,opt.datatype,opt.outdir
if len(sys.argv) == 1:
    opt_parser.print_help()
    sys.exit(1)

cmdline = 'Rscript combined_methods.R '+ ' '.join([cgps_home,opt.expfile, opt.phefile, opt.datatype, opt.species, opt.outdir])

subprocess.Popen(cmdline)


