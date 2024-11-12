import argparse
import os
import sys
import copy
import math
import ROOT
import time
   
parser = argparse.ArgumentParser(description='run')

parser.add_argument('--none', action='store_true'  , help = 'none')
parser.add_argument('--dryrun', action='store_true'  , help = 'dry run')
parser.add_argument('--ntoys', dest = 'ntoys'  , type = int,  default=1, help='')
parser.add_argument('--tag',   default='SmearRealistic' , help = 'algo')
parser.add_argument('--niter', dest = 'niter'  , type = int,  default=1, help='')
parser.add_argument('--forceIter', dest = 'forceIter'  , type = int,  default=-1, help='')

args = parser.parse_args()

def loop_onetoy(seed, toy):    

    tag = args.tag+'_toy'+str(toy)
    cmd_histo_iter0 = './massscales --lumi=-1 --firstIter=0 --lastIter=2 --skipReco '+\
        ' --tag='+tag+' '+\
        ' --run=Iter0 '+\
        ' --nRMSforGausFit=-1.0 '+\
        ' --biasResolution=0.1 ' +\
        ' --minNumEvents=100 --minNumEventsPerBin=30 '+\
        ' --minNumMassBins=4 '+\
        ' --rebin=2 '+\
        ' --fitNorm --fitWidth '+\
        ' --seed='+str(seed)
    if not args.forceIter>0:
        print(cmd_histo_iter0)
    if not (args.dryrun or args.forceIter>0):
        os.system(cmd_histo_iter0)
    cmd_fit_iter0 = './massfit --ntoys=1 --bias=-1 '+\
        '--tag='+tag+' '+\
        '--run=Iter0 '
    if not args.forceIter>0:
        print(cmd_fit_iter0)
    if not (args.dryrun or args.forceIter>0):
        os.system(cmd_fit_iter0)
    cmd_resol_iter0 = './resolfit --ntoys=1 --bias=-1 '+\
        ' --tag='+tag+' '+\
        ' --run=Iter0 '+\
        ' --maxSigmaErr=0.1 '
    if not args.forceIter>0:
        print(cmd_resol_iter0)
    if not (args.dryrun or args.forceIter>0):
        os.system(cmd_resol_iter0)

    for iter in range(1, args.niter+1):
        if args.forceIter>0 and iter!=args.forceIter:
            continue
        cmd_histo_iteri = cmd_histo_iter0.replace('--run=Iter0', '--run=Iter'+str(iter))
        cmd_histo_iteri += ' --usePrevFit '+\
            ' --tagPrevFit='+tag+' '+\
            ' --runPrevFit=Iter'+str(iter-1)+' '
        cmd_histo_iteri += ' --useSmearFit '+\
            ' --tagSmearFit='+tag+' '+\
            ' --runSmearFit=Iter'+str(iter-1)+' '
        print(cmd_histo_iteri)
        if not args.dryrun:
            os.system(cmd_histo_iteri)
        cmd_fit_iteri = cmd_fit_iter0.replace('--run=Iter0', '--run=Iter'+str(iter))
        print(cmd_fit_iteri)
        if not args.dryrun:
            os.system(cmd_fit_iteri)
        cmd_resol_iteri = cmd_resol_iter0.replace('--run=Iter0', '--run=Iter'+str(iter))
        print(cmd_resol_iteri)
        if not args.dryrun:
            os.system(cmd_resol_iteri)
    return


if __name__ == '__main__':
    iseed = 4357
    start = time.time()
    for itoy in range(0, args.ntoys):
        iseed += itoy*2
        print('Running toy with seed '+str(iseed))
        loop_onetoy(seed=iseed,toy=itoy)
    if args.ntoys>1:
        for iter in range(0, args.niter+1):
            if args.forceIter>0 and iter!=args.forceIter:
                continue
            cmd_hadd = 'hadd -f massfit_'+args.tag+'_merged_Iter'+str(iter)+'.root massfit_'+args.tag+'_toy*_Iter'+str(iter)+'.root'
            print(cmd_hadd)
            if not args.dryrun:
                os.system(cmd_hadd)
    end = time.time()
    print(args.ntoys, 'toys run in', (end - start)/60., 'min.')
