import ROOT
from array import array
import json

ROOT.gInterpreter.ProcessLine(".O3")
ROOT.ROOT.EnableImplicitMT()
ROOT.gInterpreter.Declare('#include "resolutionhistmaker.cc"')
import os
#from os import listdir
import time
import sys

import subprocess
import glob
import random
import pathlib
import socket
import XRootD.client

def is_zombie(file_path):
    # Try opening the ROOT file and check if it's a zombie file
    #if os.stat(file_path).st_size == 0 :
    #    print(f"WARNING! Found empty file: {file_path}")
    #    return True 
    f = ROOT.TFile.Open(file_path)
    if not f or f.IsZombie() :
        print(f"WARNING! Found zombie file: {file_path}")
        return True
    f.Close()
    return False

def buildFileListPosix(path):
    outfiles = []
    for root, dirs, fnames in os.walk(path):
        for fname in fnames:
            if fname.lower().endswith(".root"):
                if os.stat(f"{root}/{fname}").st_size > 0:
                    outfiles.append(f"{root}/{fname}")
    return outfiles

def buildFileList(path):
    return buildFileListPosix(path)

def makeFilelist(paths, checkFileForZombie=False, maxFiles=0):
    filelist = []
    expandedPaths = []
    for path in paths:
        expandedPaths.append(path)
        print("Reading files from path {p}".format(p=path))
        files = buildFileList(path)
        filelist.extend(files)

    if checkFileForZombie:
        filelist = [p for p in paths if not is_zombie(p)]

    if maxFiles > 0 and len(filelist) > maxFiles:
        filelist = filelist[:maxFiles]

    print(f"Length of list is {len(filelist)} for paths {expandedPaths}")
    if len(filelist) == 0:
        print()
        print("WARNING! 0 files selected, please check")
        print()
        exit(0)

    return filelist

import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i","--input_path", nargs='+', default=[], help="path of the input root files")

parser.add_argument("-o","--output_file", help="name of the output root file",
                    type=str)

parser.add_argument("--maxFiles", help="Maximum number of files, for tests (default is all)",
                        type=int, default=0)

args = parser.parse_args()
tstart = time.time()
cpustrat = time.process_time()

outdir = os.path.dirname(os.path.abspath(args.output_file))
if not os.path.exists(outdir):
    print()
    print(f"Creating folder {outdir} to store outputs")
    os.makedirs(outdir)    
    print()

files=[]

files = makeFilelist(args.input_path, maxFiles=args.maxFiles)

filenames = ROOT.std.vector('string')()

for name in files: filenames.push_back(name)

d = ROOT.RDataFrame("Events", filenames)

d = d.Define("GenMuonBare", "GenPart_status == 1 && (GenPart_statusFlags & 1 || GenPart_statusFlags & (1<<5)) && abs(GenPart_pdgId) == 13")
d = d.Define("GenMuonBare_pt", "GenPart_pt[GenMuonBare]")
d = d.Define("GenMuonBare_eta", "GenPart_eta[GenMuonBare]")
d = d.Define("GenMuonBare_phi", "GenPart_phi[GenMuonBare]")
d = d.Define("GenMuonBare_pdgId", "GenPart_pdgId[GenMuonBare]")
d = d.Define("GenMuonBare_charge", "makecharge(GenMuonBare_pdgId)")

d = d.Define("goodMuons", "Muon_cvhidealPt > 0 && Muon_isGlobal") #MAYBE YOU WANT TO CHANGE THIS
d = d.Define("goodMuons_pt", "Muon_cvhidealPt[goodMuons]")
d = d.Define("goodMuons_eta", "Muon_cvhidealEta[goodMuons]")
d = d.Define("goodMuons_phi", "Muon_cvhidealPhi[goodMuons]")

d = d.Define("isGenMatchedMuonIdx", "hasGenMatchIdx(GenMuonBare_eta, GenMuonBare_phi, goodMuons_eta, goodMuons_phi)")

d = d.Define("resolution", "makeResolution(GenMuonBare_pt, goodMuons_pt, isGenMatchedMuonIdx)")

d = d.Redefine("GenMuonBare_pt","GenMuonBare_pt[isGenMatchedMuonIdx!=-1]")
d = d.Redefine("GenMuonBare_charge","GenMuonBare_charge[isGenMatchedMuonIdx!=-1]")
d = d.Redefine("GenMuonBare_pt","GenMuonBare_charge*GenMuonBare_pt")
d = d.Redefine("GenMuonBare_eta","GenMuonBare_eta[isGenMatchedMuonIdx!=-1]")
d = d.Redefine("GenMuonBare_phi","GenMuonBare_phi[isGenMatchedMuonIdx!=-1]")

binning_eta = array('d',[round(-2.4 + i*0.2,2) for i in range(25)])
binning_pt = array('d',[round(-100. + i*1.,2) for i in range(201)])
binning_reso = array('d',[-0.5+i*0.0005 for i in range(2001)])

model = ROOT.RDF.TH3DModel("histo", "reso-qpt-eta",
                           len(binning_reso)-1, binning_reso,
                           len(binning_pt)-1, binning_pt,
                           len(binning_eta)-1, binning_eta)

f_out = ROOT.TFile(args.output_file, "RECREATE")
histogram = d.Histo3D(model, "resolution", "GenMuonBare_pt", "GenMuonBare_eta")
histogram.Write()
f_out.Close()
