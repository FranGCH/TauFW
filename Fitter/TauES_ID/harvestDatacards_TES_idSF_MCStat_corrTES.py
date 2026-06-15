#! /usr/bin/env python
"""
Description : 
Script to generate datacards for the mutau channel and for each regions defined in the config file. 
The tes is defined as a POI for each "tesRegions" defined in the config file. Horizontal morphing is used
to interpolate between the template genrated in "TESvariations" in config file. 
The tid SF is defined as rateParameter for each "tid_SFRegions" defined in the config file.
The autoMCstat function is used to have bin-by-bin uncertainties for the sum of all backgrounds
"""



import ROOT; ROOT.PyConfig.IgnoreCommandLineOptions = True
import os, sys, re
import yaml
import CombineHarvester.CombineTools.ch as ch
from CombineHarvester.CombineTools.ch import CombineHarvester, MassesFromRange, SystMap, BinByBinFactory, CardWriter, SetStandardBinNames, AutoRebin
import CombineHarvester.CombinePdfs.morphing as morphing
from CombineHarvester.CombinePdfs.morphing import BuildRooMorphing
from CombineHarvester.CombinePdfs.morphing import BuildCMSHistFuncFactory

import ROOT
from ROOT import RooWorkspace, TFile, RooRealVar

def check_integral(filename, procs, name, region):
    from ROOT import TFile, TH1
    file = TFile(filename)
    ret = []
    name = name.replace('$BIN', region)
    for proc in procs:
        histname = region + '/' + proc + '_' + name
        print(histname, '\t', filename)
        histup = file.Get(histname+'Up')
        histdn = file.Get(histname+'Down')
        # Check existence and type before using .Integral()
        if not histup or not hasattr(histup, 'Integral'):
            print(f"ERROR: Histogram '{histname}Up' not found or not a valid TH1 in file {filename}")
            continue
        if not histdn or not hasattr(histdn, 'Integral'):
            print(f"ERROR: Histogram '{histname}Down' not found or not a valid TH1 in file {filename}")
            continue
        if histup and histdn:
            ret.append(proc)
    print('returning: ', name, '\t', ret)
    return ret

def nonempty_regions(input_file, regions, proc="data_obs", min_yield=0.0):
    """Drop regions whose `proc` histogram integral is <= min_yield (e.g. empty DMrest
    at tight WPs, which crashes CombineHarvester on all-zero shape templates).
    Fail-safe: keep regions if the file/dir/hist can't be read, and never drop all."""
    if not input_file or not os.path.exists(input_file):
        return list(regions)
    f = TFile(input_file)
    if not f or f.IsZombie():
        return list(regions)
    keep = []
    for r in regions:
        d = f.Get(r)
        h = d.Get(proc) if d else None
        integ = h.Integral() if (h and h.InheritsFrom("TH1")) else 0.0
        if integ > min_yield:
            keep.append(r)
        else:
            print(">>> [skip-empty] region %r: %s integral=%.4g <= %g -- skipping" % (r, proc, integ, min_yield))
    f.Close()
    return keep if keep else list(regions)

def harvest(setup, year, obs, **kwargs):
    """Harvest cards."""

    channel = setup["channel"].replace("mu","m").replace("tau","t")
    
    tag         = kwargs.get('tag',        ""               )
    extratag    = kwargs.get('extratag',   ""               )
    era         = kwargs.get('era',        '%s-13TeV'%year  )
    analysis    = kwargs.get('analysis',   'ztt'            )
    indir       = kwargs.get('indir',      'input_%s'%year  )
    # corrTES variant: route per-region datacards to the *_corrTES output tree
    # so they don't clobber (or get clobbered by) the uncorrelated workflow.
    if 'input_pt_less_region' in indir:
        outdir  = indir.replace('input_pt_less_region', 'output_pt_less_region_corrTES')
    else:
        outdir  = indir.replace('input', 'output')
    outdir      = os.path.join(outdir, year)
    multiDimFit = kwargs.get('multiDimFit')
    verbosity   = kwargs.get('verbosity')
    outtag      = tag+extratag
   
    filename = "%s/%s_%s_tes_%s.inputs-%s%s.root"%(indir,analysis,channel,obs,era,tag)

    # Auto-skip zero-yield regions (e.g. empty DMrest at tight WPs). Empty all-zero shape
    # templates crash CombineHarvester; dropping them keeps the chain robust per-WP.
    for _key in ("fitRegions", "scanRegions"):
        if _key in setup["observables"][obs]:
            _orig = list(setup["observables"][obs][_key])
            _kept = nonempty_regions(filename, _orig)
            if len(_kept) != len(_orig):
                print(">>> [skip-empty] %s: %d -> %d regions (dropped %s)"
                      % (_key, len(_orig), len(_kept), [r for r in _orig if r not in _kept]))
            setup["observables"][obs][_key] = _kept

    # For each region = DM
    # each variable can have a subset of regions in which it is fitted defined in config file under this variable entry
    if "fitRegions" in setup["observables"][obs]:
      for region in setup["observables"][obs]["fitRegions"]:
        cats = []
        icat = 0
        icat += 1
        cats.append((icat, region))
        # if not given, assume all defined regions should be fitted (be careful with potential overlap!)
        print("region: %s" % (cats))

        signals = [] # ZTT is the signal
        backgrounds = []
        # Make ZTT, ZL and ZJ signals, but only ZTT will get TES and tid_SF
        ztt_signals = [proc for proc in setup["processes"] if "ZTT" in proc]
        other_z_signals = [proc for proc in setup["processes"] if any(x in proc for x in ["ZL","ZJ"])]
        signals = ztt_signals + other_z_signals
        backgrounds = [proc for proc in setup["processes"] if (proc not in signals and "data" not in proc)]
        print("Signals: %s" % signals)
        print("Backgrounds: %s" % backgrounds)

        if("TESvariations" in setup):
          print("Take TESvariations as defined in the config file")
          tesshifts = [ "%.3f"%tes for tes in setup["TESvariations"]["values"] ]
        else:
          print("No TESvariations")
          tesshifts = [ "1.00" ]


        # Creation of the CombineHarvester
        harvester = CombineHarvester()

        # Change flag causing bug : 
        harvester.SetFlag("workspaces-use-clone", True)

        # Add Observation and processes.
        harvester.AddObservations(['*'], [analysis], [era], [channel], cats)
        harvester.AddProcesses(['*'], [analysis], [era], [channel], backgrounds, cats, False)
        # Add ZTT with mass points (TES morphing), and ZL/ZJ as signals without mass points
        harvester.AddProcesses(tesshifts, [analysis], [era], [channel], ztt_signals, cats, True)
        if other_z_signals:
            harvester.AddProcesses(['*'], [analysis], [era], [channel], other_z_signals, cats, True)

        print(green("\n>>> defining nuissance parameters ..."))
  
        if "systematics" in setup:
          for sys in setup["systematics"]:
            sysDef = setup["systematics"][sys]
            scaleFactor = 1.0  
            if "scaleFactor" in sysDef:
              scaleFactor = sysDef["scaleFactor"]
            if "name" in sysDef: sysDef["processes"] = check_integral(filename, sysDef["processes"], sysDef["name"], region)
            harvester.cp().process(sysDef["processes"]).AddSyst(harvester, sysDef["name"] if "name" in sysDef else sys, sysDef["effect"], SystMap()(scaleFactor))
            #print sysDef

        # Adding id SF as a rate parameter affecting only ZTT
        listbin = region.split("_")
        tid_name = "tid_SF"
        # Recommended to define tid_SFRegions in the config file 
        if ("tid_SFRegions" in setup): 
          if region in setup["tid_SFRegions"]:
            tid_name = "tid_SF_%s" %(region)
          else:
              found_match = False
              for bin in listbin:
                #print("bin = %s" %(bin))
                if bin in setup["tid_SFRegions"]:
                  tid_name = "tid_SF_%s" % bin
                  found_match = True
                  break
              if not found_match:
                print('region: ', region)
                print('ERROR : wrong definition of tid_SFRegions for in the config file')
        else: 
          if len(listbin) == 1: # Example : DM
            tid_name = "tid_SF_%s"%(listbin[0]) # tid_SF_DM
          else: # Example : DM_pt
            tid_name = "tid_SF_%s"%(listbin[1]) # tid_SF_pt
        
        print("tid : %s" % (tid_name))
        # Add SF only to ZTT signals
        harvester.cp().process(ztt_signals).AddSyst(harvester, tid_name,'rateParam', SystMap()(1.00))

        
        # Add W+Jets SF as a free parameter 
        if not "norm_wj" in setup["systematics"]:
          print("W+Jets SF as a free parameter ")
          sf_W = "sf_W_%s"%(region)
          harvester.cp().process(['W']).AddSyst(harvester, sf_W,'rateParam', SystMap()(1.00))
          print(">>>Add sf_W : %s" % (sf_W))
        
        # Add DY cross section as a free parameter. Don't forgot to add Zmm CR !
        if not "xsec_dy" in setup["systematics"]:
          print("DY cross section as a free parameter")
          harvester.cp().process(['ZTT','ZL','ZJ']).AddSyst(harvester, "xsec_dy" ,'rateParam', SystMap()(1.00))
        def scaleProcess(process,scale): 
          """Help function to scale a given process."""
          process.set_rate(process.rate()*scale)

        # if "scaleFactors"  in setup and "xsec_dy" in setup["scaleFactors"]:
        #   print("DY cross section from config file")
        #   xsec_def = setup["scaleFactors"]["xsec_dy"]
        #   # rate = hist.GetBinContent(1)
        #   harvester.cp().scaleProcess(xsec_def["processes"], xsec_def["value"])

        # Add DY cross section as a free parameter. Don't forgot to add Zmm CR !
        
        # print("muon fake rate free parameter")
        # harvester.cp().process(['ZL', 'TTL']).AddSyst(harvester, "muonFakerate" ,'rateParam', SystMap()(1.00))



        # EXTRACT SHAPES
        print(green(">>> extracting shapes..."))
        print(">>>   file %s" % (filename))
        # Extract shapes: backgrounds normal, ZTT with TES mass templates, ZL/ZJ normal templates
        harvester.cp().channel([channel]).backgrounds().ExtractShapes(filename, "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC")
        # ZTT: use TES templates if present
        
	############################
	#HERE MIGHT BE PROBLEMATIC
        if ztt_signals:
            if("TESvariations" in setup):
                harvester.cp().process(ztt_signals).ExtractShapes(filename, "$BIN/$PROCESS_TES$MASS", "$BIN/$PROCESS_$SYSTEMATIC")  #$BIN/$PROCESS_TES$MASS_$SYSTEMATIC  #$BIN/$PROCESS_TES$MASS_$SYSTEMATIC
            else:
                harvester.cp().process(ztt_signals).ExtractShapes(filename, "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC")
		# ZL/ZJ: regular templates (no TES)
        if other_z_signals:
            harvester.cp().process(other_z_signals).ExtractShapes(filename, "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC")
# ...existing code...
        # backgrounds: normal templates (no TES)
        # ZTT: use TES templates for the nominal, but read systematics from non‑TES names
        # harvester.cp().process(ztt_signals).ExtractShapes(
        #     filename,
        #     "$BIN/$PROCESS_TES$MASS",   # e.g. DM0/ZTT_TES1.000
        #     "$BIN/$PROCESS_$SYSTEMATIC" # e.g. DM0/ZTT_shape_dy_DM0Up
        # )
        # # ZL/ZJ: regular templates (no TES)
        # harvester.cp().process(other_z_signals).ExtractShapes(filename, "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC")
# ...existing code...
############################"$BIN/$PROCESS_TES$MASS", "$BIN/$PROCESS_$SYSTEMATIC")


       # ROOVAR
        workspace = RooWorkspace(analysis,analysis)
        
#print analysis

        # Adding TES as a POI the signal ZTT
        # CORRELATED-TES variant: drop the pT suffix so all pT bins of one DM share one TES.
        # e.g. region "DM0_pt1" -> tes_name "tes_DM0"
        dm_part = region.split('_')[0]
        tes_name = "tes_%s" % dm_part
        print("tes: %s"%(tes_name))

        if("TESvariations" in setup):
          #print("TESvariations")
          tes = RooRealVar(tes_name,tes_name, min(setup["TESvariations"]["values"]), max(setup["TESvariations"]["values"]))
        else:
          tes = RooRealVar(tes_name,tes_name, 1.000, 1.000)

        #tes = RooRealVar(tes_name,tes_name, min(setup["TESvariations"]["values"]), max(setup["TESvariations"]["values"]))
        tes.setConstant(True)

    
        # MORPHING (only for ZTT)
        print(green(">>> morphing..."))
        BuildCMSHistFuncFactory(workspace, harvester, tes, "ZTT")
    
        #workspace.Print()
        workspace.writeToFile("workspace_py.root")


        # EXTRACT PDFs
        print(green(">>> add workspace and extract pdf..."))
        harvester.AddWorkspace(workspace, False)
        harvester.ExtractPdfs(harvester, "ztt", "$BIN_$PROCESS_morph", "")  # Extract all processes (signal and bkg are named the same way)
        
        harvester.ExtractData("ztt", "$BIN_data_obs")  # Extract the RooDataHist

        harvester.SetAutoMCStats(harvester, 0, 1, 1) # Set the autoMCStats line (with -1 = no bbb uncertainties)



        # NUISANCE PARAMETER GROUPS
        # To do: export to config file
        print(green(">>> setting nuisance parameter groups..."))
        harvester.SetGroup('all', [ ".*"           ])
        harvester.SetGroup('sys', [ "^((?!bin).)*$"]) # everything except bin-by-bin
        harvester.SetGroup( 'bin',      [ ".*_bin.*"        ])
        harvester.SetGroup( 'lumi',     [ ".*lumi"           ])
        harvester.SetGroup( 'eff',      [ ".*eff_.*"         ])
        harvester.SetGroup( 'jtf',      [ ".*jTauFake.*"     ])
        harvester.SetGroup( 'ltf',      [ ".*mTauFake.*"     ])
        harvester.SetGroup( 'zpt',      [ ".*shape_dy.*"     ])
        harvester.SetGroup( 'shape_ttbar', [ ".*shape_ttbar.*"  ])
        harvester.SetGroup( 'xsec',     [ ".*xsec.*"         ])
        harvester.SetGroup( 'norm',     [ ".*(lumi|Xsec|Norm|norm_qcd).*" ])
        harvester.SetGroup( 'tid',      [ ".*tid.*"          ])
        harvester.SetGroup( 'tes',      [ ".*tes.*"          ])


        #PRINT
        if int(verbosity) > 0:
            print(green("\n>>> print observation...\n"))
            harvester.PrintObs()
            print(green("\n>>> print processes...\n"))
            harvester.PrintProcs()
            print(green("\n>>> print systematics...\n"))
            harvester.PrintSysts()
            print(green("\n>>> print parameters...\n"))
            harvester.PrintParams()
            print("\n")
    
        # WRITER
        print(green(">>> writing datacards..."))
        datacardtxt  = "$TAG/$ANALYSIS_$CHANNEL_%s-%s%s-$ERA.txt"%(obs,region,outtag)
        datacardroot = "$TAG/$ANALYSIS_$CHANNEL_%s-%s%s-inputs.$ERA.root"%(obs,region,outtag)
        writer = CardWriter(datacardtxt,datacardroot)
        writer.SetVerbosity(verbosity)
        writer.SetWildcardMasses([ ])
        writer.WriteCards(outdir, harvester)
    
        # REPLACE bin ID by bin name
        for region, DM in cats:
          oldfilename = datacardtxt.replace('$TAG',outdir).replace('$ANALYSIS',analysis).replace('$CHANNEL',channel).replace('$BINID',str(region)).replace('$ERA',era)
          newfilename = datacardtxt.replace('$TAG',outdir).replace('$ANALYSIS',analysis).replace('$CHANNEL',channel).replace('$BINID',DM).replace('$ERA',era)
          if os.path.exists(oldfilename):
            os.rename(oldfilename, newfilename)
            print('>>> renaming "%s" -> "%s"' % (oldfilename, newfilename))
          else:
            print('>>> Warning! "%s" does not exist!' % (oldfilename))
        

  
def setYield(process,file,dirname,scale=1.):
    """Help function to get yield from file."""
    histname = "%s/%s"%(dirname,process.process()) if dirname else process.process()
    hist = file.Get(histname)
    if not hist:
        print('setYield: Warning! Did not find histogram "%s" in "%s"' % (histname, file.GetName()))
    if hist.GetXaxis().GetNbins()>1:
        print('setYield: Warning! Histogram "%s" has more than one bin!' % (histname))
    rate = hist.GetBinContent(1)
    process.set_rate(rate*scale)
  
def green(string,**kwargs):
    return kwargs.get('pre',"")+"\x1b[0;32;40m%s\033[0m"%string

def ensureDirectory(dirname):
    """Make directory if it does not exist."""
    if not os.path.exists(dirname):
        os.makedirs(dirname)
        print(">>> made directory " + dirname)
    return dirname



def main(args):

    ## Open and import information from config file here to be publicly accessible in all functions
    print("Using configuration file: %s" % args.config)
    with open(args.config, 'r') as file:
        setup = yaml.safe_load(file)

    verbosity = 1 if args.verbose else 0
    observables = []
    for obs in setup["observables"]:
        observables.append(obs)
    
    indir = args.input_dir
    if args.multiDimFit:
        args.extratag += "_MDF"

    tag = setup["tag"] if "tag" in setup else ""
    print("producing datacards for %s" % (args.year))
    for obs in observables:
        print("producing datacards for %s" % (obs))
        harvest(setup,args.year,obs,tag=tag,extratag=args.extratag,indir=indir,multiDimFit=args.multiDimFit,verbosity=verbosity)
    



if __name__ == '__main__':

  from argparse import ArgumentParser
  argv = sys.argv
  description = '''This script makes datacards with CombineHarvester.'''
  parser = ArgumentParser(prog="harvesterDatacards_TES",description=description,epilog="Succes!")
  parser.add_argument('-y', '--year', dest='year', type=str, default=2018, action='store', help="select year")
  parser.add_argument('-c', '--config', dest='config', type=str, default='TauES/config/defaultFitSetupTES_mutau.yml', action='store', help="set config file containing sample & fit setup")
  parser.add_argument('-e', '--extra-tag', dest='extratag', type=str, default="", action='store', metavar='TAG', help="extra tag for output files")
  parser.add_argument('-M', '--multiDimFit', dest='multiDimFit', default=False, action='store_true', help="assume multidimensional fit with a POI for each DM")
  parser.add_argument('-v', '--verbose', dest='verbose', default=False, action='store_true', help="set verbose")
  parser.add_argument('-i', '--input_dir', dest='input_dir', type=str, help='input_dir for datacard root file')
  args = parser.parse_args()

  main(args)
  print(">>>\n>>> done harvesting\n")


