#! /usr/bin/env python
# Author: Izaak Neutelings (August 2020)
# Description: Create input histograms for datacards
#   ./createinputs.py  -y UL2017 -c <config>
import sys
from collections import OrderedDict
sys.path.append("../Plotter/") # for config.samples
# from config.samples import *
from config.samples_v15_SF import *
from TauFW.Plotter.plot.utils import LOG as PLOG, ensuredir
from TauFW.Fitter.plot.datacard import createinputs, plotinputs
from TauFW.Fitter.plot.rebinning import rebinning
import yaml
import time
import os

#nano doc: https://cms-nanoaod-integration.web.cern.ch/autoDoc/NanoAODv12/2022/2023/doc_DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2.html

map_wp_to_int = OrderedDict([('againstjet', 
				OrderedDict([('VVLoose', 2),
               ('VLoose', 3),
               ('Loose',  4),
					     ('Medium', 5),
					     ('Tight',  6),
					     ('VTight', 7)])),
	                    ('againstelectron',
				OrderedDict([('VVLoose', 2),
                                             ('Tight', 6),
                                             ('Loose', 4 )]))
                           ])

map_PNetscores_to_wp = OrderedDict([('againstjet',
          OrderedDict([
                       ('VVLoose',0.1774),
                       ('VLoose', 0.3810),
                       ('Loose',0.6857),
                       ('Medium', 0.8347),
                       ('Tight', 0.9059),
                       ('VTight', 0.9494)])),
                       ('againstelectron',
          OrderedDict([('VVLoose',0.1266),
                      #  ('VLoose',0.6997),
                      #  ('Loose',0.9354),
                      #  ('Medium',0.9791),
                       ('Tight',0.9897)]))
                       ])



def main(args):
  start_time = time.time()  # Start timing
  eras      = args.eras
  parallel  = args.parallel
  verbosity = args.verbosity
  againstjet = args.againstjet
  againstelectron = args.againstelectron
  setupConfFile = args.config
  DM = args.DM
  inclusive = args.inclusive
  plot      = True
  outdir    = ensuredir(f"/eos/user/f/fcasalin/CMSSW_14_1_0_pre4/src/TauFW/Fitter/input_pt_less_region/againstjet_{againstjet}/againstelectron_{againstelectron}")
  plotdir   = ensuredir(outdir,"plots")
  analysis  = 'ztt'

  print("Using configuration file: %s"%setupConfFile  )
  with open(setupConfFile, 'r') as file:
    setup = yaml.safe_load(file)
    
    # FIX: Handle missing 'regions' key (e.g. for mumu config)
    if "regions" not in setup:
        # If fitRegions are used in observables, we might need them. 
        # For mumu config provided, it uses "baseline".
        print("WARNING: 'regions' key missing in config. Adding default 'baseline' region.")
        setup["regions"] = { 'baseline': {'definition': '1', 'title': 'Baseline'} }

    if DM and 'regions' in setup or inclusive:
       newreg = {}
       for reg in setup['regions']:
          if not inclusive:
             if DM+'_' not in reg: continue
             assert(reg not in newreg)
             newreg[reg] = setup['regions'][reg]
          else:
             if 'pt' in reg: continue
             assert(reg not in newreg)
             newreg[reg] = setup['regions'][reg] 
       setup['regions'] = newreg
  channel  = setup["channel"]
  tag      = "13TeV"
  if "tag" in setup:
    tag += setup["tag"]

  for era in eras:
      
    ###############
    #   SAMPLES   #
    ###############
      
    # GET SAMPLESET
    sname     = setup["samples"]["filename"]
    sampleset = getsampleset(channel,era,fname=sname,join=setup["samples"]["join"],split=[],table=False,rmsf=setup["samples"].get("removeSFs",[]),addsf=setup["samples"].get("addSFs",[]))

    # Potentially split up samples in several processes
    if "split" in setup["samples"]:
      for splitSample in setup["samples"]["split"]:
        print("Splitting sample %s into %s"%(splitSample,setup["samples"]["split"][splitSample]))
        sampleset.split(splitSample, setup["samples"]["split"][splitSample])

    # Rename processes according to custom convention
    if "rename" in setup["samples"]:
      for renamedSample in setup["samples"]["rename"]:
        print("Renaming sample %s into %s"%(renamedSample,setup["samples"]["rename"][renamedSample]))
        sampleset.rename(renamedSample,setup["samples"]["rename"][renamedSample])

    # On-the-fly reweighting of specific processes -- do after splitting and renaming!
    if "scaleFactors" in setup:
      for SF in setup["scaleFactors"]:
        SFset = setup["scaleFactors"][SF]
        if not era in SFset["values"]: continue
        print("Reweighting with SF -- %s -- for the following processes: %s"%(SF, SFset["processes"]))
        for proc in SFset["processes"]:
          weight = "( q_1*q_2<0 ? ( "
          for cond in SFset["values"][era]:
            weight += cond+" ? "+str(SFset["values"][era][cond])+" : ("
          weight += "1.0)"
          for i in range(len(SFset["values"][era])-1):
            weight += " )"
          weight +=  ") : 1.0 )"
          print("Applying weight: %s"%weight)
          sampleset.get(proc, unique=True).addextraweight(weight)

    # Name of observed data 
    sampleset.datasample.name = setup["samples"]["data"]

    ################################################
    #   APPLY WORKING POINTS TO BASELINE CUTS      #
    ################################################
    # Map string WP to integer and replace in baseline cuts
    # Assumes config has 'idDeepTau2018v2p5VSjet_2>=5' (Medium) and 'idDeepTau2018v2p5VSe_2>=2' (VVLoose)
    
    if "baselineCuts" in setup:
        # Only apply WP replacement if we are in a channel that likely uses Taus (mutau, etau, etc)
        # or if the cuts are actually present.
        if "tau" in channel or "rawPNet" in setup["baselineCuts"]:
            jetcut = map_PNetscores_to_wp["againstjet"][againstjet]
            electroncut = map_PNetscores_to_wp["againstelectron"][againstelectron]
            
            print(f"Updating baseline cuts for WP: VSjet {againstjet} (idx {jetcut}), VSele {againstelectron} (idx {electroncut})")
            
            # Replace VSjet cut (Default Medium=5)
            if 'rawPNetVSjet_2>=0.8347' in setup["baselineCuts"]:
                setup["baselineCuts"] = setup["baselineCuts"].replace('rawPNetVSjet_2>=0.8347', f'rawPNetVSjet_2>={jetcut}')
            else:
                print("WARNING: Could not find standard VSjet cut 'rawPNetVSjet_2>=0.8347' in baselineCuts to replace!")

            # Replace VSele cut (Default VVLoose=2)
            if 'rawPNetVSe_2>=0.1266' in setup["baselineCuts"]:
                setup["baselineCuts"] = setup["baselineCuts"].replace('rawPNetVSe_2>=0.1266',   f'rawPNetVSe_2>={electroncut}')
            else:
                print("WARNING: Could not find standard VSele cut 'rawPNetVSe_2>=0.1266' in baselineCuts to replace!")
                
            print(f"New baselineCuts: {setup['baselineCuts']}")
        else:
            print(f"Channel '{channel}' does not seem to use Tau ID working points. Skipping WP replacement in baseline cuts.")

     
    ###################
    #   OBSERVABLES   #
    ###################
      
    # Organize observables and regions
    # Structure: [ (Var_obj, [Sel_obj, Sel_obj, ...]) ]
    obs_region_groups = []
    
    for obsName in setup["observables"]:
        obs_config = setup["observables"][obsName]
        default_binning = obs_config["binning"]
        
        # Get target regions for this observable
        target_region_names = obs_config.get("fitRegions", setup.get("regions", {}).keys())
        # Filter to only those present in setup['regions']
        target_region_names = [r for r in target_region_names if r in setup["regions"]]
        
        # Group regions by binning
        regions_by_binning = {} # Key: tuple(binning), Value: list of region names
        
        if not target_region_names:
             print(f"WARNING: No target regions found for observable {obsName}")
             print(f"  fitRegions: {obs_config.get('fitRegions')}")
             print(f"  setup['regions'].keys(): {list(setup['regions'].keys())}")

        for regName in target_region_names:
            reg_config = setup["regions"][regName]
            # Check if region has specific binning
            binning = tuple(reg_config.get("binning", default_binning))
            
            if binning not in regions_by_binning:
                regions_by_binning[binning] = []
            regions_by_binning[binning].append(regName)
            
        # Create Var and Sel objects for each group
        for binning, regNames in regions_by_binning.items():
            # Create Var object
            obsExpr = obs_config.get('variable', obsName)
            extra = obs_config.get("extra", {}).copy()
            if 'filename' not in extra:
                 extra['filename'] = obsName
            
            var_obj = Var(obsExpr, binning[0], binning[1], binning[2], **extra)
            
            # Create Sel objects
            sel_objs = []
            for regName in regNames:
                reg_config = setup["regions"][regName]
                sel_objs.append(
                    Sel(regName, reg_config['title'], setup["baselineCuts"]+" && "+reg_config["definition"])
                )
            
            obs_region_groups.append((var_obj, sel_objs))


    #######################
    #   DATACARD INPUTS   #
    #######################
    # histogram inputs for the datacards
      
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SMTauTau2016
    chshort = channel.replace('tau','t').replace('mu','m') # abbreviation of channel

    fname   = "%s/%s_%s_tes_$OBS%s.inputs-%s-%s.root"%(outdir,analysis,chshort,DM,era,tag)

    print("Nominal inputs")
    
    # Helper function to run createinputs for all groups
    def run_createinputs_all(fname_in, sampleset_in, groups_in, **kwargs):
        for var_obj, sel_objs in groups_in:
             createinputs(fname_in, sampleset_in, [var_obj], sel_objs, **kwargs)

    run_createinputs_all(fname, sampleset, obs_region_groups, filter=setup["processes"], dots=True, parallel=parallel)

    if "TESvariations" in setup:
      for var in setup["TESvariations"]["values"]:
        print("Variation: TES = %f"%var)

        newsampleset = sampleset.shift(setup["TESvariations"]["processes"], ("_TES%.3f"%var).replace(".","p"), "_TES%.3f"%var, " %.1d"%((1.-var)*100.)+"% TES", split=True,filter=False,share=True)
        run_createinputs_all(fname, newsampleset, obs_region_groups, filter=setup["TESvariations"]["processes"], dots=True, parallel=parallel)

    if "systematics" in setup:
      for sys in setup["systematics"]:

        sysDef = setup["systematics"][sys]
        if sysDef["effect"] != "shape":
          continue
        print("Systematic: %s"%sys)

        # Iterate through the variations of the systematic
        for iSysVar in range(len(sysDef["variations"])):

          # Extract relevant parameters for modifying the sample
          sampleAppend = sysDef["sampleAppend"][iSysVar] if "sampleAppend" in sysDef else ""
          weightReplaced = [sysDef["nomWeight"],sysDef["altWeights"][iSysVar]] if "altWeights" in sysDef else ["",""]
          # Create a new sample set with systematic variations
          newsampleset_sys = sampleset.shift(sysDef["processes"], sampleAppend, "_"+sysDef["name"]+sysDef["variations"][iSysVar], sysDef["title"], split=True,filter=False,share=True)
          run_createinputs_all(fname, newsampleset_sys, obs_region_groups, filter=sysDef["processes"], replaceweight=weightReplaced, dots=True, parallel=parallel)

          # Check for overlap with TES variations in setup #### HERE these should be removed
          # if "TESvariations" in setup:
          #   overlap_TES_sys = list( set(sysDef["processes"]) & set(setup["TESvariations"]["processes"]) )
          #   # If overlap exists, apply TES variations
          #   if overlap_TES_sys:
          #     for var in setup["TESvariations"]["values"]:
          #       print("Variation: TES = %f"%var)
          #       newsampleset_TESsys = sampleset.shift(overlap_TES_sys, ("_TES%.3f"%var).replace(".","p")+sampleAppend, "_TES%.3f"%var+"_"+sysDef["name"]+sysDef["variations"][iSysVar], " %.1d"%((1.-var)*100.)+"% TES" + sysDef["title"], split=True,filter=False,share=True)
          #       createinputs(fname,newsampleset_TESsys, observables, bins, filter=overlap_TES_sys, replaceweight=weightReplaced, dots=True, parallel=parallel)

      ############
      #   PLOT   #
      ############
      # control plots of the histogram inputs
      
      if plot:
        pname  = "%s/%s_$OBS_%s-$BIN-%s$TAG%s.png"%(plotdir,analysis,chshort,era,tag)
        text   = "%s: $BIN"%(channel.replace("mu","#mu").replace("tau","#tau_{h}"))
        groups = [ ] #(['^TT','ST'],'Top'),]

        if "mumu" in channel:
            varprocs = OrderedDict([
                       ('Nom',      ['ZL','ZTT', 'ZJ','W','ST','TT','QCD','data_obs'])])
        elif "mutau"in channel:
            varprocs = OrderedDict([
                       ('Nom',      ["ZTT","ZL","ZJ","W","VV","ST","TTT","TTL","TTJ","QCD","data_obs"])])

        # Helper function for plotting per observable
        def run_plotinputs_all(fname_in, varprocs_in, groups_in, **kwargs):
            for var_obj, sel_objs in groups_in:
                plotinputs(fname_in, varprocs_in, [var_obj], sel_objs, **kwargs)

        if obs_region_groups:
            run_plotinputs_all(fname,varprocs,obs_region_groups,text=text,
                       pname=pname,tag=tag,group=groups, parallel=parallel)
            rebinning(fname, obs=obs_region_groups[0][0].filename, tag=tag) 

            original_dir = os.path.dirname(fname)
            original_basename = os.path.basename(fname)
            fname = os.path.join(original_dir, "rebinning", original_basename)
            # fname = fname.split('/')[0] + '/rebinning/' + fname.split('/')[-1]
        else:
            print("WARNING: No observable/region groups found. Skipping plotting.")

        plotdir   = ensuredir(plotdir,"rebinning")
        pname  = "%s/%s_$OBS_%s-$BIN-%s$TAG%s.png"%(plotdir,analysis,chshort,era,tag)
        run_plotinputs_all(fname,varprocs,obs_region_groups,text=text,
                   pname=pname,tag=tag,group=groups, parallel=parallel)
        pname  = "%s/%s_$OBS_%s-$BIN-%s$TAG%s_wqcd_subtracted.png"%(plotdir,analysis,chshort,era,tag)
        varprocs = OrderedDict([
                       ('Nom',      ['ZTT', 'data_obs_nonztt_subtratced'])])
        run_plotinputs_all(fname,varprocs,obs_region_groups,text=text,
                   pname=pname,tag=tag,group=groups, parallel=parallel, mean=True) 
  end_time = time.time()    # End timing
  elapsed = end_time - start_time
  print(f"\n>>> Done. Total runtime: {elapsed:.2f} seconds.")

if __name__ == "__main__":
  from argparse import ArgumentParser
  description = """Create input histograms for datacards"""
  parser = ArgumentParser(prog="createInputs",description=description,epilog="Good luck!")
  parser.add_argument('-y', '--era',     dest='eras', nargs='*', default=['UL2017'], action='store', help="set era" )
  parser.add_argument('-c', '--config', dest='config', type=str, default='TauES/config/defaultFitSetupTES_mutau.yml', action='store',
                                         help="set config file containing sample & fit setup" )
  parser.add_argument('-s', '--serial',  dest='parallel', action='store_false',
                                         help="run Tree::MultiDraw serial instead of in parallel" )
  parser.add_argument('-v', '--verbose', dest='verbosity', type=int, nargs='?', const=1, default=0, action='store',
                                         help="set verbosity" )
  parser.add_argument('-d', '--decaymode', dest='DM', type=str, default='', help="which decay mode" )
  parser.add_argument('-i', '--inclusive', dest='inclusive', action='store_true', default=False, help="which pt region" )
  parser.add_argument('-j', '--jet', dest='againstjet', default='Medium', help="against jet cut")
  parser.add_argument('-e', '--electron', dest='againstelectron', default='VVLoose', help="against electron cut")
  args = parser.parse_args()
  # LOG.verbosity = args.verbosity
  PLOG.verbosity = args.verbosity
  main(args)
  print("\n>>> Done.")
