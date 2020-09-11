from TauFW.PicoProducer.storage.Sample import MC as M
from TauFW.PicoProducer.storage.Sample import Data as D
storage  = "/pnfs/psi.ch/cms/trivcat/store/user/$USER/samples/NANOAOD_2018/$PATH"
url      = None #"root://cms-xrd-global.cern.ch/"
filelist = None #"samples/files/2016/$SAMPLE.txt"
samples  = [
  
  # DRELL-YAN
  M('DY','DYJetsToLL_M-50',
    "/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM",
    store=storage,url=url,file=filelist,opts='zpt=True',
  ),
  M('DY','DY1JetsToLL_M-50',
   "/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM",
   store=storage,url=url,file=filelist,opts='zpt=True',
  ),
  M('DY','DY2JetsToLL_M-50',
   "/DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM",
   store=storage,url=url,file=filelist,opts='zpt=True',
  ),
  M('DY','DY3JetsToLL_M-50',
   "/DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM",
   store=storage,url=url,file=filelist,opts='zpt=True',
  ),
  M('DY','DY4JetsToLL_M-50',
   "/DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM",
   store=storage,url=url,file=filelist,opts='zpt=True',
  ),
  M('DY','DYJetsToLL_M-10to50',
    "/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM",
    "/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20_ext1-v1/NANOAODSIM",
    store=storage,url=url,file=filelist,opts='zpt=True',
  ),
  
  # TTBAR
  M('TT','TTTo2L2Nu',
    "/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM",
    store=storage,url=url,file=filelist,opts='toppt=True',
  ),
  M('TT','TTToSemiLeptonic',
    "/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM",
    "/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20_ext3-v1/NANOAODSIM",
    store=storage,url=url,file=filelist,opts='toppt=True',
  ),
  M('TT','TTToHadronic',
    "/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v3/NANOAODSIM",
    "/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20_ext2-v1/NANOAODSIM",
    store=storage,url=url,file=filelist,opts='toppt=True',
  ),
  
  # W+JETS
  M('WJ','WJetsToLNu',
    "/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM",
    store=storage,url=url,file=filelist,
  ),
  M('WJ','W1JetsToLNu',
    "/W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM",
    store=storage,url=url,file=filelist,
  ),
  M('WJ','W2JetsToLNu',
    "/W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM",
    store=storage,url=url,file=filelist,
  ),
  M('WJ','W3JetsToLNu',
    "/W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM",
    store=storage,url=url,file=filelist,
  ),
  M('WJ','W4JetsToLNu',
    "/W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM",
    store=storage,url=url,file=filelist,
  ),
  
  # SINGLE TOP
  M('ST','ST_tW_antitop',
    "/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20_ext1-v1/NANOAODSIM",
    store=storage,url=url,file=filelist,
  ),
  M('ST','ST_tW_top',
    "/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20_ext1-v1/NANOAODSIM",
    store=storage,url=url,file=filelist,
  ),
  M('ST','ST_t-channel_antitop',
    "/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM",
    store=storage,url=url,file=filelist,
  ),
  M('ST','ST_t-channel_top',
    "/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM",
    store=storage,url=url,file=filelist,
  ),
  
  # DIBOSON
  M('VV','WW',
    "/WW_TuneCP5_13TeV-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM",
    #"/WW_TuneCP5_PSweights_13TeV-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM",
    store=storage,url=url,file=filelist,
  ),
  M('VV','WZ',
    "/WZ_TuneCP5_13TeV-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM",
    #"/WZ_TuneCP5_PSweights_13TeV-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM",
    store=storage,url=url,file=filelist,
  ),
  M('VV','ZZ',
    "/ZZ_TuneCP5_13TeV-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM",
    store=storage,url=url,file=filelist,
  ),
  
  # SINGLE MUON
  D('Data','SingleMuon_Run2018A',"/SingleMuon/Run2018A-Nano25Oct2019-v1/NANOAOD",
   store=storage,url=url,file=filelist,channels=["skim*",'mutau','mumu','emu'],
  ),
  D('Data','SingleMuon_Run2018B',"/SingleMuon/Run2018B-Nano25Oct2019-v1/NANOAOD",
    store=storage,url=url,file=filelist,channels=["skim*",'mutau','mumu','emu'],
  ),
  D('Data','SingleMuon_Run2018C',"/SingleMuon/Run2018C-Nano25Oct2019-v1/NANOAOD",
   store=storage,url=url,file=filelist,channels=["skim*",'mutau','mumu','emu'],
  ),
  D('Data','SingleMuon_Run2018D',"/SingleMuon/Run2018D-Nano25Oct2019-v1/NANOAOD", # /SingleMuon/Run2018D-Nano25Oct2019_ver2-v1/NANOAOD ???
   store=storage,url=url,file=filelist,channels=["skim*",'mutau','mumu','emu'],
  ),
  
  # SINGLE ELECTRON
  D('Data','EGamma_Run2017A',"/EGamma/Run2018A-Nano25Oct2019-v1/NANOAOD",
    store=storage,url=url,file=filelist,channels=["skim*",'etau','ee'],
  ),
  D('Data','EGamma_Run2017B',"/EGamma/Run2018B-Nano25Oct2019-v1/NANOAOD",
    store=storage,url=url,file=filelist,channels=["skim*",'etau','ee'],
  ),
  D('Data','EGamma_Run2017C',"/EGamma/Run2018C-Nano25Oct2019-v1/NANOAOD",
    store=storage,url=url,file=filelist,channels=["skim*",'etau','ee'],
  ),
  D('Data','EGamma_Run2017D',"/EGamma/Run2018D-Nano25Oct2019-v1/NANOAOD",
    store=storage,url=url,file=filelist,channels=["skim*",'etau','ee'],
  ),
  
  # TAU
  D('Data','Tau_Run2017A',"/Tau/Run2018A-Nano25Oct2019-v1/NANOAOD",
    store=storage,url=url,file=filelist,channels=["skim*",'tautau'],
  ),
  D('Data','Tau_Run2017B',"/Tau/Run2018B-Nano25Oct2019-v1/NANOAOD",
    store=storage,url=url,file=filelist,channels=["skim*",'tautau'],
  ),
  D('Data','Tau_Run2017C',"/Tau/Run2018C-Nano25Oct2019-v1/NANOAOD",
    store=storage,url=url,file=filelist,channels=["skim*",'tautau'],
  ),
  D('Data','Tau_Run2017D',"/Tau/Run2018D-Nano25Oct2019_ver2-v1/NANOAOD",
    store=storage,url=url,file=filelist,channels=["skim*",'tautau'],
  ),
  
  ## LQ
  #M('LQ','SLQ_single_M1100_L1p0_old',
  #  "/LQ_Single_M1000/LegacyRun2_2018_deepTauIDv2p1/USER",
  #  store=storage,url=url,file=filelist,nfilesperjob=30,
  #),
  #M('LQ','SLQ_single_M1100_L1p0',
  #  "/SingleScalarLQ_InclusiveDecay_M-1100_L-1p0_TuneCP2_13TeV-madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM",
  #  #store=storage,url=url,file=filelist,
  #),
  #M('LQ','SingleVectorLQ_InclusiveDecay_M-800_L-1p5_K-0',
  #  "/SingleVectorLQ_InclusiveDecay_M-800_L-1p5_K-0_TuneCP2_13TeV-madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM",
  #  #store=storage,url=url,file=filelist,
  #),
  
]
