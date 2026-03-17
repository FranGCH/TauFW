from postfit_TES import drawpostfit
import yaml

def main(args):
    configs   = args.configs
    againstjet = args.againstjet
    againstelectron = args.againstelectron

    for config in configs:
        if not config.endswith(".yml"): # config = channel name
            config = "config/setup_%s.yml"%(config) # assume this file name pattern
        print(">>> Using configuration file: %s"%config)
        with open(config, 'r') as file:
            setup = yaml.safe_load(file)
        tag = setup.get('tag',"")

    for region in setup["regions"]:
        if region == "baseline": continue

        else:
            print(">>>   Region: %s"%(region))
            era = "%s" % args.year  ## Use the year from command line arguments
            # Define the parameters
            fname = './postfit_pt_less_region/againstjet_%s/againstelectron_%s/%s/PostFitShape_%s__mutau_%s.root' %(againstjet, againstelectron, era, era,region)
            # bin = 'DM0'  # This should match the bin name in your ROOT file
            procs = setup["processes"]  # Replace with the actual processes in your file
            # procs = ["ZTT","ZL","ZJ","W","VV","ST","TTT","TTL","TTJ","QCD","data_obs"]  # Replace with the actual processes in your file
            text = setup["regions"][region]["title"]
            print(">>>   Title: %s"%(text))


            # Call the function
            drawpostfit(fname, region, procs,
                         outdir='output_plots/jet_%s_ele_%s/'%(args.againstjet,args.againstelectron), pname='$FIT.png', ratio=True, era=era, text=text)
            if args.include_cr:
                fname = './postfit_pt_less_region/againstjet_%s/againstelectron_%s/%s/PostFitShape_%s__mutau_%s.root' %(againstjet, againstelectron, era, era,region)    
                
                procs = ['ZL', 'ZTT', 'ZJ', 'W','VV','ST', 'TT','QCD','data_obs']

                # Try CR directory names that include the DM region so each DM gets its own CR plot
                import ROOT
                froot = ROOT.TFile.Open(fname)
                cr_candidates = [f"{args.cr_name}_{region}", f"{region}_{args.cr_name}", args.cr_name]
                found = False
                for cr_region in cr_candidates:
                    # drawpostfit expects directories named "<region>_prefit" / "<region>_postfit"
                    if froot and not froot.IsZombie() and froot.Get(f"{cr_region}_postfit"):
                        pname_cr = f"{region}_{cr_region}_$FIT_CR.png"
                        drawpostfit(fname, cr_region, procs,
                                     outdir='output_plots/jet_%s_ele_%s/'%(args.againstjet,args.againstelectron), pname=pname_cr, ratio=True, era=era, text="Z#rightarrow#mu#mu CR")
                        found = True
                        break
                if not found:
                    # fallback: draw with the plain CR name (still include region in filename)
                    drawpostfit(fname, args.cr_name, procs,
                                 outdir='output_plots/jet_%s_ele_%s/'%(args.againstjet,args.againstelectron), pname=f"{region}_{args.cr_name}_$FIT_CR.png", ratio=True, era=era, text="Z#rightarrow#mu#mu CR")
                if froot:
                    froot.Close()
if __name__ == "__main__":
    from argparse import ArgumentParser, RawTextHelpFormatter
    description = """Simple plotting script for postfit plots"""

    parser = ArgumentParser(prog="plot",description=description,epilog="Good luck!")

    parser.add_argument('-c', '--config', '--channel',
                                         dest='configs', type=str, nargs='+', default=['config/setup_mutau.yml'], action='store',
                                         help="config file(s) containing channel setup for samples and selections, default=%(default)r" )
    parser.add_argument('--include-cr', dest='include_cr', action='store_true', default=False,
                                         help="also draw control-region (Zmm) prefit/postfit plots" )
    parser.add_argument('--cr-name', dest='cr_name', type=str, default='Zmm',
                                         help="control-region directory prefix in ROOT file (default='Zmm')" )
    parser.add_argument('-j', '--jet', dest='againstjet', default='Tight', help="against jet WP")
    parser.add_argument('-e', '--electron', dest='againstelectron', default='Tight', help="against electron WP")
    parser.add_argument('-y', '--year', dest='year', default='2024', help="year for plotting")

    args = parser.parse_args()
  
    main(args)
    print("\n>>> Done.")