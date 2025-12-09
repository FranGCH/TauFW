from postfit_TES import drawpostfit
import yaml

def main(args):
    configs   = args.configs
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
            era = "2024" ## Hardcoded
            # Define the parameters
            fname = '/eos/home-f/fcasalin/TauFW_230425/Fitter_out/postfit_pt_less_region/againstjet_Medium/againstelectron_VVLoose/%s/PostFitShape_2024__mutau_%s.root' %(era,region)
            # bin = 'DM0'  # This should match the bin name in your ROOT file
            procs = setup["processes"]  # Replace with the actual processes in your file
            # procs = ["ZTT","ZL","ZJ","W","VV","ST","TTT","TTL","TTJ","QCD","data_obs"]  # Replace with the actual processes in your file
            text = setup["regions"][region]["title"]
            print(">>>   Title: %s"%(text))


            # Call the function
            drawpostfit(fname, region, procs,
                         outdir='/eos/home-f/fcasalin/TauFW_230425/Fitter_out/output_plots', pname='$FIT.png', ratio=True, era=era, text=text)
            if args.include_cr:
                fname = '/eos/home-f/fcasalin/TauFW_230425/Fitter_out/postfit_pt_less_region/againstjet_Medium/againstelectron_VVLoose/%s/PostFitShape_2024__mutau_%s.root' %(era,region)    
                
                procs = ['ZL', 'ZTT', 'ZJ', 'W','VV','ST', 'TT','QCD','data_obs']

                drawpostfit(fname, args.cr_name, procs,
                             outdir='/eos/home-f/fcasalin/TauFW_230425/Fitter_out/output_plots', pname=f"$FIT'+'_CR_{region}.png", ratio=True, era=era, text="Z#rightarrow#mu#mu CR")
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
     
    args = parser.parse_args()
  
    main(args)
    print("\n>>> Done.")