import ROOT
import argparse
import ctypes
import pandas as pd 
import re
# Open the ROOT file

def SampleName(HistName):
    process_name = re.search(r'.+_(.+)',HistName)
    return process_name.group(1) if process_name else HistName

def GetTaggerName(filepath):
    tagger_name = re.search(r'.+-(.+)-(20).+',filepath)
    if tagger_name:
        return tagger_name.group(1)
    else:
        print("Error: Unable to extract tagger name from the file path")
        exit()

def GetEraName(filepath):
    era_name = re.search(r'.+-(20.+).root',filepath)
    if era_name:
        return era_name.group(1)
    else:
        print("Error: Unable to extract era name from the file path")
        exit()


def getNEvents_mvis(filepath):
    
    file = ROOT.TFile.Open(filepath)

    tagger_dic = {}
    # needs to be tagger = {'QCD': value, 'TTbar': value, ...}


    # Check if the file is successfully opened
    if not file or file.IsZombie():
        print("Error: Unable to open ROOT file")
        exit()

    # Get the TCanvas saved in the ROOT file
    canvas = file.Get("canvas")

    # Check if the canvas is found
    if not canvas:
        print("Error: Canvas not found in the ROOT file")
        file.Close()
        exit()

    # Get the list of primitives (histograms, functions, etc.) drawn on the canvas
    primitives = canvas.GetListOfPrimitives()

    # Loop through the list to access each histogram
    for obj in primitives:
        #print(obj.GetName())
        print("first loop")
        subprims = obj.GetListOfPrimitives()
        for subprim in subprims:
            # print("subprim name",subprim.GetName())
            if subprim.GetName() == "stack_mt_1":
                histograms = subprim.GetHists()
                for hist in histograms:
                    print(">>>Histogram_Name:", hist.GetName())
                    #print("Integral: ", hist.Integral())
                    I_error=ctypes.c_double()
                    I_value=hist.IntegralAndError(0,hist.GetNbinsX()+1,I_error,"")
                    print("IntegralError: ",I_value, "±" ,I_error.value )

                    tagger_dic[SampleName(hist.GetName())] = (I_value, I_error.value)

                
            if isinstance(subprim, ROOT.TH1): 
                print("--------")
                print("Histogram_Name:", subprim.GetName())
                #print("Integral: ", subprim.Integral())
                I_error=ctypes.c_double()
                I_value=subprim.IntegralAndError(0,subprim.GetNbinsX()+1,I_error,"")
                print("IntegralError: ",I_value, "±" ,I_error.value )

                if 'hframe' in subprim.GetName(): continue
                if 'ratio_mt_1' in subprim.GetName():
                    tagger_dic[subprim.GetName()] = (I_value, I_error.value)
                    continue
                tagger_dic[SampleName(subprim.GetName())] = (I_value, I_error.value)

                # Access other properties/methods of the histogram as needed

    # Close the ROOT file
    file.Close()
    return tagger_dic

def main():
    parser = argparse.ArgumentParser(prog='openCanvas',
                                 description='calculate the integral of the histogram')

    parser.add_argument('--path',
                        action='extend',
                        nargs="+",
                        type=str,
                        required=True,
                        help='absolute path to the root file'
                        )
    parser.add_argument('-n','--name',
                        type=str,
                        help="Name for the table",
                        default="test",
                        # required=True
                        )

    args=parser.parse_args()
    paths = args.path


    table = {}

    region_name = []
    for path in paths:
        era_name = GetEraName(path)
        region_name.append(GetTaggerName(path))

        if len(region_name) > 1:
            if region_name[-2] !=  region_name[-1]:
                exit("Error: Different regions found in the provided file paths.")

        # if era_name not in eras:
        #     exit(f"Error: Era name {era_name} extracted from the file path does not match any of the provided eras.")

        print("Processing file:", path)
        table[era_name]=getNEvents_mvis(path)
    
    print(table)

    df = pd.DataFrame.from_dict(table, orient='index')
    df.index.name = region_name[0]
    df.reset_index(inplace=True)

    def format_pm(value):
        if isinstance(value, tuple):
            return f"{value[0]:.2f} \\pm {value[1]:.2f}"
        return value
    
    df = df.applymap(format_pm)
    df_T = df.T
    print(df_T)

    for era, values in table.items():
        print("--------")
        print(f"Eras: {era}")
        integral_expected = 0
        integral_wjets = 0
        integral_nowjets = 0
        for process, (value, error) in values.items():
            # print(f"{process}: {value} ± {error}")
            if "WJ" in process: integral_wjets += value
            elif ('Muon' in process or 'Run20' in process) and 'ratio' not in process: integral_expected += value
            # elif 'QCD' in process: continue
            else: integral_nowjets += value

            factor = ( integral_expected - integral_nowjets ) / integral_wjets if integral_wjets !=0 else 1
        # print("Expected : ", integral_expected)
        # print("no Wjets: ", integral_nowjets )
        # print("W+jets): ", integral_wjets)
        print(f"Scale factor for W+jets in {era}: {factor:.4f}")
            

if __name__ == "__main__":
    main()