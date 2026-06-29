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
    tagger_name = re.search(r'.+-(.+)_(region|hpscond).+',filepath)
    if tagger_name:
        return tagger_name.group(1)
    else:
        print("Error: Unable to extract tagger name from the file path")
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
            if subprim.GetName() == "stack_mvis":
                histograms = subprim.GetHists()
                for hist in histograms:
                    # print(">>>Histogram_Name:", hist.GetName())
                    #print("Integral: ", hist.Integral())
                    I_error=ctypes.c_double()
                    I_value=hist.IntegralAndError(0,hist.GetNbinsX()+1,I_error,"")
                    # print("IntegralError: ",I_value, "±" ,I_error.value )

                    tagger_dic[SampleName(hist.GetName())] = (I_value, I_error.value)

                
            if isinstance(subprim, ROOT.TH1): 
                # print("Histogram_Name:", subprim.GetName())
                #print("Integral: ", subprim.Integral())
                I_error=ctypes.c_double()
                I_value=subprim.IntegralAndError(0,subprim.GetNbinsX()+1,I_error,"")
                # print("IntegralError: ",I_value, "±" ,I_error.value )

                if 'hframe' in subprim.GetName(): continue
                if 'ratio_mvis' in subprim.GetName():
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
    parser.add_argument('--taggers',
                        action='store',
                        nargs='+',
                        type=str,
                        default=['DeepTau','PNet','UParT'],
                        help='list of taggers to process',
                        required=False
                        )
    parser.add_argument('-n','--name',
                        type=str,
                        help="Name for the table",
                        default="test",
                        )

    args=parser.parse_args()
    paths = args.path
    taggers = args.taggers

    tagger_table = {}

    if len(paths) != len(taggers):
        exit("Error: The number of paths and taggers must match.")

    for path, tagger in zip(paths, taggers):
        taggername_path = GetTaggerName(path)

        if taggername_path.lower() != tagger.lower():
            exit("Error: Tagger name in the file path does not match the provided tagger name.")

        print("Processing file:", path)
        tagger_table[tagger]=getNEvents_mvis(path)

    print(tagger_table)

    df = pd.DataFrame.from_dict(tagger_table, orient='index')
    df.index.name = 'Tagger'
    df.reset_index(inplace=True)

    def format_pm(value):
        if isinstance(value, tuple):
            return f"{value[0]:.2f} \\pm {value[1]:.2f}"
        return value
    
    df = df.applymap(format_pm)
    df_T = df.T
    print(df_T)
    latex = df_T.to_latex(index=True, column_format="l|c|c|c|" ,escape=False)  # escape=False keeps the \pm
    # print(latex)

    with open(f'table_output_{args.name}.tex', "w") as f:
        f.write(latex)
    
    ratio_upart_deepTau = tagger_table['UParT']['ZTT'][0] / tagger_table['DeepTau']['ZTT'][0]
    ratio_pnet_deepTau = tagger_table['PNet']['ZTT'][0] / tagger_table['DeepTau']['ZTT'][0]
    print(f"Ratio UParT/DeepTau: {ratio_upart_deepTau:.4f}")
    print(f"Ratio PNet/DeepTau: {ratio_pnet_deepTau:.4f}")

if __name__ == "__main__":
    main()