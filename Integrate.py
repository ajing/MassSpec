'''
    Integrate tRNA unique mass and find corresponding peaks
'''

from RNAMassGen import GetUniqueMass
from GetPeak import GetPeakbyMZRange

def main():
    ms_file   = "./Data"
    mass_dict = GetUniqueMass("./Data/tRNAseq.txt")
    mass_list = []
    for rna in mass_dict:
        mass_list += [x[1] for x in mass_dict[rna]]
    max_int_dict = GetPeakbyMZRange(ms_file, mass_list)
    for header in mass_dict:
        print header
        for mass in mass_dict[header]:
            print mass, "max intensity:", max_int_dict[mass]


if __name__ == "__main__":
    main()
