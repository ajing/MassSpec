'''
    Integrate tRNA unique mass and find corresponding peaks
'''

from RNAMassGen import GetUniqueMass
from GetPeaks import GetPeaksWithRTEstimate

def main():
    ms_file   = "./E165ug.mzML"
    mass_dict = GetUniqueMass("./Data/tRNAseq.txt")
    print mass_dict
    mass_list = []
    mass_all_possible = dict()
    for rna in mass_dict:
        mass_all_possible[rna] = dict()
        for seg_mass in mass_dict[rna]:
            all_mass = AllMassPossibility(seg_mass[1])
            mass_all_possible[rna][seg_mass[0]] = all_mass
            mass_list += all_mass
    print mass_all_possible
    max_int_dict = GetPeaksWithRTEstimate(ms_file, mass_list)
    print max_int_dict
    for header in mass_all_possible:
        print header
        for seg in mass_all_possible[header]:
            print seg
            for mass in mass_all_possible[header][seg]:
                print mass
                print "max intensity:", max_int_dict[mass]


if __name__ == "__main__":
    main()
