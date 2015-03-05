'''
    Generate tRNA mass and corresponding sequence
'''
import re
import csv

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def MassListParser(infile):
    # take mass list and return dictionary with seq as key and mass as value
    mass_dict = dict()
    input_file = csv.DictReader(open(infile), quotechar='|', quoting=csv.QUOTE_MINIMAL)
    for row in input_file:
        for eachkey in row:
            if is_number(row[eachkey]) and eachkey != "symbol":
                row[eachkey] = float(row[eachkey])
        mass_dict[row["symbol"]] = row
    return mass_dict

def GetMassForBase(base, mass_dict):
    if len(base) == 1:
        return mass_dict[base]["N"]
    elif base.endswith("p"):
        return mass_dict[base[0]]["Np"]

def GetMassForSeq(seq, mass_dict):
    total = 0
    seg = re.findall(".[p]?", seq)
    for each in seg:
        #print each,GetMassForBase(each, mass_dict)
        total += GetMassForBase(each, mass_dict)
    total += (len(seq.replace("p","")) - 1) * 63.98050
    return total

def FivePrimePho(seq):
    # five prime phosphate
    return seq + "p"

def ThreePrimePho(seq):
    # three prime phosphate
    return "p" + seq

def RNAGen(seq):
    # take tRNA sequence as input and generate all unique segment aftersplicing
    # RNase T1 digest
    seg = seq.split("G")
    #seg = re.findall(seq, "[^G]*G")
    # the last element in list is not spliced by digest protein
    seg_mod = [x + "Gp" for x in seg[:-1]]
    seg_mod.append(seg[-1])
    final_seg = list(set(seg_mod))
    final_seg.sort(key=len)
    return final_seg

def MassList(seq, mass_dict):
    seq_list = RNAGen(seq)
    masslist = []
    for each in seq_list:
        masslist.append([each, GetMassForSeq(each, mass_dict)])
    print masslist

def Test():
    massdict =  MassListParser("modified_bases_list.txt_mod")
    #GetMassForSeq("GGGp", massdict)
    print GetMassForSeq("CD", massdict)
    #assert(RNAGen("GGGGCUAUAGCUCAGCD") == ["CD", "Gp" , "CUCAGp", "CUAUAGp"])
    MassList("GGGGCUAUAGCUCAGCD", massdict)

if __name__ == "__main__":
    Test()
    RNAGen("GGGGCUAUAGCUCAGCDGGGAGAGCGCCUGCUUVGCACGCAGGAG")
