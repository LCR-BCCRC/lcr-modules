import sys
import re
import pandas as pd

keys_to_check = ["SVTYPE", "SVLEN", "END","CIPOS","CIEND","CIGAR","EVENT", "LEFT_SVINSSEQ","RIGHT_SVINSSEQ","BND_DEPTH","MATE_BND_DEPTH","SOMATIC","SOMATICSCORE","JUNCTION_SOMATICSCORE", "IMPRECISE","MATEID", "HOMLEN", "HOMSEQ", "SVINSLEN","SVINSSEQ"]
vcf_array = []

#print(vcf_file)
def parse_vcf(vcf_file):
    ###Open your file
    with open(vcf_file, 'r') as vcf_f:
        for line in vcf_f:
            ###Skip metadata lines
            if line[0] != '#':
                ###Split the line by "tab" to keep info field (for me it's the 8th column, so choose the index = 7, I don't remember if it's always the 8th, change this if needed)
                info_field_line = line.split("\t")[7]
                ###Split the info line by ";"
                info_field_line_array = info_field_line.split(";")
                pairs = info_field_line.strip().split(";")
                ###Find missing keys from VCF field
                missing_keys = {key for key in keys_to_check if key not in {p.split('=')[0] for p in pairs}}
                #print(missing_keys)
                ###For each line of your VCF, create a dictionary with this array key : info, value : value of this info
                dict_info={}
                dict_info["chrom"] = line.split("\t")[0]
                dict_info["pos"] = line.split("\t")[1]
                dict_info["ID"] = line.split("\t")[2]
                dict_info["REF"] = line.split("\t")[3]
                dict_info["ALT"] = line.split("\t")[4]
                for i in info_field_line_array:
                    ###Just looking for line with "=" character (as key = value)
                    if "=" in i:
                        ###Left from equal sign is key 
                        key = i.split("=")[0]
                        ###Right from equal sign is value 
                        value = i.split("=")[1]
                        ###Put them in a dictionary
                        dict_info[key]=value
                    for key in missing_keys:
                        dict_info[key] = "NA"
                ###This is the result for the first line, you can save the data in array to use it later or process line by line as you wish
                #print(dict_info)
                vcf_array.append(dict_info.copy())

    df = pd.DataFrame(vcf_array)
    df.to_csv(snakemake.output["tsv"], sep='\t')

if __name__ == '__main__':
    parse_vcf(snakemake.input["vcf"])