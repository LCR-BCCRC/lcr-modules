import os
import sys
cwd = os.getcwd()
cwd
fileIN = open( cwd + "/results/battenberg-1.1/00-inputs/reference/" + sys.argv[0] + "/impute_info.txt", 'r')
filedata = fileIN.read()
fileIN.close()

newdata = filedata.replace("<path_to_impute_reference_files>", cwd + "/results/battenberg-1.1/00-inputs/reference/" + sys.argv[0] + "/battenberg_impute_v3")

fileOut = open(cwd + "/results/battenberg-1.1/00-inputs/reference/" + sys.argv[0] + "/impute_info.txt", 'w')
fileOut.write(newdata)
fileOut.close()

