##### ATTRIBUTION #####


# Original Author:  Lakshay Sethi

'''Comments
    getcwd function of os package was used to retrieve the working directory.
    argv function of sys package was used to uilize the parameters passed with the file.
'''
import os
import sys
cwd = os.getcwd()

fileIN = open( cwd + "/results/battenberg-1.1/00-inputs/reference/" + sys.argv[1] + "/impute_info.txt", 'r')
filedata = fileIN.read()
fileIN.close()

newdata = filedata.replace("<path_to_impute_reference_files>", cwd + "/results/battenberg-1.1/00-inputs/reference/" + sys.argv[1] + "/battenberg_impute_v3")

fileOut = open(cwd + "/results/battenberg-1.1/00-inputs/reference/" + sys.argv[1] + "/impute_info.txt", 'w')
fileOut.write(newdata)
fileOut.close()

