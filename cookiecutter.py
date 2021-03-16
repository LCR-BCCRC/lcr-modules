import pyfiglet
import os
import subprocess

inputs = []


def welcome(text):
    return pyfiglet.figlet_format(text, font="small")


def modName(inputs):
    while True:
        os.system("clear")
        print(welcome("Enter Module Name:"))
        print(
            """
        Important:
        This field should only consist of lowercase alphanumerical characters or underscores (i.e. no spaces).
        """
        )
        choice = input("\n=> ")
        inputs.append(choice)
        if choice != "1":
            modAuthor(inputs)
        os.system("clear")


def modAuthor(inputs):
    while True:
        os.system("clear")
        print(welcome("Enter Module Author:"))
        choice = input("\n=> ")
        inputs.append(choice)
        if choice != "1":
            origAuthor(inputs)
        os.system("clear")


def origAuthor(inputs):
    while True:
        os.system("clear")
        print(welcome("Enter Original Author:"))
        print(
            """
        Important:
        This field should contain the  full name of the person who originally wrote the Snakefile or script that is 
        being used to inform the module.If the module is being written from scratch, this field can be set to N/A.
        """
        )
        choice = input("\n=> ")
        inputs.append(choice)
        if choice != "1":
            input_Output_File(inputs)
        os.system("clear")


def input_Output_File(inputs):
    while True:
        os.system("clear")
        print(welcome("Enter Input File Type:"))
        print(
            """
        Important:
        -> If there is more than one input file type, just list one of them for now. The same applies for the output file type. 
        -> Each of these should only consist of lowercase alphanumerical characters or underscores (i.e. no spaces).
        """
        )
        choice1 = input("\n=> ")
        inputs.append(choice1)
        os.system("clear")
        print(welcome("Enter Output File Type:"))
        print(
            """
        Important:
        -> If there is more than one input file type, just list one of them for now. The same applies for the output file type. 
        -> Each of these should only consist of lowercase alphanumerical characters or underscores (i.e. no spaces).
        """
        )
        choice2 = input("\n=> ")
        inputs.append(choice2)
        if choice1 != " ":
            modRun(inputs)
        os.system("clear")


def modRun(inputs):
    while True:
        os.system("clear")
        print(welcome("Enter Module run per"))
        print("\n1 - sample ")
        print("2 - tumour ")
        print("""Choose from 1, 2 : """)
        choice = input("\n=> ")
        inputs.append(choice)
        if choice == "1":
            seqType(inputs)
        elif choice == "2":
            seqType(inputs)
        else:
            exit()

        os.system("clear")


def seqType(inputs):
    while True:
        os.system("clear")
        seqType = ["Genome", "Capture", "MRNA"]
        print(welcome("Enter Sequence Type"))
        for x in seqType:
            print("\n" + x + ":")
            if inputs[5] == "1":
                print("\n1 - omit ")
                print("2 - unpaired")
                print("""Choose from 1, 2 : """)
                choice = input("\n=> ")
                if choice == "1":
                    inputs.append(choice)
                elif choice == "2":
                    inputs.append(choice)
                else:
                    exit()
            elif inputs[5] == "2":
                print("\n1 - matched_only  ")
                print("2 - allow_unmatched ")
                print("3 - no_normal ")
                print("""Choose from 1, 2, 3 : """)
                choice = input("\n=> ")
                if choice == "1":
                    inputs.append("3")
                elif choice == "2":
                    inputs.append("4")
                elif choice == "3":
                    inputs.append("5")
                else:
                    exit()
        data = ""
        for i in range(9):
            data = data + inputs[i] + "\n"
        subprocess.run(
            ["cookiecutter", "template/", "--output-dir", "modules/"],
            stdout=subprocess.PIPE,
            text=True,
            input=data,
        )

        exit()
        os.system("clear")


def main():
    while True:
        os.system("clear")
        print(welcome("LCR-Modules"))
        print("\nWelcome to the Module builder! ")
        print(
            "\nFor furthure info please visit : https://lcr-modules.readthedocs.io/en/latest/for_developers.html"
        )
        print(
            """
        Press 1 to start"""
        )
        choice = input("\nEnter your choice : ")
        if choice == "1":
            modName(inputs)
        else:
            exit()
        os.system("clear")


if __name__ == "__main__":
    main()

