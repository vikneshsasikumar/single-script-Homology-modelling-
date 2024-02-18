from modeller import *
from modeller import *
from modeller.automodel import *
import webbrowser as wb
#from modeller import soap_protein_od
import os
import requests

#pdb file download from PDB database

pdb=input("ENTER YOUR TEMPLATE PDB ID IN CAPS : ")

# URL of the file to download
url = "https://files.rcsb.org/download/"+pdb+".pdb"

# Send a GET request to the URL
response = requests.get(url)

# Check if the request was successful
if response.status_code == 200:
    # Get the file name from the URL
    file_name = url.split('/')[-1]
    # Save the file to the current working directory
    file_path = os.path.join(os.getcwd(), file_name)
    with open(file_path, 'wb') as f:
        f.write(response.content)
    print(f"File downloaded to {file_path}")
else:
    print("Failed to download file")

#query alignment file creating for template and query alignment
    
que_name=input("MAKE YOUR QUERY ALIGNMENT FILE NAME : ")
f= open(que_name+".ali","w+")
sequence=input("""ENTER YOUR QUERY SEQUENCE WITHOUT HEADER FILE: """)
with open(que_name+".ali","w") as file:
    file.write(">P1;"+que_name+"\n")
    file.write("sequence:"+que_name+":::::::0.00: 0.00\n")
    file.write(sequence+"*")

template=pdb.lower()
temp=template+".pdb"
que=que_name+".ali"
#Homology model building 
num=int(input("HOW MANY MODELS YOU WANT TO DO : "))

env = environ()
aln = alignment(env)
mdl = model(env, file=template, model_segment=('FIRST:A','LAST:A'))
aln.append_model(mdl, align_codes=template, atom_files=temp)
aln.append(file=que, align_codes=que_name)
aln.align2d()
aln.write(file=que_name+"-"+template+'.ali', alignment_format='PIR')
aln.write(file=que_name+"-"+template+'.pap', alignment_format='PAP')


env = environ()
a = automodel(env, alnfile=que_name+"-"+template+'.ali',
              knowns=template, sequence=que_name,
              assess_methods=(assess.DOPE,
                              #soap_protein_od.Scorer(),
                              assess.GA341))
a.starting_model = 1
a.ending_model = num
a.make()
