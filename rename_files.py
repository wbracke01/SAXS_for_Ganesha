import os
import csv
import tkinter as tk
from tkinter import filedialog
import re
root = tk.Tk()
root.withdraw()
folder = tk.filedialog.askdirectory()

files=os.listdir(folder)

for file in files:
    print(file)
    if file.endswith('.csv'):
        RunList=file
        RunDict=csv.DictReader(open(folder+'/'+RunList))


#
filenum=re.compile(r"\d\d\d\d\d\d\d")
batch_info_dicts=[]

for dict in RunDict:
    print(dict)
    batch_info_dicts.append(dict)

for file in files:
    if file.endswith('.tiff'):
        match=filenum.search(file)
        ID=match.group()
        for dict in batch_info_dicts:
            if ID in dict['Filename']:
               newName = dict['Sample Description'][12:] + '.tiff'
               os.rename(folder+'/'+file, folder+'/'+newName)
               print(dict['Filename']+' renamed to '+newName)

            else:
                pass


