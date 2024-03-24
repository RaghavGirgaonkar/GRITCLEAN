#!/usr/bin/env python
# coding: utf-8

# In[3]:


import requests
import os
import wget
import json
import h5py
from urllib.parse import urlparse
from pathlib import Path


# In[4]:


#Get all json filenames in specified directory
path_to_json = './'
json_files = [pos_json for pos_json in os.listdir(path_to_json) if pos_json.endswith('.json')]
print(json_files) 


# In[7]:


urlfile = open('urls.txt','a')
for file in json_files:
    f = open(file)
    jsonObj = json.load(f)
    strains = jsonObj['strain']
    for i in range(len(strains)):
        if strains[i]['format'] == 'hdf5' and strains[i]['duty_cycle'] == 100:
            fileURL = strains[i]['url']
            urlfile.write(fileURL + '\n')
            print(fileURL)
#             path = '/Users/kqa493/Desktop/GLITCH/IntervalFiles/LargeIntervalFiles/'
#             respFile = wget.download(fileURL, out = path)


        else:
            continue
    f.close()
    
urlfile.close()


# In[ ]:





# In[ ]:




