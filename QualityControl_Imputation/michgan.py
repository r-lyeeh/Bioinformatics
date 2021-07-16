#!/bin/bash
import requests
import json
import sys

#usage: michgan.py -d file_directory/filename_
folder=sys.argv[1]

# imputation server url
url = 'https://imputationserver.sph.umich.edu/api/v2'
token = 'eyJjdHkiOiJ0ZXh0XC9wbGFpbiIsImFsZyI6IkhTMjU2In0.eyJtYWlsIjoic2hpY2hhby5wYW5nQHR1bS5kZSIsImV4cGlyZSI6MTYyOTAyODUxOTE3MiwibmFtZSI6IlNoaWNoYW8gUGFuZyIsImFwaSI6dHJ1ZSwidXNlcm5hbWUiOiJtb3N0aW1hIn0.EzljJXNT6RKT1mzGXMTqYPu2095EnROljhRjWu9Hd8c';

# add token to header (see Authentication)
headers = {'X-Auth-Token' : token }
#mode=imputation
#phasing=eagle
data = {
  #'refpanel': '1000g-phase-3-v5',
  'refpanel': 'hrc-r1.1',
  'population': 'eur',
  'r2Filter': '0.3'
}

# merge string
def str_join(*args):
    return ''.join(map(str, args))

# submit new job
vcf = str_join(folder,1,'.recode.vcf.gz');
for i in range(2,23):
    # calculate value
    globals()["vcf"+str(i-1)] = str_join(folder,i,'.recode.vcf.gz')

files = [('files',open(vcf,'rb'))]
for i in range(2,23):
    files.append(('files',open(globals()["vcf"+str(i-1)],'rb')))
#files = [('files', open(vcf, 'rb')), ('files', open(vcf1, 'rb')),('files', open(vcf2, 'rb')), ('files', open(vcf3, 'rb')),('files', open(vcf4, 'rb')), ('files', open(vcf5, 'rb')),('files', open(vcf6, 'rb')), ('files', open(vcf7, 'rb')),('files', open(vcf8, 'rb')), ('files', open(vcf9, 'rb')),('files', open(vcf10, 'rb')), ('files', open(vcf11, 'rb')),('files', open(vcf12, 'rb')), ('files', open(vcf13, 'rb')),('files', open(vcf14, 'rb')),('files', open(vcf15, 'rb')), ('files', open(vcf16, 'rb')), ('files', open(vcf17, 'rb')),('files', open(vcf18, 'rb')), ('files', open(vcf19, 'rb')), ('files', open(vcf20, 'rb')),('files', open(vcf21, 'rb'))]

r = requests.post(url + "/jobs/submit/minimac4", files=files, data=data, headers=headers)
if r.status_code != 200:
  print(r.json()['message'])
  raise Exception('POST /jobs/submit/minimac4 {}'.format(r.status_code))

# print message
print(r.json()['message'])
print(r.json()['id'])
