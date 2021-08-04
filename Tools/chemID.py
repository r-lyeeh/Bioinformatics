# -*- coding:utf-8 -*-
import urllib
import urllib2
import json
import re
import requests
from bs4 import BeautifulSoup
import lxml
#headers={"User-Agent": "Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/51.0.2704.103 Safari/537.36"}
#request=urllib2.Request(url,header=headers)
#response=urllib2.urlopen(request)
#cookies={''}
url='https://chem.nlm.nih.gov/chemidplus/ProxyServlet?objectHandle=Search&actionHandle=searchChemIdLite&nextPage=jsp%2Fchemidheavy%2FChemidDataview.jsp&responseHandle=JSP&QV8=LD50&QF8=ToxTestType&QV6=mouse&QF6=ToxSpecies&QV7=intraperitoneal&QF7=ToxRoute&QV9=LUNGS%2C+THORAX%2C+OR+RESPIRATION&QF9=ToxEffect'
chem=[]

for i in range(1,99):
    j=1+25*i
    payload={'DT_START_ROW':j,'DT_ROWS_PER_PAGE':25}
    html=requests.get(url,params=payload)
    res=html.text
    soup=BeautifulSoup(res,"lxml")
    chemlist = soup.find_all('span',attrs={'class':"chem-name"})
    print "page"+str(i)
    for i in range(1,len(chemlist)):
        chem.append(chemlist[i].string)
print "spider network finished"
f = open ("chem-name.txt",'w')
for i in chem:
    f.write(i+"\n")
f.close()
