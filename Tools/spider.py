# -*-coding:utf8 -*-
import cookielib
import urllib2
import re
import gzip
#import pandas as pd
import time
from StringIO import StringIO

#cookie = cookielib.CookieJar()
#handler = urllib2.HTTPCookieProcessor(cookie)
#opener = urllib2.build_opener(handler)
#urllib2.install_opener(opener)
user_agent = 'Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/49.0.2623.110 Safari/537.36'
headers = {'Accept':'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
           'Accept-Encoding':'gzip, deflate, sdch',
           'Accept-Language':'zh-CN,zh;q=0.8,ja;q=0.6',
           'Content-Encoding':'gzip',
           'Connection':'keep - alive',
           'Cookie':'test=cookies; PHPSESSID=eobdmmlkhv4alcphdn2irp81l2;__utmc=124367126; visit=2',
           'Host':'metlin.scripps.edu',
           'User-Agent':user_agent}


  	
def getdisease(ddclass,disease):
		url = "http://www.malacards.org/categories/"+ddclass
    request = urllib2.Request(url, headers=headers)
    time.sleep(10)
    response = urllib2.urlopen(request)
    buf = StringIO(response.read())
    f = gzip.GzipFile(fileobj=buf)
    content = f.read()
    pattern1 = re.compile('<a href="/card/(.*)".*</a>', re.S)
    dise = re.findall(pattern1, content)
    disease.extend(dise)
    #print disea
    return disease
    
def getgene(ddclass,gene):
		disease=getdisease(ddclass,disease)
		for i in disease:
			url = "http://www.malacards.org/card/"+i
    	request = urllib2.Request(url, headers=headers)
    	time.sleep(10)
    	response = urllib2.urlopen(request)
   	 	buf = StringIO(response.read())
   	 	f = gzip.GzipFile(fileobj=buf)
   	 	content = f.read()
    	pattern1 = re.compile('<a href="http://genecards.org/cgi-bin/carddisp.pl?gene=(.*)" target="_blank".*</a>', re.S)
    	qe = re.findall(pattern1, content)
    	ge = disease,ge
    	gene.extend(qe)
    	print ge
    	return gene    
    
#def mapp(q,i):
#    mola = []
#    mapp0 = q.iat[i, 0], getmolecule(q.iat[i, 1], mola)
#    return mapp0


#for i in range(len(q)):
#    tota.append(mapp(q, i))
#tota.append(mapp(q,0))
#tota.append(mapp(q,1))
#print tota


url = "http://www.malacards.org/categories"
request = urllib2.Request(url, headers=headers)
time.sleep(10)
response = urllib2.urlopen(request)
buf = StringIO(response.read())
f = gzip.GzipFile(fileobj=buf)
content = f.read()
pattern1 = re.compile('<.*class>(.*)</a>', re.S)
dclass = re.findall(pattern1, content)
ddclass.extend(dclass)   
print dclass

#for i in ddclass:
#	gene.extend(getgene(ddclass,gene)


