import os,sys
import json, re, time
from tkinter.tix import WINDOW
import urllib.request, urllib.parse, urllib.error
import traceback
from datetime import datetime
import pickle
import subprocess
from operator import itemgetter
import random
import math
import importlib
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

from xml.dom.minidom import parseString

import ssl
ssl._create_default_https_context = ssl._create_unverified_context

#AUTO-load classifiers
#a trick to get the current module
_modname = globals()['__name__']
_this_mod = sys.modules[_modname]

# _ppath = "/".join(_this_mod.__file__.split("/")[:-1])
_ppath = os.getcwd()

import getGEOSamples_byType_gse
import scrna_parser_detail_gse
import pubmed
# ## build mysql environment
# import env
# models = env.models


""" ========== main script ========== """

def getSyncLog(infoStr):
    """ouput the record to DoneGsmXml.log file
    """
    os.system('echo "[%s] %s"' % (time.strftime('%H:%M:%S'), infoStr))


### GDS interface
def getGDSSamples(date_region=False):
    """Will run the predefined query and return a list of GDS ids
    NOTE: this returns ALL GDS samples which are of SRA type i.e.
    ALL CHIP-SEQ, RNA-SEQ, etc.
    """
    #expireDate = now - 30 days in seconds
    #ref: http://stackoverflow.com/questions/7430928/python-comparing-date-check-for-old-file
    # _expireDate = time.time() - 60 * 60 * 24 * 30

    ret = []
    # #TRY: to read a file first -- IF IT IS NOT STALE
    path = os.path.join(_ppath, "gdsSamples.txt")
    
    URL = """http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=SRA[Sample%20Type]%20AND%20gse[Entry%20Type]%20AND%20(homo%20sapiens[Organism]%20OR%20mus%20musculus[Organism])&retmax=100000&usehistory=y"""
    
    if date_region:
        maxTime = date_region.split('-')[1]
        minTime = date_region.split('-')[0]

        URL = """http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=SRA[Sample%20Type]%20AND%20gse[Entry%20Type]%20AND%20(homo%20sapiens[Organism]%20OR%20mus%20musculus[Organism])&mindate={0}&maxdate={1}&datetype=pdat&retmax=100000&usehistory=y""".format(minTime, maxTime)
    try:
        getSyncLog("getGDSSample: %s" % URL) # output record
        f = urllib.request.urlopen(URL)
        root = ET.fromstring(f.read())
        f.close()

        #Get the IDList
        tmp = root.findall("IdList/Id")
        ret = [i.text for i in tmp] # list of gds IDs

        #write to disk all the IDs - to gdsSamples.txt
        f = open(path, "w")
        for l in ret:
            f.write("%s\n" % l)
        f.close()
    except:
        print("Exception in user code:")
        print('-' * 60)
        traceback.print_exc(file=sys.stdout)
        print('-' * 60)
    return ret
def gse_idToAcc(gdsId):
    """Given a GDS id, e.g. 300982523, tries to give a GDS accession, e.g.
    GSM982523

    NOTE: there is an algorithm: acc = "GSM"+gdsId[1:] (strip leading 0s)
    """
    cut = gdsId[1:].lstrip("0")
    return "GSE%s" % cut

def proxyInstead(link, using=False):
    """using proxy to avoid forbidden
    """
    context = ''
    if using:
        #using proxy first
        try: # using proxy first, or using the read ip
            agent = [x.rstrip() for x in open('./pickle_file/proxy.txt')]
            proxy = {'http':'http://%s'%random.sample(agent, 1)[0]}
            urlf = urllib.request.urlopen(link, proxies = proxy)
            getSyncLog('.')
        except:
            urlf = urllib.request.urlopen(link)
            proxy = {'proxy':'local IP'}
            getSyncLog('.') # use for record, so that we can know what happened if error occured
        # check whether we get the correct inf
        context = urlf.read()
        urlf.close()
        if ('404 - File or directory not found' in context) or ('ERR_ACCESS_DENIED' in context) or (context.strip() == ''):
            urlf = urllib.request.urlopen(link)
            context = urlf.read()
            urlf.close()
            proxy = {'proxy':'local IP'}
            getSyncLog('.')
        getSyncLog('%s: %s'%(list(proxy.values())[0], link))
        context = context.decode(encoding='utf-8',errors='ignore')
        return context
    try:
        # time.sleep(0.3)
        urlf = urllib.request.urlopen(link)
        print("opening url")
        context = urlf.read()
        context = context.decode(encoding='utf-8',errors='ignore')
        urlf.close()
        getSyncLog('local IP: %s'%link) 
        return context
    except: 
        pass
        print("still did not pass")
        print('link problem: %s'%link)
    return None

def isXML(doc):
    """
    TEST if it is a valid geo XML record
    NOTE: first lines are-
    <?xml version="1.0" encoding="UTF-8" standalone="no"?>
    """
    
    if doc == None:
        return None
    f = doc.split("\n")
    
    return f[0].strip() == """<?xml version="1.0" encoding="UTF-8" standalone="no"?>"""

  
def getGeoXML(accession, path):
    """
    HANDLES GEO XML records--i.e. our GEO XML librarian!
    Given a GEO ACCESSION ID, return the xml record for it
    (making the urllib call)
    """
    #path pattern: EXAMPLE-GSE1126513 geo/GSE1126/GSE1126513
    #path = os.path.join(_ppath, ddir)

    # empty directory of xmls files not needed anymore 
    # cmd = "rm -rf " + path + "/*"
    # os.system(cmd)

    if not os.path.exists(path):
        os.mkdir(path)
    subdir = os.path.join(path, accession[:7])
    if not os.path.exists(subdir):
        os.mkdir(subdir)
    path = os.path.join(subdir, "%s.xml" % accession)
    if os.path.exists(path):
        print("opening path")
        f = open(path)
        docString = f.read()
        f.close()
    else:
        #target:family; view:full, form=xml
        URL = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&view=full&form=xml&targ=all" % accession
        try:
            #print "getGeoXML: %s" % URL
            #signal.alarm(180)
            docString = proxyInstead(link=URL)
            '''if not isXML(docString): # try again
                #signal.alarm(180)
                getSyncLog('.')
                docString = proxyInstead(link=URL)'''
            if isXML(docString):
                #write to file
                f = open(path, "w")
                f.write(docString)
                f.close()
                #getSyncLog(proxy.values()[0]+'\t'+accession + '\n')# output record
            if isXML(docString) == None:
            #else:
                print(accession)
                print("ERROR: accession is NOT xml. (The accession may be deleted from GEO repository)")
                f1 = open('gsm_notXML.txt', 'a')
                f1.write(accession + '\n')
                f1.close()
                return None
        except:
            print("Exception in user code:")
            print('-' * 60)
            traceback.print_exc(file=sys.stdout)
            print('-' * 60)
            docString = None
    return docString

#new parser
def _sync_gse(fsave, a_or_w, add_index_or_not,key_file_input,key_logic, xmlPath, DataType=False, dateRegion = False, refresh=False, exludeFile = False):
    ## get all GDS ids of GSE from API
    gdsSamples = getGDSSamples(dateRegion) # get all geo dataset (gds) id
    print(gdsSamples)
    getSyncLog('start: There are %s GDS Samples in sum'%(len(gdsSamples)))#
    exclude = [] # the gse maybe existed
    if exludeFile:
        exclude = [x.rstrip() for x in open(exludeFile)]
    # open a file to save information
    # if a_or_w == 'a':
        # f = open(fsave, 'a')
    # else:
        # f = open(fsave, 'w')
    if add_index_or_not=='add':
        fsave=fsave.split('.csv')[0]+'_'+str(dateRegion).replace('/','_')+'.txt'

    
    # if not os.path.exists('fsave'):
    #     f = open(str(fsave), "x")
    #     f.close()

    f = open(fsave, str(a_or_w)) 
    f.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%('GSE', 'Species',
        'PMID', 'Paper', 'Title', 'CellType', 'Tissue', 'Disease',
        'Cell_Pop', 'Release_Date', 'Last_Update_Date', 'SRA', 'Cell_Line', 'GSMs','fields_dataType', 'Key_Match')+'\n')
    f.close()
    # if a_or_w == 'a':
        # out = open(fsave+'_others.txt', 'a')
    out = open(fsave+'_others.txt', str(a_or_w))
    out.write('')
    out.close()
    # convert to gds id to GSE and download XML file
    cnt = 0
    # one_percent = int(len(gdsSamples)/100) #
    for gds in gdsSamples:
        cnt += 1
        # print(cnt,'#',one_percent)
        if cnt % math.ceil(len(gdsSamples)/100) == 0:
            getSyncLog("%s%%"%(100*cnt/len(gdsSamples)))
            time.sleep(3)
        
        gseid = gse_idToAcc(gds)
        print(gseid)
        if gseid in exclude: # if existed, don't parse again
            continue

        gseXML = getGeoXML(gseid, xmlPath)

        getType = getGEOSamples_byType_gse.getGeoSamples_byTypes(path = "repository_samples.pickle", datatype = ['crispr'],
                                                                 gseids=[gseid], refresh=refresh, ddir = xmlPath,key_file_input=key_file_input,key_logic=key_logic,xml_path=xmlPath)

        if getType:
            for i in getType.keys():
                # parse sample annotation
                list_sample = scrna_parser_detail_gse.update_one_sample(gseid=gseid, ddir=xmlPath)
                try:
                    list_sample.append(str(getType[i])) # add matched key words
                    print("writing to file")
                    f = open(fsave, 'a')
                    f.write('\t'.join(list_sample)+'\n')
                    f.close()
                except:
                    getSyncLog("Error when writing in table: %s" % s)
        else:
            out = open(fsave+'_others.txt', 'a')
            out.write(''.join(gseid)+'\n')
            out.close()
        # out = open('geo_gse_collection.txt', 'a')
        # out.write(gseid+'\n')
        # out.close() 
        # time.sleep(0.03) # sleep to avoid IP blocking

    getSyncLog('done!')#



if __name__ == '__main__':
    _sync_gse()



