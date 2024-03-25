import os,sys
import json, re, time
import urllib.request, urllib.parse, urllib.error
import traceback
from datetime import datetime
import pickle
import subprocess
from operator import itemgetter
import random
import importlib
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

from xml.dom.minidom import parseString
import os.path
from os import path
import pandas as pd
import re
import numpy as np

##
#AUTO-load classifiers
#a trick to get the current module
_modname = globals()['__name__']
_this_mod = sys.modules[_modname]

_ppath = "/".join(_this_mod.__file__.split("/")[:-1])

# from django.utils.encoding import smart_str
import env
# models = env.models

import scrna_parser_from_gse

def readGeoXML(path,xml_path,docString=None):
    """
    Input: a file path or a string--default is to use the path

    Tries to read in the geo xml record,
    **KEY: REMOVES the xmlns line
    Returns the xml record text WITHOUT the xmlns line!!!
    """
    if docString:
        f = docString.split("\n")
    else:
        
        # print("readGEOXML_printing path:", path) 
        #f = ''
        if path:
            if os.path.exists(path):
                f = open(path)
            else:
                return None    
        else:
            f = ''
    tmp = []
    try:
        for l in f:
            if l.find("xmlns=\"http://www.ncbi.nlm.nih.gov/geo/info/MINiML\"") == -1:
                tmp.append(l)
    except: # sometimes, the xml is not saved as utf-8 in local, do getGeoXml again
        gsmid = os.path.basename(path).rstrip('.xml')
        ag = scrna_parser_from_gse.getGeoXML(accession=gsmid,xml_path=xml_path)
        f = open(os.path.join(xml_path+gsmid[:7]+'/'+gsmid+'.xml'))
        # f = open(os.path.join('./geo/'+gsmid[:7]+'/'+gsmid+'.xml'))
        tmp = []
        for l in f:
            if l.find("xmlns=\"http://www.ncbi.nlm.nih.gov/geo/info/MINiML\"") == -1:
                tmp.append(l)
    if not docString and not isinstance(f, str):
        f.close()
    return "".join(tmp)


def _getFieldXML(sample_path,xml_path, fields = ['Sample/Library-Strategy', 'Sample/Description', 'Sample/Data-Processing',
    'Sample/Channel/Extract-Protocol', "Sample/Title", 'Sample/Channel/Source', 'Sample/Channel/Characteristics', 'Sample/Relation' ]):
    """
    get text of items in fields list from XML
    we need these info for matching key words
    """
    text = readGeoXML(sample_path,xml_path)
    ## parse XML content
    rt = ET.fromstring(text) #<Element 'MINiML' at 0x7fb5bb37f2f0> 

    info = {}
    field = 'Sample/Relation'
    tmp = rt.findall(field)
    

    for field in fields:
        tmp = rt.findall(field)
        if field == 'Sample/Relation':
            if tmp:
                srx_list = []
                for x in tmp: 
                    if x.attrib['type'] == 'SRA':
                        url = x.attrib['target']
                        srx_loc = url.find("=")
                        srx_id = url[srx_loc+1:len(url)]
                        srx_list.append(srx_id)
                info[field] = '; '.join([x for x in srx_list])
                #info[field] = '; '.join([x.attrib.strip() for x in tmp])  
            else:
                pass      
        elif field == 'Sample/Channel/Characteristics':
            if tmp:
                tx_list = []
                for x in tmp: 
                    if x.attrib['tag'] == 'treatment':
                        tx = x.text.replace("\n", "")
                        tx_list.append(tx)
                info[field] = '; '.join([x for x in tx_list])
                # print(info[field]) 
                #info[field] = '; '.join([x.attrib.strip() for x in tmp])  
            else:
                pass                
        else:
            if tmp:
                info[field] = '; '.join([x.text.strip() for x in tmp])
            else:
                pass   
    return info

def string_found(string1, string2):
    if re.search(r"\b" + re.escape(string1) + r"\b", string2,re.IGNORECASE): 
        return True
    return False


def _matchKeyWord(xmlContent, key, fileds=False):
    ## match key words in a specific xml content
    ## return the match words
    res = {}
    if not fileds:
        fileds = list(xmlContent.keys()) # use all fields if not specified 
    for field in xmlContent.keys(): 
        if field in fileds:
            if string_found(key, xmlContent[field].replace('-', '').replace('_', '').replace('(', '').replace(')', '')):
                tmp = list(set(re.findall(r'%s'%key, xmlContent[field].replace('-', '').replace('_', '').replace('(', '').replace(')', ''),re.IGNORECASE))) 
                if tmp:
                    res[field]= tmp # a list in dict               
    return res

def return_match_res(str,df,xmlContent,fields,key_index=''):
    match_res_luo={}
    li_scale = list(df[str])
    icb_list3=[]
    for term in li_scale:
        if not pd.isnull(term):
            term_li2 = [x.strip().replace('-', '').replace('_', '').replace('(', '').replace(')', '') for x in re.split(';',term)]
            icb_list3 += term_li2   
    icb_list3 = [i for i in icb_list3 if i != '']
    
    for key3 in icb_list3:
        tmp = _matchKeyWord(xmlContent, key = key3, fileds = fields) # result of key word matching, a list in dict
        if tmp:
            # print("found %s key in fields:"%(str), key3)
            for i in tmp.keys():
                match_res_luo[i+key_index].extend(tmp[i]) if i in match_res_luo.keys() else match_res_luo.update({i+key_index:tmp[i]}) 
    return match_res_luo


def _match_crispr(key_file_input,key_logic,xmlContent, fields=False):
    """
    match key words, like crispr, immune, and related cell name, etc
    return a dict, contains geo items and matched keys
    """
    if not xmlContent:
        os.system('echo "No XML content"')
        return {}
    # filter 
    ## crispr apears in sample summary, title

    # match_res={}
    # print(key_file_input)
    # df = np.loadtxt(key_file_input, dtype=str, delimiter=',') #pd.read_table(key_file_input, delimiter = ",") #pd.concat(map(pd.read_csv, key_file_input))# # \t
    # df= pd.DataFrame(df[1:,],columns=df[0,])
    df=pd.read_csv(key_file_input,delimiter=',')

    key_column_names=df.columns
    # print(key_column_names)
    key_match_dict={}
    for i in range(len(key_column_names)):
        # print(i,'#',key_column_names[i])
        key_match_dict.update({'match_res'+str(i):return_match_res(key_column_names[i],df,xmlContent,fields,str(i))})
    
    # match_res=return_match_res('crispr',df,xmlContent,fields)
    # match_res1=return_match_res('scale',df,xmlContent,fields,'1')
    # match_res2=return_match_res('immune',df,xmlContent,fields,'2')
    # match_res3=return_match_res('tumor',df,xmlContent,fields,'3')
    

    if eval(key_logic): #(match_res and match_res1) and (match_res2 or match_res3):
        # match_resf={}
        # match_resf.update(match_res) 
        match_res=key_match_dict
        # for i in range(len(key_column_names)):
            # match_res.update(key_match_dict['match_res'+str(i)])
        # match_res.update(match_res2)
        # match_res.update(match_res3)
    else:
        match_res={}
    if match_res == {}:
        print("no desired keyword was matched")
        return {}
    else:
        print('You got it.')

    # print(match_res)
    return match_res   


def _match_scRNAseq(xmlContent, fields=False):
    """
    match key words, like single cell, single cell RNA-seq, and sequencing platform,etc
    return a dict, contains geo items and matched keys
    """
    if not xmlContent:
        os.system('echo "No XML content"')
        return {}
    #filter bulk RNA-seq in single cell RNA-seq GSE, uploder mentioned single cell but actually bulk RNA-seq 
    ## bulk RNA-seq apears in sample description, title, and source name
    for key in ['bulk rnaseq']:
        bulk1 = _matchKeyWord(xmlContent, key = 'bulk rnaseq', fileds = fields)
        if bulk1:
            os.system('echo "bulk rnaseq"')
            return {}
    ## or library stratey : bulk RNA-seq appears in any fields 
    bulk2 = _matchKeyWord(xmlContent, key = "library strategy: bulk rnaseq")
    if bulk2:
        os.system('echo "bulk rnaseq"')
        return {} # return empty, if bulk RNA-seq appears in Title 
    match_res = {}
    ## 1. match with single cell words, remove special characters, like '-'
    for key1 in ['single cell', 'scrnaseq', 'singlecell rnaseq']:
        tmp = _matchKeyWord(xmlContent, key1)
        if tmp:
            for i in tmp.keys():
                match_res[i].extend(tmp[i]) if i in match_res.keys() else match_res.update({i:tmp[i]})
    # print(match_res)
    tmp = []
    ## 2. match with platform words
    keys = {}
    with open('platform.txt') as f:
        for line in f:
            line = line.rstrip().split('\t')
            keys[line[0]] = line[1] # platforms of sequecing, like 10X Genomics, Smart-seq2
    for key2 in keys.keys():
        tmp = _matchKeyWord(xmlContent, keys[key2]) # result of key word matching, a list in dict
        if tmp:
            for i in tmp.keys():
                # print(tmp)
                match_res[i].extend(tmp[i]) if i in match_res.keys() else match_res.update({i:tmp[i]})
    return match_res


def _match_scATACseq(xmlContent, fields=False):
    """
    match key words, like single cell, single cell RNA-seq, and sequencing platform,etc
    return a dict, contains geo items and matched keys
    """
    if not xmlContent:
        os.system('echo "No XML content"')
        return {}
    #filter bulk RNA-seq in single cell RNA-seq GSE, uploder mentioned single cell but actually bulk RNA-seq 
    ## bulk RNA-seq apears in sample description, title, and source name
    for key in ['bulk atacseq']:
        bulk1 = _matchKeyWord(xmlContent, key = key, fileds = fields)
        if bulk1:
            os.system('echo "%s"' % key)
            return {}
    ## or library stratey : bulk RNA-seq appears in any fields 
    bulk2 = _matchKeyWord(xmlContent, key = "library strategy: bulk atacseq")
    if bulk2:
        os.system('echo "bulk atacseq"')
        return {} # return empty, if bulk RNA-seq appears in Title 
    match_res = {}
    ## 1. match with single cell words, remove special characters, like '-'
    for key1 in ['single cell', 'scatacseq', 'singlecell atacseq', 'singlecell accessiblity']:
        tmp = _matchKeyWord(xmlContent, key1)
        if tmp:
            for i in tmp.keys():
                match_res[i].extend(tmp[i]) if i in match_res.keys() else match_res.update({i:tmp[i]})
    # print(match_res)
    tmp = []
    ## 2. match with platform words
    keys = {}
    with open('platform.txt') as f:
        for line in f:
            line = line.rstrip().split('\t')
            keys[line[0]] = line[1] # platforms of sequecing, like 10X Genomics, Smart-seq2
    for key2 in keys.keys():
        tmp = _matchKeyWord(xmlContent, keys[key2]) # result of key word matching, a list in dict
        if tmp:
            for i in tmp.keys():
                # print(tmp)
                match_res[i].extend(tmp[i]) if i in match_res.keys() else match_res.update({i:tmp[i]})
    return match_res

def _checkSuperSeries(sample_path,xml_path):
    # remove SuperSeries
    xml_string = readGeoXML(sample_path,xml_path)
    if xml_string == None:
        return True
    root = ET.fromstring(xml_string)
    if 'SuperSeries of' in [child.attrib['type'] for child in root.findall('Series/Relation')]:
        return True
    else:
        return False

def _checkType(acc, sample_path, type_need,key_file_input,key_logic,xml_path):
    """
    check whether it is crispr or single cell RNAseq or single cell ATAC-seq
    based on key words matching
    """ 
    ## read in XML to get sample description
    ret = {}
    fields = ['Series/Title', 'Series/Summary', 'Series/Type',
            'Series/Overall-Design'] 
    if sample_path and os.path.isfile(sample_path):
        #NOTE: we need readGeoXML to process
        xmlContent = _getFieldXML(sample_path,xml_path, fields = fields)
        # parse single cell
        if 'crispr' in type_need:
            # single cell RNA-seq, Library-strategy nmush be RNA-Seq
            if ('Expression profiling by high throughput sequencing' in xmlContent['Series/Type']) or ('Genome binding/occupancy profiling by high throughput sequencing' in xmlContent['Series/Type']):
            # if 'genomewide crispr screen' in xmlContent['Series/Type']:
                res = _match_crispr(key_file_input,key_logic, xmlContent, fields)
                if res:
                    ret[acc] = res
                    return ret
            else:
                os.system('echo "%s: neither Expression profiling by high throughput sequencing nor Genome binding/occupancy profiling by high throughput sequencing"'%acc)
    if sample_path and not os.path.isfile(sample_path):
        xml = scrna_parser_from_gse.getGeoXML(acc,xml_path)
        ret = _checkType(acc = acc, sample_path = sample_path, type_need = type_need,key_file_input=key_file_input,key_logic=key_logic,xml_path=xml_path) if xml else None
        return ret

    return None




def getGeoSamples_byType(key_file_input,key_logic,ddir, xml_path,ttype=["crispr"], unique_ids=False, refresh=False):
    """A filter for our Geo model; searcehes our db for the specific sample
    type.
    NOTE: hones in on Library-Strategy tag

    Returns a dictionary of samples fitting the specified

    NOTE: building this up takes time, around 10 secs almost 1 minute!
    TRY: caching the result, reading from a cached file takes only 1 minute
    Store them in files by the .ttype--in the local dir
    """

    ret = {}
    if not unique_ids: #qury all local samples
        #NEED to generate the file, and make the call recursively
        #actually, just one level of recursion b/c geo is pretty flat
        p = os.path.join(ddir)
        ls = os.listdir(p)
        ls = [x for x in ls if x.startswith('GSE') and x != 'GSE'] # for real gsm ids
        for df in ls:
            path = os.path.join(p, df)
            if os.path.isfile(path): #it's a file--check if it's ChIP-Seq
                print("print path of gse.xml", path)
                typo = _checkType(acc=df.split(".")[0], sample_path=path, type_need=ttype,key_file_input=key_file_input,key_logic=key_logic,xml_path=xml_path) # check whether the seq type is chip-seq, atac, or dnase
                if typo:
                    ret.update(typo)
            else:
                newd = os.path.join(ddir, df)
                newdict = getGeoSamples_byType(ddir=newd, ttype=ttype, refresh=refresh,key_file_input=key_file_input,key_logic=key_logic,xml_path=xml_path)
                if newdict and (type(newdict) == type(ret)):
                    ret = dict(ret, **newdict)
    elif unique_ids: # query just for the specified gsm
        for gseid in unique_ids:
            p = os.path.join(ddir+'/'+gseid[:7]+'/'+gseid+'.xml')
            typo = _checkType(acc=gseid, sample_path=p, type_need=ttype,key_file_input=key_file_input,key_logic=key_logic,xml_path=xml_path)
            if typo:
                ret.update(typo)
    else:
        pass
    return ret


def getGeoSamples_byTypes(path, ddir,key_file_input,key_logic,xml_path,datatype=False, gseids=False, refresh=False): #ttypes = ["ATAC-Seq"]): #"ChIP-Seq", "DNase-Hypersensitivity"]):
    ret = []
    if not refresh and os.path.exists(path):
        ret = pickle.load(open(path))
        return ret
    #for t in ttypes:
    if datatype and gseids:
        ret = getGeoSamples_byType(ddir=ddir, ttype=datatype, unique_ids=gseids, refresh=refresh,key_file_input=key_file_input,key_logic=key_logic,xml_path=xml_path)
    elif datatype and not gseids:
        ret = getGeoSamples_byType(ddir=ddir, ttype=datatype, refresh=refresh,key_file_input=key_file_input,key_logic=key_logic,xml_path=xml_path)
    elif gseids and not datatype:
        ret = getGeoSamples_byType(ddir=ddir, unique_ids=gseids, refresh=refresh,key_file_input=key_file_input,key_logic=key_logic,xml_path=xml_path)
    else:
        ret = getGeoSamples_byType(ddir=ddir, refresh=refresh,key_file_input=key_file_input,key_logic=key_logic,xml_path=xml_path)
    # pickle.dump(ret, open(path, "w"))
    return ret
