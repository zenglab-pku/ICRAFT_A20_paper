import argparse
from optparse import OptionParser
import os, sys
# import django
import pickle as p
import env
from scrna_parser_from_gse import _sync_gse
import time
import datetime


def getSyncLog(infoStr):
    """ouput the record to DoneGsmXml.log file
    """
    os.system('echo "[%s] %s"' % (time.strftime('%H:%M:%S'), infoStr))


def convertTime(timeRegion):
    """
    check input time region and split per 31 days
    """
    minString = timeRegion.split('-')[0]
    maxString = timeRegion.split('-')[1]
    t1, t2 = minString.split('/'), maxString.split('/')
    mintime = datetime.datetime(int(t1[0]), int(t1[1]), int(t1[2]))
    maxtime = datetime.datetime(int(t2[0]), int(t2[1]), int(t2[2]))
    deltTime = maxtime - mintime
    getSplitTime = []
    if deltTime > datetime.timedelta(days=31):
        cnt = int(str(deltTime).split(' ')[0]) / 31
        cnt = int(cnt)
        for i in range(1, cnt+2):
            if i == 1:
                start = mintime
                mintime = start + datetime.timedelta(days=31)
                getSplitTime.append('%s-%s'%(start.strftime("%Y/%m/%d"), mintime.strftime("%Y/%m/%d")))
            elif i > 1 and i <= cnt:
                start = mintime + datetime.timedelta(days=1)
                mintime = start + datetime.timedelta(days=31)
                getSplitTime.append('%s-%s'%(start.strftime("%Y/%m/%d"), mintime.strftime("%Y/%m/%d")))
            else:
                start = mintime + datetime.timedelta(days=1)
                getSplitTime.append('%s-%s'%(start.strftime("%Y/%m/%d"), maxtime.strftime("%Y/%m/%d")))
         # output record to log file, so that user known what happened
        getSyncLog("# the input date region include %s"%str(deltTime))
        getSyncLog("# split date region into %d:"%len(getSplitTime))
        for i in getSplitTime:
            getSyncLog(i)
        return getSplitTime
    return [timeRegion]

def main():
    try:
        parser = argparse.ArgumentParser(description="dataset parser from GEO")
        sub_parsers = parser.add_subparsers(help = "sub-command help", dest = "sub_command")
        new_parser = sub_parsers.add_parser("parser",  help = "parse new sample",description = "parse new data from GEO with option of your wanted feature.")
        new_parser.add_argument('-d', dest='date_region', type=str, required = False, help='Parser will get the pubic samples in this given date region, Please use the format: 2016/01/01-2017/01/01. Default is the recent 100000 entries in GEO.')
        new_parser.add_argument('-o', dest='fsave', type=str, required = True,help='csv file with delimiter "|", the table you want to save the new sample information.')
        new_parser.add_argument('-aw', dest='a_or_w', type=str, required = False, default='w',help="Upon retrieving the data, the method for writing to the output file can either be append or overwrite. 'a' for append and 'w' for overwrite.")
        new_parser.add_argument('-ai', dest='add_index_or_not', type=str, required = False, default='add',help="Determines if the output file name should include the time of data retrieval as an index. To exclude the index, use -ai ''")
        new_parser.add_argument('-el', dest='exclude',metavar="FILE", required=False, help="a one-column file contains Accession number that has been parsed")
        new_parser.add_argument('-p', dest="path_of_xml",type = str, required =True, help = 'The folder path contain all the xml files, eg: "./geo", the xml storage format should be: "./geo/GSE1000/GSE1000102.xml".')
        new_parser.add_argument('-ki', dest="key_words_input_file",type = str, required =True, help = 'Path to the file containing keywords. Within each column, keywords should be delineated with semicolons, whereas commas should be used to demarcate keywords between different columns.')
        new_parser.add_argument('-kl', dest="key_column_list_logical_combination",type = str, required =True, help = "Logic between columns of keywords. For example: if the keyword file consists of three columns: 'crispr', 'cancer', and 'immune', the first keyword column should be referred to as 'match_res0', and so on. Then this parameter could be written as '(key_match_dict['match_res0']) and (key_match_dict['match_res1'] or key_match_dict['match_res2'])'.")

        args = parser.parse_args()
        if args.sub_command == 'parser':
            dregion, file_save,a_or_w,add_index_or_not,fill_or_not, excludes, path_of_xml, key_file_input, key_logic = args.date_region, args.fsave,args.a_or_w,args.add_index_or_not, args.fill, args.exclude, args.path_of_xml, args.key_words_input_file,args.key_column_list_logical_combination

            if not dregion:
                getSyncLog("No date region setting!")
                dregion = False
            else:
                checkTime = convertTime(dregion)

            for oneTime in checkTime:
                dregion = oneTime
                getSyncLog("#+++++++++ New Collection in %s"%dregion)
                getSyncLog('parse new sample and add in database')
                _sync_gse(file_save,a_or_w,add_index_or_not, key_file_input,key_logic, path_of_xml,fill_or_not, dateRegion = dregion, exludeFile=excludes, refresh=True)                

            
    except KeyboardInterrupt:
        sys.stderr.write("User interrupted me!\n")
        sys.exit(0)

if __name__ == '__main__':
    main()
    #signal.signal(signal.SIGALRM, handler)
