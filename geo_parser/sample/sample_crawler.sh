# conda activate py3_geo
cd /home/apps/geo_parser/sample
python ../scrna_parser_runner.py parser -d 2022/12/01-2023/01/01 -el ./problem_GSE.txt -o sample_collection.csv -p geo_xml -ki 'sample_key_words.txt' -kl "(key_match_dict['match_res0']) and (key_match_dict['match_res1'] or key_match_dict['match_res2'])"

# An example of a GSE ID that requires manual verification for compliance and exclusion from the crawler.
# GSE126374
# still did not pass
# link problem: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126374&view=full&form=xml&targ=all
# GSE126374
# ERROR: accession is NOT xml. (The accession may be deleted from GEO repository)
# opening url
# [09:34:10] local IP: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126374&view=full&form=xml&targ=all
# IsADirectoryError: [Errno 21] Is a directory: 'geo_xml/GSE2122/GSE212246.xml'