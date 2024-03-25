# ICRAFT_GEO_parser
crispr screen data parser from GEO
## Clone the repository and do as follow:
1. install a miniconda3, then create a environment by doing:
    ```shell
    conda env create -f crawler_environment.yml
    conda activate py3_geo
    ```
2. test your environment deploy, it is ok, if no error report when you run:
    ```shell
    python env.py
    ```
3. run parser:
    ```shell
    python scrna_parser_runner.py -h
    ```
## Input
| Parameter   | Description |
| ----------- | ----------- |
| -d    | Optional. Parser will get the pubic samples in this given time interval, Please use the format: 2016/01/01-2017/01/01. Default is the recent 100000 entries in GEO. |
| -o    | Required. The table where you wish to store the new sample information, should have a file name ending with '.csv'. By default, we use '|' as the default delimiter to prevent accidental splitting of the content. The example file is "crispr_collection.csv". |
| -aw   | Optional. Upon retrieving the data, the method for writing to the output file can either be append or overwrite. 'a' for append and 'w' for overwrite. The default setting is 'w'.|
| -ai   | Optional. Determines if the output file name should include the time of data retrieval as an index. By default, the index is added. To exclude the index, use -ai ''.|
| -el   | Optional. Path to the file listing the GSE IDs to be excluded. The format should be a single-column file, with each row containing one GSE ID. This is useful in scenarios where certain GSE IDs might cause deadlocks in the scraping process and need to be pre-excluded. |
| -p    | Required. The path where the XML files are stored.|
| -ki   | Required. Path to the file containing keywords. Within each column, keywords should be delineated with semicolons, whereas commas should be used to demarcate keywords between different columns. The keyword matching is case-insensitive.|
| -kl   | Required. Logic between columns of keywords. For example: if the keyword file consists of three columns: "crispr", "cancer", and "immune", the first keyword column should be referred to as "match_res0", and so on. Then this parameter could be written as "(key_match_dict['match_res0']) and (key_match_dict['match_res1'] or key_match_dict['match_res2'])".|

## Usage
    ```shell
    python geo_parser/scrna_parser_runner.py parser -d 2022/12/01-2023/01/01 -el geo_parser/sample/problem_GSE.txt -o ~/sample_collection.csv -p ~/geo_xml -ki '~/sample_key_words.txt' -kl "(key_match_dict['match_res0']) and (key_match_dict['match_res1'] or key_match_dict['match_res2'])"
    ```
    
examples can be found in geo_parser/sample directory, particularly in the sample_crawler.sh file.

## Acknowledgement
This crawler was inspired by and based upon [singleCell_gse_parser](https://github.com/zhengrongbin/singleCell_gse_parser).

## Citation
ICRAFT_A20_paper
