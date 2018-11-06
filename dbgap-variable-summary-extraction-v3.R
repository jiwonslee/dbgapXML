library(xml2)
library(RCurl)
library(stringr)
library(argparser)
require(dplyr)
#setwd('C:/Users/jq650/Desktop/work/phenotype_harmonization/dbgap-dict-pulls')

argp <- arg_parser("Parse dbgap xml var reports and data dicts")
argp <- add_argument(argp, "study", help="enter phs accession number of study")
argp <- add_argument(argp, "--table", help="enter pht accession number of table)", type="character")
argp <- add_argument(argp, '--dir', help='enter directory where files should be output')
argv <- parse_args(argp)

study = argv$study
table = argv$table
message(paste0('searching for study ', study))
if (!is.na(table)){
  message(paste0('searching for table ', table))
}else{
  message('warning: no --table set. Information of all tables are being extracted')
}
#args = commandArgs(trailingOnly=TRUE)
#study = args[1]
#table = args[2]
#study = 'phs000007'
#table = 'pht000024'
dir = argv$dir
if (!is.na(dir)){
  setwd(dir)
}

base_url = 'ftp://ftp.ncbi.nlm.nih.gov/dbgap/studies/'
study_url = paste0(base_url, study)
print(paste("study url", study_url))

###### needs suffix / for the next line to work...
folders <- unlist(strsplit(RCurl::getURL(paste0(study_url, '/'), ftp.use.epsv = TRUE, dirlistonly = TRUE), "\n"))

max_v = -1
max_p = -1
for (i in 1:length(folders)) {
    current = unlist(strsplit(folders[i], "\\."))
    
    if (current[1] == study) {
        current_v = as.numeric(substring(current[2], 2))
        current_p = as.numeric(substring(current[3], 2))
        
        if (current_v > max_v) {
           max_v = current_v
           max_p = current_p
        } else if (current_v == max_v && current_p > max_p) {
           max_p = current_p
        }
    }
    current = NULL
}
latest_study = paste0(study, ".v", toString(max_v), ".p", toString(max_p))

main_url = paste(study_url, latest_study, 'pheno_variable_summaries/', sep="/")
#print(paste("latest study url", main_url))

xml_path_suffixes <- unlist(strsplit(RCurl::getURL(main_url, ftp.use.epsv = TRUE, dirlistonly = TRUE), "\n"))
xml_path_suffixes <- trimws(xml_path_suffixes)

if (!is.na(table)){
  xml_path_suffixes <- xml_path_suffixes[as.vector(sapply(xml_path_suffixes, str_detect, table))]
}

message(paste0('extracting information from ', length(xml_path_suffixes), ' files'))

#for (i in 1:length(xml_path_suffixes)) {
for (i in 1:length(xml_path_suffixes)) {
    xml_url = paste0(main_url, xml_path_suffixes[i])
    
    # data_dict and var_report do and don't contain 'pNN' patch levels of table versions:
    # phs000007.v29.pht004361.v2.l_bamyl_ex02_3b_0818s.data_dict.xml
    # phs000007.v29.pht004361.v2.p10.l_bamyl_ex02_3b_0818s.var_report.xml
    
    current = unlist(strsplit(xml_path_suffixes[i], "\\."))
    # don't parse meta data xml
    if (length(current) > 5) {
      study=current[1]
      study_version=current[2]
        if (current[6] == 'data_dict') {
            table = current[3]
            table_version = current[4]
            table_patch = 'NA'
            table_name = current[5]
            doc='data_dict'
        } else if (current[7] == 'var_report') {
            table = current[3]
            table_version = current[4]
            table_patch = current[5]
            table_name = current[6]
            doc='var_report'
        } else {
            table = paste('FIXME', xml_path_suffixes[i])
            table_version = paste('FIXME', xml_path_suffixes[i])
            table_patch = paste('FIXME', xml_path_suffixes[i])
            table_name = paste('FIXME', xml_path_suffixes[i])
            print(table)
            message('error in table name format - check and modify code')
            quit()
        }

        xml = read_xml(RCurl::getURLContent(xml_url))

        recs <- xml_find_all(xml, "//variable")
        vals = trimws(xml_text(recs))
        table_parts <- as.vector(str_split(xml_path_suffixes[i], '\\.', simplify=TRUE))
        # extract and clean (if needed) the area names
        variable_ids <- trimws(xml_attr(recs, "id"))
        descriptions <- trimws(xml_text(xml_find_all(xml, "//variable//description")))
        if (doc=='var_report'){
          
          variable_name <- trimws(xml_attr(recs, "var_name"))
          type <- trimws(xml_attr(recs, "calculated_type"))
          units <- trimws(xml_attr(recs, 'units'))
          stats <- xml_find_all(recs, '//total//stats')
          stat <- xml_find_first(stats, 'stat')
          n <- trimws(xml_attr(stat, 'n'))
          nulls <- trimws(xml_attr(stat, 'nulls'))
          mean <- trimws(xml_attr(stat, 'mean'))
          sd <- trimws(xml_attr(stat, 'sd'))
          median <- trimws(xml_attr(stat, 'median'))
          min <- trimws(xml_attr(stat, 'min'))
          max <- trimws(xml_attr(stat, 'max'))
          df <- cbind(variable_ids, variable_name, type , units, descriptions, n, nulls, mean, sd, median, min, max, study, study_version, table, table_version, table_name, table_patch)
        }
        else if (doc=='data_dict'){
          phv_parts <- unlist(sapply(variable_ids, str_split, pattern='\\.'))
          phv <- phv_parts[seq(1, length(phv_parts)), by=2]
          phv_version <- phv_parts[seq(2, length(phv_parts)), by=2]
          variable_name <- trimws(xml_text(xml_find_first(recs, 'name')))
          units <- trimws(xml_text(xml_find_first(recs, 'unit')))
          type <- trimws(xml_text(xml_find_first(recs, 'type')))
          df <- cbind(phv, phv_version, variable_name, units, descriptions, study, study_version, table, table_version, table_name)
        }
        fname <- str_replace(xml_path_suffixes[i], ".xml", ".txt")
        write.table(df, sep='\t', file=fname, row.names=FALSE, quote=FALSE)
        
    }
}
