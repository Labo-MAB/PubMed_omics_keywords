#! /usr/bin/gawk -f
# Author Savvas Paragkamian
#
# Run :  gunzip -c path-to-pubmed/*.tsv.gz | ./scripts/search_engine.awk keywords.txt - > results.tsv
# Note the - , it takes the input from the pipe |
#

BEGIN {
    FS="\t"
}

# Treat keywords

(ARGIND==1){
    keywords[$1]=tolower($1)
}

# Treat corpus. The input comes from the - from the execution of the script
# old col order: PMID; Authors; Journal.volume:pages; year; title; abstract
# new col order: PMID; Authors; Journal.volume:pages; year; title; keywords; abstract
(ARGIND==2){

    if ( $4 ~ /[[:alnum:]]/ && $5 ~ /[[:alnum:]]/ && $6 ~ /[[:alnum:]]/) {

        line=tolower($5$6$7) # combine title, keywords and abstract and make all text lowercase

        sub("\\|.*$","",$1) # remove the doi addresses from all the PMIDs
        
        for (k in keywords){
            
            if (match(line,"[^[:alpha:]]" keywords[k] "[^[:alpha:]]")){
                print $4 "\t" $1 "\t" keywords[k] "\t" $6 # print year  PMID	keywords(of article)    keyword(of given keywords)
            }
        }
    }
}
