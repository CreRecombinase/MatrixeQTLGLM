#!/bin/bash
#Usage import-index.sh [-i (import) [-s (snps)| -g (genes)]|-x (index) -f <file> -d <database>
import=false
index=false
snps=false
genes=false
database=



while getopts "ixsgf:d:" opt; do
    case $opt in 
	i)

	    import=true;
	    ;;
	x)

	    index=true;
	    ;;
	s)
	    snps=true

	    ;;
	g)
	    genes=true

	    ;;
	f)
	    file=$OPTARG
	    ;;
	d) 

	    database=$OPTARG
	    ;;
	
	\?)
	    echo "Invalid option -$OPTARG" >&2
	    ;;
    esac
done
if $index; then
    if $snps; then
	sqlite3 $database "pragma main.page_size=4096;pragma main.cache_size=10000; pragma main.locking_mode=EXCLUSIVE;pragma main.journal_mode=WAL;pragma main.cache_size=5000; pragma temp_store=1;pragma temp_store_directory='.';create index ss on snps(Snp,Sample)"
	else
	sqlite3 $database "pragma main.page_size=4096;pragma main.cache_size=10000; pragma main.locking_mode=EXCLUSIVE;pragma main.journal_mode=WAL;pragma main.cache_size=5000; pragma temp_store=1;pragma temp_store_directory='.';create index gs on gene(Gene,Sample)"
    fi
    #sqlite3 $database "pragma journal_mode=memory; pragma synchronous=0; pragma cache_size=500000; create 
elif $import; then
    if $snps; then
	echo -e "pragma journal_mode=memory; pragma synchronous=0; pragma cache_size=250000;create table snps(Snp TEXT, Sample TEXT, Value INTEGER);\n.separator \"\\t\"\n.headers on\n.import $file snps" > snp_metafile.sql
	sqlite3  $database <snp_metafile.sql
	rm snp_metafile.sql
	else
	echo -e "pragma journal_mode=memory; pragma synchronous=0;\n pragma cache_size=250000;\n create table gene(Gene TEXT, Sample TEXT, Value REAL);\n.separator \"\\t\"\n.headers on\n.import $file gene" >gene_metafile.sql
	sqlite3 $database <gene_metafile.sql
	rm gene_metafile.sql
	fi
fi
