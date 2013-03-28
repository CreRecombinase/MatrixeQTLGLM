#!/bin/bash
#Usage import-index.sh [-i (import) [-s (snps)| -g (genes)]|-x (index) -f <file> -d <database>

import=false
index=false
snps=false
genes=false
eqtls=false
database=



while getopts "ixsgef:d:" opt; do
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
	e)
	    eqtls=true
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
	echo -e ".timeout 20000\npragma main.page_size=4096;pragma main.cache_size=10000; pragma synchronous=0;pragma main.journal_mode=WAL; pragma temp_store=1;pragma temp_store_directory='.';create index ss on snps(Snp,Sample);" > index_metafile_snps.sql
	sqlite3 $database < index_metafile_snps.sql
	rm index_metafile_snps.sql
	else
	if $genes; then
	    echo -e ".timeout 20000 \npragma main.page_size=4096;pragma main.cache_size=10000; pragma synchronous=0;pragma main.journal_mode=WAL;pragma main.cache_size=5000; pragma temp_store=1;pragma temp_store_directory='.';create index gs on gene(Gene,Sample);" > index_metafile_gene.sql
	    sqlite3 $database < index_metafile_gene.sql
	    rm index_metafile_gene.sql
	    else
	    if $eqtls; then
		echo -e ".timeout 20000\npragma main.page_size=4096;pragma main.cache_size=10000;pragma synchronous=0; pragma main.journal_mode=WAL;pragma main.cache_size=5000;pragma temp_stire=1;pragma temp_store_directory='.';create index sgkc on eqtls(Gene,SNP,Kfold);" > index_eqtlfile.sql
		sqlite3 $database < index_eqtlfile.sql
		rm index_eqtlfile.sql
	    fi
	fi
    fi
elif $import; then
    if $snps; then
	echo -e ".timeout 20000\npragma journal_mode=memory; pragma synchronous=0; pragma cache_size=250000;create table snps(Snp TEXT, Sample TEXT, Value INTEGER);\n.separator \"\\t\"\n.headers on\n.import $file snps" > snp_metafile.sql
	sqlite3  $database <snp_metafile.sql
	rm snp_metafile.sql
	else
	echo -e ".timeout 20000\npragma journal_mode=memory; pragma synchronous=0;\n pragma cache_size=250000;\n create table gene(Gene TEXT, Sample TEXT, Value REAL);\n.separator \"\\t\"\n.headers on\n.import $file gene" >gene_metafile.sql
	sqlite3 $database <gene_metafile.sql
	rm gene_metafile.sql
	fi
fi
