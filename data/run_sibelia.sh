#!/bin/bash

SIB_DIR=$HOME/Bioinf/Sibelia/distr
GENOMES=$HOME/Bioinf/refass/data/references

DIR=`pwd`
cd $SIB_DIR
bin/Sibelia -s loose -o $DIR/$2 $GENOMES/DH1.fasta $GENOMES/UMNK88.fasta $GENOMES/Xuzhou21.fasta $GENOMES/W.fasta $GENOMES/MG1655-K12.fasta $DIR/$1
cd $DIR/$2
rm -r circos
rm coverage_report.txt
rm d3_blocks_diagram.html