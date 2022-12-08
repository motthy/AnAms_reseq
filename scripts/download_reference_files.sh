#!/bin/bash

mkdir reference
cd reference
echo Downloading reference files from Cats-I...
curl -LO https://cat.annotation.jp/download/AnAms1.0/AnAms1.0.genome.fa.gz
# curl -LO https://cat.annotation.jp/download/AnAms1.0/AnAms1.0r1.0.2.gff.gz
curl -LO https://cat.annotation.jp/download/AnAms1.0/AnAms1.0r1.0.2.prot.fa.gz
echo Done.
echo Uncompressing files...
gunzip *.gz
echo Done.