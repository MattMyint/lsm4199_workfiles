#!/bin/bash


## unlike SRA, EGA provides urls for simple download
## random wait ensures that server doesn't detect automation and reject the request
echo Starting Embryo Fastq
wget --random-wait -i embryo_sample_urls.txt -P embryo/
