#!/bin/bash

## to improve
 
 meryl count k=21 hifi_reads.fq output meryl_db
cd /opt/merqury-1.3
export MERQURY=$PWD 
./merqury.sh read-db.meryl ../Tethysbaena_scabra.fa merq



 
