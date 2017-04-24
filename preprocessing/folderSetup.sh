#!/bin/bash

mkdir cosmic_exac
mkdir cosmic_exac/coinciding_somatic
mkdir cosmic_exac/coinciding_germline
mkdir cosmic_exac/unique_somatic
mkdir cosmic_exac/unique_germline
cp -r cosmic_exac cosmic_oneKG
cp -r cosmic_exac cosmic_dbSNP
cp -r cosmic_exac icgc_exac
cp -r cosmic_exac icgc_oneKG
cp -r cosmic_exac icgc_dbSNP
mkdir chromFA

echo "new folders created"
