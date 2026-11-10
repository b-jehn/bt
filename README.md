# Bachelors Thesis - Bjarne Jehn
Instructios to set up the repository and run the computation of stats needed for running and evaluating the MINLP

## Clone My "bt" Repository:
```
git clone https://github.com/b-jehn/bt.git
```
Structure:
```
.
├── TexStudio               #contains everything tex related
│   └── Thesis              #main document
└── code                    #contains everything code related
    └── lp_work             #should be MINLP work
        └── minimizer_stats     #minimizer distributions for one level - outdated
        └── leveledDistribution #leveled minimizer distributions
```

## Prerequisites for Experiments

### Clone the "RSHash" Repository Into ./code 
```
cd bt/code
git clone --recursive https://github.com/jonsmcode/rshash
```
"RSHash is a compressed, exact data structure for k-mers (strings of length k over the DNA alphabet {A,C,G,T}) based on bitvectors with Rank and Select support to Hash k-mers." - (https://github.com/jonsmcode/rshash/blob/main/README.md)

### Replace rshash.hpp, main.cpp and build.cpp
Replace the three files in the rshash directory
```
bt/code/rshash/source/rshash.hpp
bt/code/rshash/source/main.cpp
bt/code/rshash/source/build.cpp
```
with those needed for the bruteforce stats computation 
```
bt/code/lp_work/bruteForceStats/rshash_changes/rshash.hpp
bt/code/lp_work/bruteForceStats/rshash_changes/main.cpp
bt/code/lp_work/bruteForceStats/rshash_changes/build.cpp
```

### Build RSHash
```
cd bt/code/rshash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release .. -D CMAKE_CXX_COMPILER=g++-14
make
```


### (Install bcalm in ./code Directory) - only if necessary to load the datasets
```
cd bt/code
git clone --recursive https://github.com/GATB/bcalm.git

cd ./bcalm
mkdir -p build
cd ./build
cmake ..
make -j"$(nproc)"
```
"BCALM 2 is a bioinformatics tool for constructing the compacted de Bruijn graph from sequencing data." - (https://github.com/GATB/bcalm/blob/master/README.md)

### (Set up UST in ./code Directory) - only if necessary to load the datasets
```
cd bt/code
git clone https://github.com/medvedevgroup/UST.git

cd ./UST
make
```
"UST is a bioinformatics tool for constructing a spectrum-preserving string set (SPSS) representation from sets of k-mers." - (https://github.com/medvedevgroup/UST/blob/master/README.md)

### Download and Preprocess Datasets
Get datasets used by Pibiri (2022) with
```
cd bt/code/rshash
./benchmarks/download_and_preprocess_datasets.sh
```
**In order to run the experiments, rshash needs to have been build in the regarding directory and the original and preprocessed datasets for cod, falcon and human genome have to be stored in** 
```
bt/code/rshash/DNA_datasets/Gadus_morhua.gadMor3.0.dna.toplevel.fa
bt/code/rshash/DNA_datasets/Gadus_morhua.gadMor3.0.dna.toplevel.fa.unitigs.fa

bt/code/rshash/DNA_datasets/Falco_tinnunculus.FalTin1.0.dna.toplevel.fa
bt/code/rshash/DNA_datasets/Falco_tinnunculus.FalTin1.0.dna.toplevel.fa.unitigs.fa

bt/code/rshash/DNA_datasets/Homo_sapiens.GRCh38.dna.toplevel.fa
bt/code/rshash/DNA_datasets/Homo_sapiens.GRCh38.dna.toplevel.fa.unitigs.fa
```
this should be the case after following the steps for downloading and preprocessing the datasets

## Run Experiments
The only parameter that needs to be specified for the computation of stats is the number of threads to be used.
```
cd bt/code/lp_work/bruteForceStats
python3 bruteForceStats.py <number of threads>
```
Once the script has run, the results should be contain in 
```
bt/code/lp_work/bruteForceStats/allStats_NONStreamingLookup.zip
```
That file is all I need :) 
