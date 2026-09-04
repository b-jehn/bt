# Bachelors Thesis - Bjarne Jehn

## Repository Structure:
```
.
├── TexStudio               #contains everything tex related
│   └── Thesis              #main document
│       ├── MINLP           #chapter: MINLP     
│       ├── figures 
│       └── subsampling     #chapter: subsampling
└── code                    #contains everything code related
    └── lp_work             #should be MINLP work
        └── minimizer_stats #minimizer distributions
            ├── cod 
            └── homo_sapiens
```

## Preperation
### Set up RSHash in the ./code Directory
```
cd code
git clone --recursive https://github.com/jonsmcode/rshash

mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release .. -D CMAKE_CXX_COMPILER=g++-14
make
```
"RSHash is a compressed, exact data structure for k-mers (strings of length k over the DNA alphabet {A,C,G,T}) based on bitvectors with Rank and Select support to Hash k-mers." - (https://github.com/jonsmcode/rshash/blob/main/README.md)

### Install bcalm in the ./code Directory
```
cd code
git clone --recursive https://github.com/GATB/bcalm.git

cd bcalm
mkdir -p build
cd build
cmake ..
make -j"$(nproc)"
```
"BCALM 2 is a bioinformatics tool for constructing the compacted de Bruijn graph from sequencing data." - (https://github.com/GATB/bcalm/blob/master/README.md)

### Set up UST in the ./code Directory
```
cd code
git clone https://github.com/medvedevgroup/UST.git

cd UST
make
```
"UST is a bioinformatics tool for constructing a spectrum-preserving string set (SPSS) representation from sets of k-mers." - (https://github.com/medvedevgroup/UST/blob/master/README.md)

### Benchmarks
Get datasets used by Pibiri (2022) with
```
cd code/rshash
./benchmarks/download_and_preprocess_datasets.sh
```

## Run RSHash
```
./build/source/rshash build -i ./DNA_datasets/Gadus_morhua.gadMor3.0.dna.toplevel.fa -d ./DNA_datasets/Gadus_morhua.rshash -k 31 --m1 18 --t1 64 --level 1
./build/source/rshash query -d ./DNA_datasets/Gadus_morhua.rshash -q ./DNA_datasets/Gadus_morhua.fastq.gz

./build/source/rshash build -i ./DNA_datasets/Gadus_morhua.gadMor3.0.dna.toplevel.fa -d ./DNA_datasets/Gadus_morhua.rshash -k 31 --m1 18 --t1 64 --m2 18 --t2 64 --level 2
```