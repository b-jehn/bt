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
### Set up RSHash
Clone the RSHash repository into the code directory:
```
cd ./code
git clone --recursive https://github.com/jonsmcode/rshash
```

Follow the setup steps for RSHash - (also documented in the repositories README):
```

```

### Install bcalm
```
cd ./code
git clone --recursive https://github.com/GATB/bcalm.git

cd bcalm
mkdir -p build
cd build
cmake ..
make -j"$(nproc)"
```
### Set up UST