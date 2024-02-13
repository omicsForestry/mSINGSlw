# mSINGSlw
MSI calling from sequence capture experiments


A lightweight tool to replicate [mSINGS](https://bitbucket.org/uwlabmed/msings/src/master/) without any installation or virtual environment requirements.
It is written to work as it is with minimal dependencies.

## Introduction

mSINGSlw is a tool to call MSI/MSS in tumour sequencing data (specifically from small capture experiments), based on repeat length distributions. The rational being that repeat sites in MSI samples are frequently mutated, resulting a a wider distribution of lengths. If enough reginos are samples, the proportion of mutant sites is higher in MSI samples than MSS.
The original paper describing this is available [here](https://pubmed.ncbi.nlm.nih.gov/24987110/). We have written the basic algorithm into this new tool in order to minimise the requirement of the user to change permissions, create virtual environments or require administrator privileges. We make no claims of novelty, just usefulness.

## Scope

mSINGSlw will provide a score, based on the proportion of mutant repeat regions. For this, you need a bam file for every sample, and a list of regions to test (see below). You will also need a list of known MSS samples as controls, or a precomputed control file, made from MSS samples. No matched normal for the tumour sample is required.

## Installation and requirements

mSINGSlw requires python3, with the numpy and pysam modules installed. We leave it up to the user to install these as they see fit, either with pip, conda or other means.
No installation is required, just download the mSINGSlw.py python script manually, or by:
```
wget site_to_be_added
```

You can even copy/paste the code from this site into your own file. If you want to, you can make the script executable.
```
chmod +x mSINGSlw.py
```
You can move into a directory in your PATH, or make the current directory part of your PATH
```
export PATH=`pwd`:$PATH
```

However, none of this is necessary. The script should run as it is now, but you will have to specify its path when calling, and possibly that it requires python.


## Usage

```
python mSINGSlw.py [-h] [-c CONTROLS] [-p PREBUILT] [-b BUILD] [-s SAMPLES]
                   [-r REGIONS] [-o OUTPUT] [-d DEPTH]

Lightweight rewrite of MSIngs

optional arguments:
  -h, --help            show this help message and exit
  -c CONTROLS, --controls CONTROLS
                        Input list of control bam files
  -p PREBUILT, --prebuilt PREBUILT
                        Prebuilt control file, created by -b option. Either -c
                        or -p must be specified.
  -b BUILD, --build BUILD
                        Filename to save control information from files
                        provided by -c option.
  -s SAMPLES, --samples SAMPLES
                        Sample bam files to test in tab separated: file,
                        sampleID format.
  -r REGIONS, --regions REGIONS
                        Regions to test. Must be specified if -c option is
                        used.
  -o OUTPUT, --output OUTPUT
                        Output file to save if -s option is used to specify
                        samples.
  -d DEPTH, --depth DEPTH
                        Minimum depth required for processing repeat
```

## Options

* **Control samples**
The `CONTROLS` option specifies a file with a list of bam files from known MSS samples eg:
```
pathToControlBam1.bam
pathToControlBam2.bam
pathToControlBam3.bam
...
```
* **Prebuilt controls**
If you want to use the same controls multiple times, you can save the control data with the `BUILD` option. If you have already saved control data with the `BUILD` option, you can reload it with the `PREBUILT` option, and leave out the `CONTROLS` option. The more control samples you provide, the better your calling will be.

* **Test samples**
The `SAMPLES` option specifies a file with a list of bam files to test, and their sample names eg:
```
path2TestBam1.bam	sample1
path2TestBam2.bam	sample2
path2TestBam3.bam	sample3
...
```
If you are merely loading `CONTROLS` samples and saving with the `BUILD` option, then you do not need to specify `SAMPLES`.

* **Test regions**
To either make a prebuilt control using the `BUILD` option or testing `SAMPLES` against `CONTROLS`, you will need to specify the `REGIONS` to test. This option is a filename with information about repeat regions. We have used the format described by [MSIsensor-pro](https://github.com/xjtu-omics/msisensor-pro/) using their scan option:
```
chromosome	location	repeat_unit_length	repeat_unit_binary	repeat_times	left_flank_binary	right_flank_binary	repeat_unit_bases	left_flank_bases	right_flank_bases
chr1	155407852	1	0	12	943	1004	A	TGGTT	TTGTA
chr1	167411405	1	3	11	958	81	T	TGTTG	ACCAC
chr1	182177045	1	3	11	802	23	T	TAGAG	AACCT
```
We provide an example `REGIONS` file in the "examples" folder here, but you will need to make your own depending on your sequencing experiment.

* **Output**
The `OUTPUT` option specifies an output text file to be made. This will provide a list of filenames, and a proportion of mutated sites. Our experience suggests that a score of above 0.2 is generally indicative of MSI, but that will depend on your regions.

* **Depth**
The `DEPTH` option specifies sequence depth below which to not process a region. The default is 30, but you may wish to change this to be more or less stringent. If your data is fairly good quality, lower depth may be fine. If your data is messy, you might want to impose a higher threshold.


## Example usage
We provide example MSS and MSI samples, in bam format, alongside a precomputed control file in the "examples" folder. The test bam files are listed in the `samples.txt` file.
To run mSINGSlw on these samples:
```
python mSINGSlw.py -p preCompControls.txt -s samples.txt -o msi_output.txt
```
The following output file is made:
```
sample	score
sample1	54.166666666666664
sample2	0.0
```
This file indicates that the MSI sample1 has 54% of its repeat regions mutated, whilst the MSS sample2 has none of its regions mutated.
