## Introduction

ALLHiC_adjuster is a toolset contains several scripts for adjusting ALLHiC result.

## Dependencies

### Software

[jcvi](<https://github.com/tanghaibao/jcvi>)

[ALLHiC](<https://github.com/tangerzhang/ALLHiC>)

### Python Modules

matplotlib

## Installation

```bash
cd /path/to/install
git clone https://github.com/sc-zhang/ALLHiC_adjuster.git
pip install -r requirements.txt
chmod +x ALLHiC_adjuster/ALLHiC_adjuster.py
echo 'export PATH=/path/to/install/ALLHiC_adjuster:$PATH' >> ~/.bash_profile
source ~/.bash_profile
```

## Usage

### Main program

```bash
usage: ALLHiC_adjuster.py [-h] {locator,extractor,convertor,adjuster,builder} ...

options:
  -h, --help            show this help message and exit

sub commands:
  {locator,extractor,convertor,adjuster,builder}
    locator             Visualize each block of chromosome
    extractor           Extract sequences
    convertor           Convert files
    adjuster            Adjust tour or group files
    builder             Build chromosome-level fasta from tour files
```

### Step 1. Locate each contig block with jcvi anchors

```bash
usage: ALLHiC_adjuster.py locator [-h] -q QUERY -r REFERENCE -c ANCHORS -a AGP [-s RESOLUTION] -o OUTPIC

options:
  -h, --help            show this help message and exit
  -q QUERY, --query QUERY
                        Query bed file
  -r REFERENCE, --reference REFERENCE
                        Reference bed file
  -c ANCHORS, --anchors ANCHORS
                        Query.Reference.anchors file, generated by jcvi
  -a AGP, --agp AGP     AGP file of query genome
  -s RESOLUTION, --resolution RESOLUTION
                        Resolution means 1/resolution of chromosome length, default=20
  -o OUTPIC, --outpic OUTPIC
                        Output picture
```

Example:

```bash
ALLHiC_adjuster.py locator -q query.bed -r reference.bed -c query.reference.anchors -a query.agp -s 100
```

The picture will figure out each block, and show start contig and end contig of each block, and a .block.txt file with
same name with picture will save contig list in each block.

### Step 2. Adjust

```bash
usage: ALLHiC_adjuster.py adjuster [-h] {merge,reverse,split} ...

options:
  -h, --help            show this help message and exit

Adjuster commands:
  {merge,reverse,split}
    merge               Merge tour files
    reverse             Reverse tour file
    split               Split tour file or txt file
```

This three sub commands of adjuster sub command can be use with tour files, for detail, you can use commands below for
more help

```bash
ALLHiC_adjuster.py adjuster merge -h
ALLHiC_adjuster.py adjuster reverse -h
ALLHiC_adjuster.py adjuster split -h
```

Example:

```bash
# merge tours
ALLHiC_adjuster.py adjuster merge -i <groupX.tour,groupY.tour> -o <groupM.tour>

# reverse tours
ALLHiC_adjuster.py adjuster reverse -i <groupX.tour>

# split tour
ALLHiC_adjuster.py adjuster split -i <groupX.tour> -c <contig1,contig2>
```

### Step 3. Build

Make sure the folder for input only contain necessary .tour files then run:

```bash
usage: ALLHiC_adjuster.py builder [-h] -r REF -i INPUT [-o OUTPUT]

options:
  -h, --help            show this help message and exit
  -r REF, --ref REF     Input contig-level fasta file
  -i INPUT, --input INPUT
                        Input directory contain all required tour files
  -o OUTPUT, --output OUTPUT
                        Output directory of groups.asm.fasta file and groups.agp file, default='.'
```

## Otherwise

Extract sequences

```bash
usage: ALLHiC_adjuster.py extractor [-h] {tour,list} ...

options:
  -h, --help   show this help message and exit

extractor commands:
  {tour,list}
    tour       Extract sequences with tour
    list       Extract sequences with list
```

Convert files

```bash
usage: ALLHiC_adjuster.py convertor [-h] {agp2tour,tour2txt,tours2cluster,txt2cluster,anchors2circos,ragoo2agp,ragoo2tour} ...

options:
  -h, --help            show this help message and exit

convertor commands:
  {agp2tour,tour2txt,tours2cluster,txt2cluster,anchors2circos,ragoo2agp,ragoo2tour}
    agp2tour            Convert AGP file to tour files
    tour2txt            Convert tour to txt file
    tours2cluster       Convert tour files to cluster file
    txt2cluster         Convert txt files to cluster file
    anchors2circos      Convert anchors to circos link file
    ragoo2agp           Convert RaGOO ordering files to AGP file
    ragoo2tour          Convert single RaGOO ordering file to tour file
```

For more information, use command below:

```bash
ALLHiC_adjuster.py subcommand1 -h
ALLHiC_adjuster.py subcommand1 subcommand2 -h
```
