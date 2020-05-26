# Eukaryotyping: an ensemble similarity-based classification algorithm   

This method constitutes a type of unsupervised [machine learning](https://en.wikipedia.org/wiki/Machine_learning), and was developed to address three issues:

 1. The common occurrence of missing data in genotyping datasets (e.g. such as a situation where 2 or 3 out of 6 multi-locus sequence typing fail to amplify for a subset of your specimens).

 2. The use of MLST methods in the context of sexually reproducing populations, where even closely related individuals may not possess the same genotype (and may be heterozygous) due to chromosomal crossover and random reassortment of chromosomes as occurs during meiosis.

 3. The issue of analyzing specimens which may be extremely complex, potentially representing mixed populations of individuals. Ever try to construct a phylogeny or generate a cluster dendrogram in a situation where for one MLST marker you detect one haplotype, at another you detect three, at another you detect four and another you detect two - in the same specimen? This is essentially what we deal with when we attempt to genotype *Cyclospora cayetanensis* directly     from human stool. It gets extremely complicated.

_Please cite the following [manuscript](https://www.cambridge.org/core/journals/parasitology/article/genotyping-genetically-heterogeneous-cyclospora-cayetanensis-infections-to-complement-epidemiological-case-linkage/0C51FBFFB172DF50357C1D171E9B8657) if you use the Eukaryotyping code in your work:_

```
Barratt, JLN, S Park, FS Nascimento, J Hofstetter, M Plucinski, S Casillas, RS Bradbury, MJ Arrowood, Y Qvarnstrom, E Talundzic (2019) Genotyping genetically heterogeneous Cyclospora cayetanensis infections to complement epidemiological case linkage. Parasitology:1–9 doi:10.1017/S0031182019000581
```

## Getting Started

>These instructions will help you set up a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

This code was developed and tested using a Mac running OSX Catalina 10.15.3. Subsequent instructions are provided _only_ for installing it on an OSX system.

First create a local copy of this repository:

`git clone git@github.com:Joel-Barratt/Eukaryotyping.git` 



### Prerequisites for OSX Catalina

>What prerequisites you need to install this code

#### Xcode Command Line Tools

Install Xcode

```bash
xcode-select --install
```
Check it installed in your $PATH (e.g., /Library/Developer/CommandLineTools)

```bash
xcode-select -p
```

#### Local R package

Go to [CRAN](https://cran.r-project.org/bin/macosx/), download and install `R-3.6.3.pkg (notarized, for Catalina)`

Check that R installed

```bash
R -h   # you should get a help print out with options for R  
```

### Running this code

>While in the folder with all the files from the cloned Eukaryotyping github run:

```bash
Rscript run.r
```
> This will run about 100 samples and produce matrices which can be used for further downstream analysis  


## Additional Information

For additional detailed inforamation about how the Eukaryotyping algorithm works:


Please read [BACKGROUND](background.md) for details on our code of conduct, and the process for submitting pull requests to us.



## Deployment

<!-- need to update once on SciComp and CDCgov github -->

Add additional notes about how to deploy this on a live system


## Contributing
 <!-- need to update @Joel -->


## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc

## Public Domain
This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

## License
The repository utilizes code licensed under the terms of the Apache Software
License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This soruce code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

The source code forked from other open source projects will inherit its license.


## Privacy
This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
Surveillance Platform [Disclaimer](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md)
and [Code of Conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/privacy.html](http://www.cdc.gov/privacy.html).

## Contributing
Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page are subject to the [Presidential Records Act](http://www.archives.gov/about/laws/presidential-records.html)
and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

## Records
This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).

## Notices
Please refer to [CDC's Template Repository](https://github.com/CDCgov/template)
for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/master/CONTRIBUTING.md),
[public domain notices and disclaimers](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md),
and [code of conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
