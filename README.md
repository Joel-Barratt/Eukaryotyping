# Eukaryotyping
This repository contains the algorithms first described by Barratt, Park and colleagues in the manuscript that can be accessed under the following DOI: https://doi.org/10.1017/S0031182019000581
"Eukaryotyping" refers to the analysis performed using the ensemble of similarity-based classification algorithms developed by Dr Mateusz Plucinski and Dr Joel Barratt. This analysis approach was initially designed for analysis of complex Cyclospora cayetanensis genotyping datasets but we are applying these algorithms to other genotyping challenges. This approach constitutes a type of unsupervised machine learning, and was developed to tackle three major issues: (1) the common occurrence of missing data in genotyping datasets (i.e. lets say, 3 out of 6 MLST markers fail to amplify for a subset of your specimens), (2) the use of MLST methods in the context of sexually reproducing populations, where even closely related individuals may not possess the same genotype (and may be heterozygous) due to chromosomal crossover and random reassortment of chromosomes as occurs during meiosis and, (3) the issue of analyzing specimens which may be extremely complex, potentially representing mixed populations of individuals. Ever try to perform a phylogeny or generate a cluster dendrogram in a situation where for one MLST marker you detect one haplotype, at another you detect three, at another you detect four and another you detect two? This is essentially what we deal with when we attempt to genotype Cyclospora cayetanensis directly from human stool. It gets extremely complicated.

Using the "Eukaryotyping" method you do not need to exclude specimens with a partial genotype (within reason). For example, if you have a large MLST dataset where some specimens have only half of their MLST markers successfully sequenced, this method will usually cope just fine provided the markers that were successfully amplified are sufficiently diverse. If you have a complex dataset where some specimens are heterozygous or homozygous at any marker, this method will also cope just fine. Even if your dataset represents a population of mixed genotypes (i.e., such as a stool specimen containing a set of haplotypes representing a sexually reproducing population of Cyclospora), this method should cope just fine. The end result is a distance matrix where the values (all between 0 and 1) represent how closely related every possible pair of specimens is. A value of "0" indicates that the specimens are identical, and a value of "1" indicates the specimens are not related.

The R scripts provided here constitute the earliest version of the "Eukaryotyping" code that was initially prepared for analysis of our complex Cyclospora cayetanensis genotyping data. These files were written by Dr Mateusz Plucinski who is a co-author of the manuscript that first describes the algorithms: https://doi.org/10.1017/S0031182019000581

Five files are provided here (six if you include the README file). One of these files (import_data_V2.r) must be modified by the user depending on how the user wishes the haplotypes to be defined, and depending on the name of the file/data sheet containing the haplotype data. An example of how the haplotype data sheet must be formatted is provided in this file: "Example_haplotype_data_sheet.txt". This format (developed by Dr Joel Barratt) must be adhered to or the scripts will not run correctly. The names of the markers used can be changed by the end user, but the format must not change. The names of the markers (i.e., the "locinames" variable) must be modified in the "import_data_V2.r" script(line 6) according to a system defined by the user. For example, If the user has several large genes that they wish to split into smaller parts, there is the option to do so by defining them as GENE_1_PART_A, GENE_1_PART_B, GENE_1_PART_C etc.., and the base name (i.e., the locinames_base variable - line 9) of this marker is then defined as "GENE_1" (take a look at "import_data_V2.r" and you will get the idea). If the end user does not wish to divide a given marker (e.g. "MARKER_Z") into sub-markers, then the base name is simply "MARKER_Z" - the same as the "locinames" variable. Ultimately, it is the end user that develops a system for defining haplotypes that best suits their purpose.

There are some caveats that must be considered when dividing up markers in the way described above, and this would depend on the entropy of the marker and the number of markers being examined. Ultimately, if the marker is extermely hypervariable (e.g., lets say, there are 50 haplotypes in your population of 60 specimens), very few links will be identified at this locus by the algorithm. In this case, the end user may wish to include additional markers possessing fewer haplotypes (i.e., a lower entropy), or split this hypervariable locus into smaller pieces such that the number of haplotypes at each piece is drastically reduced. Alternatively, the end user may wish to exclude such hypervariable markers as loci possessing extreme variation are not ideal for genotyping in general.

The user also needs to define the ploidy of each marker in the "import_data_V2.r" script - line 11. For example, if your markers are mitochondrial, then you would indicate "1" for the ploidy as the mitochondrial genome is inherited via an extranuclear mechanism. If your organism is a sexually reproducing eukaryote with 2N chromosomes, nuclear/chromosomal loci must have their ploidy set to "2". The ploidy must be defined for all markers previously included in the "locinames" variable (and in the same order) in the import sheet (not the base names - if you do this for only the base names, the script will fail to run). The scripts provided here require *a combination* of nuclear and mitochondrial markers in order to run correctly.

Next, the user might consider modifying the "euk_bayesian_fulldataset_V2.r" script to provide an epsilon value. This value reflects the rate of missing data or how often markers or haplotypes that *should* be present in a specimen are *not* detected due to amplification and/or sequencing failures. We determined this value to be 0.3072 based on a series of deep amplicon sequencing (Illumina) experiments, performed on a large Cyclospora dataset, considering 8 genotyping markers. The user is open to run the algorithm using this epsilon value (this will probably provide a sensible answer), or the user may wish to run experiments to determine this value for their own dataset. For the sake of accuracy, we recommend that users determine an epsilon value that is specific for their dataset.

The user should also modify the filtering parameters for their genotyping dataset. This information is contained in the "cleandata" variable within the "import_data_V2.r" script - line 54. This variable will define the inclusion criteria for your specimens. For instance, if 1 out of 8 typing markers was successfully sequenced for a given specimen, you may wish to exclude that specimen as this is not enough information to make an accurate cluster assignment. You may wish to only include specimens that sequenced successfully for at least markers 1, 2 and 8, plus at least one other marker, or only include specimens that sequenced for at least markers 1 to 4. These combinations are up to the user to decide upon. When modified correctly, the script will automatically exclude specimens that fail to meet the criteria defined. Note that these loci are defined numerically, and these numbers correspond to the precise order the markers are listed in, within the "locinames_base" variable of the "import_data_V2.r" (line 9).

Once the "import_data_V2.r" and  "euk_bayesian_fulldataset_V2.r" are modified by the user to reflect their data, calculation of two matrices by each algorithm is commenced by running the script "run.r". This script first initiates the importation of data using the "import_data_V2.r" script. Next, calculation of the Bayesian matrix is initiated by "run.r" which runs the script "euk_bayesian_fulldataset_V2.r". This is the Bayesian algorithm developed by Dr Mateusz Plucinski. Once the Bayesian matrix is calculated, the "run.r" script then initiates calculation of the heuristic matrix by running the script:  "euk_heuristic_fulldataset.r". This script is the heuristic algorithm developed by Dr Joel Barratt.  The script called "euk_heuristic_fulldataset.r" does not need to be modified by the end user at any point.

These scripts will generate multiple files including the Bayesian matrix, the heuristic matrix and the ensemble matrix. In the manuscript available under DOI https://doi.org/10.1017/S0031182019000581 the output values of the Bayesian matrix and the heuristic matrix were simply averaged to generate the ensemble matrix. However, we later found that the predictive value of this method was improved using an alternative normalization approach. This improved approach involves mapping the distribution of distances generated by the Bayesian algorithm to the empiric distribution of distances generated by the heuristic algorithm. The mean of these mapped pairs is then used to generate the final ensemble matrix of pairwise distances. The scripts provided here perform the latter (i.e., improved) normalization procedure.

Finally, the ensemble matrix is clustered and visualized as a dendrogram. We prefer to use Wards clustering method and Manhattan distances, as described here: https://doi.org/10.1017/S0031182019000581

The manuscript describing this method and the scripts provided here have been distributed under the terms of the Creative Commons Attribution-NonCommercial-ShareAlike licence (http://creativecommons.org/licenses/by-nc-sa/4.0/), which permits non-commercial re-use, distribution, and reproduction in any medium, provided the same Creative Commons licence is included and the original work is properly cited.

If you make use of these files and find our work helpful, please don't forget to cite us:

Barratt, JLN, S Park, FS Nascimento, J Hofstetter, M Plucinski, S Casillas, RS Bradbury, MJ Arrowood, Y Qvarnstrom, E Talundzic (2019) Genotyping genetically heterogeneous Cyclospora cayetanensis infections to complement epidemiological case linkage. Parasitology:1–9 doi:10.1017/S0031182019000581
