NP-SMLR is a software for GpC methylation detection and nucleosome occupancy inference at single-molecule level. It takes the aligned signal-levels from Oxford Nanopore sequencing reads, and outputs the posterior probability of methylation on every GpC site together with the positions of all inferred nuclosomes on each molecule.


# Dependencies

The current version has been tested in the Linux system with:
1. GCC (version 7.3.0);
2. Perl (version 5.16.3);
3. SAMtools (version 1.9);
4. BEDtools (version 2.28.0);
5. Boost (version 1.73.0).


# Installation

You can download and compile the latest version (v1.0) as follows:

```
git clone https://github.com/imatrm/NP-SMLR.git
cd NP-SMLR
make
```

After compilation, three executable files (`Likelihood`, `Detection` and `NclsPos`) will be generated in the folder `NP-SMLR/bin`.


**Please ensure that you have added**

1. `SAMtools`
2. `BEDtools`

**to your `PATH`.**



# Before running NP-SMLR

Before running NP-SMLR,
1. the reads have be aligned to reference genome; and
2. the alignments between event signals and reference genome should be obtained using `nanopolish eventalign`. Please make sure that the flag --print-names is set.

The tutorial of `nanopolish eventalign` can be found at https://nanopolish.readthedocs.io/en/latest/quickstart_eventalign.html.



# Usage

The command is

```
./NP-SMLR.pl -b sorted_bam_file -e event_align_scale -o output_dir
```

* Parameters

```
-b   <STRING>   Name of the BAM file that records the alignment of reads. The BAM file must be sorted.
-g   <STRING>   Event alignment file generated by "nanopolish eventalign" (with flag --scale-events).
-t   <STRING>   Name of output folder.
```

* Example

The test example can be run using the command

```
./NP-SMLR.pl -b ./testdata/sort_test.bam -e ./testdata/eventalign_scale_test.txt -o output
```

The generated files are stored under `output`.



# Output

The file `ncls_pos.bed` in the output folder records the coordinates of nucleosomes detected on each molecule. The columns represent:

Column 1: Chromosome ID;

Column 2: Start position of nucleosome;

Column 3: End position of nucleosome;

Column 4: Name of molecule (read);

Column 5: Quality score of alignment (obtained from BAM file);

Column 6: Strand on which the read is aligned.


The file `detection.txt` in the output folder records the methylation score of every GpC site on each molecule, which is essentially the posterior probability of methylation. The columns represent:

Column 1: Chromosome ID;

Column 2: Position of GpC site;

Column 3: Nucleotide sequence covering the GpC site (the 6-mers on this sequence are involved in the detection of GpC methylation);

Column 4: Name of molecule (read);

Column 5: Log likelihood of event-level with respect to negative control;

Column 6: Log likelihood of event-level with respect to positive control;

Column 7: Methylation score of GpC site (the posterior probability of methylation).



# Current version

The version of the current release is v1.0.



# Contact

Please contact wanganqi18@gmail.com for any question.



# License

**Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public License**

For details, please read `NP-SMLR/License`.



# Citation

Yunhao Wang*, Anqi Wang*, Zujun Liu, Andrew L. Thurman, Linda S. Powers, Meng Zou, Yue Zhao, Adam Hefel, Yunyi Li, Joseph Zabner, Kin Fai Au. Single-molecule long-read sequencing reveals the chromatin basis of gene expression. *Genome Research*, 2019. doi: 10.1101/gr.251116.119 (*contributed equally)
