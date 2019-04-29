# CRISPR GUARD Finder
## Introduction
Given an on-target CRISPR guide this tool will identify all its off-targets and design “GUARD” sequences (short guides) to interfere with the off-target activity. This is the code referred to in Coehlo et al [1].

The tool incorporates and off-target search based on the one used in the Sanger WGE website [2, 3] enhanced with the calculation of the probability of the off-target [4]. R is used to identify and score the guard designs, and the whole is coordinated by a nextflow script.

[1] Coehlo et all, in preparation
[2] [WGE: a CRISPR database for genome engineering.  - PubMed - NCBI](https://www.ncbi.nlm.nih.gov/pubmed/25979474)
[3] [WTSI Genome Editing](https://www.sanger.ac.uk/htgt/wge/)
[4] [Repurposing CRISPR as an RNA-guided platform for sequence-specific control of gene expression.  - PubMed - NCBI](https://www.ncbi.nlm.nih.gov/pubmed/23452860)

## Installation
Create a new directory for the installation which I will denote `$guard_root`  in what follows.

Clone out the GitHub repo:

```
cd $guard_root

github clone https://github.com/MatthewACoelho/GUARDfinder.git
```

The pipeline requires R with `optparse`  and `BSgenome` packages for each of the genomes required, and `nextflow`. One way to get these if you don’t already have them is using `conda`:

```
conda create --prefix `pwd`/ext
source activate `pwd`/ext

conda install -c r r-base
conda install -c bioconda bioconductor-bsgenome.hsapiens.ucsc.hg38
conda install -c bioconda bioconductor-bsgenome.mmusculus.ucsc.mm10
conda install -c bioconda r-optparse
conda install -c bioconda nextflow

source deactivate
```

The off-target search is a C program requiring 64-bit architecture, which can be compiled as follows (only tested on Centos 7).

```
gcc -mcmodel=medium -O4 -o bin/ot -pthread src/ot.c
```

All that remains is to edit `activate.sh` to set the `guard_root` to the full path of the one you have created.

## Usage
Activate the environment, if required:

```
source $guard_root/activate.sh
```

alternatively make sure an environment variable `guard_root` is set appropriately and R and nextflow are in your `PATH`.

Create a directory for the run, and within it a  `params.nf` file like:

```
params {
	/* Name of this run - for file naming */
	id = "VEGFA"

	/* On-target guide sequence without PAM */
	guide = "GACCCCCTCCACCCCGCCTC"

  /* On-target guide location if there is > 1 with 0 mismatches */
	chr = ""        /* chr4 for example */
	start = ""
	end = ""
	strand = ""     /* + or - */

	/* Genome version - hg38 or mm10 */
	genome = "hg38"

	/* PAM sequence - must correspond to one of the off-target indexes */
	pam = "NGG"

	/* Mismatches for the guide off-target search */
	guide_mismatches = 5

	/* Filter for which off-targets to take forward to guard design */
	guide_min_pvalue = 0.1

	/* Desired length of guard */
	guard_length = 14

	/* Mismatches for the guard off-target search */
	guard_mismatches = 3
}
```

In the same directory run:

```
$guard_root/bin/find_guards.sh
```

which will run the process locally, or:

```
$guard_root/bin/find_guards.sh --submit
```

which will submit jobs to slurm.

In both cases the script will pick up cached results if they exist, so add `--force` to force recalculation.

When complete you should have a text file named `id_final.txt`, e.g. `VEGFA_final.txt`.

## Output
The output is a tab-delimited text file with the following columns:

* Guard: sequence of the GUARD
* ID: unique identifier for the off-target this guard relates to - a unique off-target number followed by the gene, it hits if there is one 
* GuardChr, GuardStart, GuardStop, GuardStrand: genomic location and strand of the GUARD
* GuardWithPAM: GUARD with PAM from the + strand view
* GuardGC: GC content of the GUARD
* ForwardGuardPAM: GUARD PAM on the GUARD strand
* ForwardGuardWithPAM: GUARD with PAM on the GUARD strand
* OffGuideOverlap, OffGuideSeedOverlap, OffGuidePAMOverlap: number of bases the GUARD overlaps with the off-target guide_seed region_PAM
* OffGuide: sequence of the off-target guide
* OffGuideGC: GC content of the off-target guide
* pCoding, nCoding: p-value and number of off-targets of the GUARD hitting coding (CDS) regions
* pNonCoding, nNonCoding: p-value and number of off-targets of the GUARD hitting non-coding regions
* nBadSeed: number of off-targets for the GUARD where the seed region is identical to the on-target
* X0, X1, X2, X3: number of genomic hits of the GUARD with 0, 1, 2, 3 mismatches.
* OnTargetHit, OnTargetMismatches, OnTargetSeedMismatches: sequence of the potential on-target hit of the GUARD
* PAMScore, SeedScore, OnTargetScore, OffTargetScore, GCScore: components of the final GUARD score
	* Score: final GUARD score

## Additional Genomes
To create an off-target index for a genome `$genome` for PAM `$pam`:

```
mkdir $guard_root/$genome/$pam
$guard_root/bin/ot index -out $guard_root/data/$genome/$pam -pam $pam $genome.fa
```

And to convert an Ensembl GTF into the format required by the off-target search (which requires perl):

```
mkdir $guard_root/$genome/info
$guard_root/bin/process_gtf.pl $guard_root/data/$genome/info < $genome.gtf
```

## License
The code is freely available under the [MIT license](http://www.opensource.org/licenses/mit-license.html) .


#ngs/crispr