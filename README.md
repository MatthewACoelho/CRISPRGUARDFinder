# CRISPR GUARD Finder
## Introduction
Given an on-target CRISPR guide this tool will identify all of its off-targets and design “GUARD” sequences (short guides) to interfere with the off-target activity. This is the code referred to in Coelho et al [1].

The tool incorporates an off-target search based on the one used in the Sanger WGE website [2, 3] enhanced with the calculation of the probability of the off-target [4]. R is used to identify and score the GUARD designs, and the whole is coordinated by a nextflow script.

1. Coelho et al, in preparation
2. [WGE: a CRISPR database for genome engineering.  - PubMed - NCBI](https://www.ncbi.nlm.nih.gov/pubmed/25979474)
3. [WTSI Genome Editing](https://www.sanger.ac.uk/htgt/wge/)
4. [Repurposing CRISPR as an RNA-guided platform for sequence-specific control of gene expression.  - PubMed - NCBI](https://www.ncbi.nlm.nih.gov/pubmed/23452860)

## Installation
Create a new directory for the installation which I will denote `$guard_root` in what follows.

Clone out the GitHub repo:

```
cd $guard_root

git clone https://github.com/MatthewACoelho/GUARDfinder.git .
```

The pipeline requires R with `optparse`  and `BSgenome` packages for each of the genomes required, and `nextflow`. One way to get these if you don’t already have them is using `conda`:

```
conda create --prefix $guard_root/ext
source activate $guard_root/ext

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

## Genomes
The genome-related files and indexes are too large to include here, so before you can run the tool you will need to generate them.

Download the genome fasta file and GTF from Ensembl for each organism in which you are interested. For example Human:

```
curl http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz > $guard_root/data/hg38.fa.gz
curl http://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.99.chr.gtf.gz > $guard_root/data/hg38.gtf.gz
```

and Mouse:

```
curl http://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz > $guard_root/data/mm10.fa.gz
curl http://ftp.ensembl.org/pub/current_gtf/mus_musculus/Mus_musculus.GRCm38.96.gtf.gz > $guard_root/data/mm10.gtf.gz
```

To create an off-target index for a genome `$genome` for PAM `$pam`:

```
mkdir -p $guard_root/data/$genome/$pam
gunzip $guard_root/data/$genome.fa.gz
$guard_root/bin/ot index -out $guard_root/data/$genome/$pam -pam $pam $guard_root/data/$genome.fa
```

You can remove some of the indexes to reduce noise in the results, we normally only keep the major chromosome indexes:

```
rm $guard_root/data/$genome/$pam/chrCHR*
rm $guard_root/data/$genome/$pam/chrGL*
rm $guard_root/data/$genome/$pam/chrK*
```

To convert an Ensembl GTF into the format required by the off-target search (which requires perl):

```
mkdir -p $guard_root/data/$genome/info
gunzip $guard_root/data/$genome.gtf.gz
$guard_root/bin/process_gtf.pl $guard_root/data/$genome/info < $guard_root/data/$genome.gtf
```

## Usage
Activate the environment, if required:

```
source activate $guard_root/ext
```

make sure an environment variable `guard_root` is set appropriately and R and nextflow are in your `PATH`.

```
export guard_root = $PATH/guard_root
```

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

	/* Maximum distance from a designed guard to its off-target guide */
	max_guard_distance = 10
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
* OffGuideStrand: Strand of the off-target guide
* GuardToOffGuidePAMDistance: base distance between the guard and off-target guide PAMs
* pCoding, nCoding: p-value and number of off-targets of the GUARD hitting coding (CDS) regions
* pNonCoding, nNonCoding: p-value and number of off-targets of the GUARD hitting non-coding regions
* nBadSeed: number of off-targets for the GUARD where the seed region is identical to the on-target
* X0, X1, X2, X3: number of genomic hits of the GUARD with 0, 1, 2, 3 mismatches.
* OnTargetHit, OnTargetMismatches, OnTargetSeedMismatches: sequence of the potential on-target hit of the GUARD
* PAMAndSeedScore, OnTargetScore, OffTargetScore, GCScore: components of the final GUARD score

## License
The code is freely available under the [MIT license](http://www.opensource.org/licenses/mit-license.html) .


#ngs/crispr
