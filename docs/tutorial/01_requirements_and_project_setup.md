# 01. Requirements, beginner installation, and project setup

[← Previous: Overview](00_overview.md) · [Tutorial contents](README.md) · [Next: Samples and metadata →](02_samples_and_metadata.md)

## Recommended environment

The copy-and-adapt commands in this tutorial target Bash on GNU/Linux, a
managed Linux HPC system, or Windows Subsystem for Linux (WSL). Bioconda also
supports macOS, but some system utilities in the command blocks need macOS-
specific equivalents. Bioconda does not support native Windows. Windows
beginners should normally use WSL unless their institution provides a suitable
Linux server.

Do not install software on a managed university or institutional system unless
the administrator or tutor permits it. Such systems often provide approved
modules or environments that should be used instead.

Core tools used by the tutorial are FastQC, MultiQC, Cutadapt, Bowtie2,
SAMtools, Picard, bedtools, MACS3, featureCounts, R, DESeq2, and apeglm. See the
[recorded case-study versions](../../environment/software_versions.tsv). The
beginner environment below installs compatible available versions; it does not
reconstruct the historical internship environment byte for byte.

## Beginner route: Windows with WSL

Skip this section if Linux, macOS, or a managed HPC environment is already
available.

Open **PowerShell as Administrator**, paste the following command, and restart
Windows when prompted:

```powershell
wsl --install
```

Open the installed Ubuntu application and create the requested Linux username
and password. All remaining commands in this tutorial belong in the Ubuntu/WSL
terminal, not PowerShell. Microsoft's current WSL instructions are available
in the [official WSL installation guide](https://learn.microsoft.com/windows/wsl/install).

Install the two basic download tools inside Ubuntu/WSL:

```bash
sudo apt update
sudo apt install -y curl git
```

## Install Miniforge on Linux, macOS, or WSL

If `conda --version` or `mamba --version` already works, skip this section.
Otherwise, paste the following commands into the Bash terminal:

```bash
set -euo pipefail

curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash "Miniforge3-$(uname)-$(uname -m).sh"
```

Read the installer prompts, accept the license, keep the suggested user-level
location, and allow shell initialization. Close and reopen the terminal, then
verify the installation:

```bash
conda --version
mamba --version
```

These commands follow the
[official Miniforge instructions](https://github.com/conda-forge/miniforge).

## Download the tutorial repository

Choose a working location with enough free space. FASTQ, BAM, reference, and
temporary files can require tens or hundreds of gigabytes.

```bash
git clone https://github.com/yazalj/bulk-atacseq-mature-memory-nk.git
cd bulk-atacseq-mature-memory-nk
```

If the repository is already downloaded, enter its directory instead and
confirm that the expected file is present:

```bash
test -f environment/tutorial_environment.yml
pwd
```

## Install the tutorial software

Create the isolated environment from the supplied file. Installation may take
several minutes:

```bash
mamba env create --file environment/tutorial_environment.yml
conda activate bulk-atacseq-tutorial
```

The environment file uses `conda-forge` and `bioconda` and excludes the
`defaults` channel. For background and alternative installation methods, see
the [official Bioconda usage guide](https://bioconda.github.io/).

At the start of every future terminal session, reactivate it with:

```bash
conda activate bulk-atacseq-tutorial
```

The environment includes the core workflow tools, DESeq2/apeglm and plotting
packages, plus optional HOMER and g:Profiler support for Chapter 10. The source
file deliberately does not pin the repository to one historical build. Record
the versions that were actually resolved for each new study.

## Verify the installation

Paste this block after activating the environment. It stops at the first
missing command or R package:

```bash
set -euo pipefail

for program in fastqc multiqc cutadapt bowtie2 bowtie2-build samtools \
  picard bedtools macs3 featureCounts R Rscript annotatePeaks.pl \
  findMotifsGenome.pl; do
  command -v "${program}" >/dev/null || {
    echo "Missing program: ${program}" >&2
    exit 1
  }
done

Rscript -e '
required <- c("DESeq2", "apeglm", "vsn", "ggplot2", "pheatmap", "gprofiler2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
cat("Required R packages are available.\n")
'

echo "Environment verification: PASS"
```

Record exact versions before analysis:

```bash
mkdir -p logs

{
  date '+%Y-%m-%dT%H:%M:%S%z'
  fastqc --version
  multiqc --version
  cutadapt --version
  bowtie2 --version | head -n 1
  samtools --version | head -n 1
  picard --version
  bedtools --version
  macs3 --version
  featureCounts -v
  R --version | head -n 1
} > logs/tutorial_software_versions.txt 2>&1

Rscript -e 'sessionInfo()' > logs/tutorial_R_sessionInfo.txt
```

If only the R packages are missing from an otherwise approved R installation,
use Bioconductor's supported installer from inside R:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("DESeq2", "apeglm", "vsn"))
install.packages(c("ggplot2", "pheatmap", "gprofiler2"))
BiocManager::valid()
```

Do not mix this fallback with a managed institutional R library unless its
administrator authorizes package installation. Bioconductor recommends
`BiocManager::install()` because R and Bioconductor releases must remain
compatible; see the [official installation page](https://bioconductor.org/install/).

## Suggested working layout

Keep large data outside Git even when they live inside the analysis directory:

```text
project/
├── config/       metadata, input manifest, and parameters
├── workflow/     scripts or copied command records under version control
├── raw/          immutable reads, outside Git
├── reference/    FASTA, indexes, annotation, and blacklist, outside Git
├── results/      stage-specific, non-overwriting outputs
├── logs/         commands, versions, warnings, and validation summaries
└── docs/         decisions, methods, limitations, and report material
```

Create the working directories and copy the blank configuration templates:

```bash
set -euo pipefail

mkdir -p config/local raw reference results logs docs
cp config/templates/samples.template.tsv config/local/samples.tsv
cp config/templates/input_manifest.template.tsv config/local/input_manifest.tsv
cp config/templates/analysis_parameters.template.yaml config/local/analysis_parameters.yaml
```

Edit all three copies and replace every example row, placeholder path,
`REVIEW_ME`, and `REVIEW_FROM_*` value. Then check for missed placeholders:

```bash
if grep -R -n -E 'REVIEW_ME|REVIEW_FROM_|/path/to/' config/local; then
  echo "Unresolved configuration placeholders remain." >&2
  exit 1
else
  echo "Configuration placeholder check: PASS"
fi
```

## Resources that software installation does not provide

The environment does **not** choose or download the biological resources for a
new study. Before alignment, obtain and record:

- the authoritative reference FASTA and its checksum;
- the Bowtie2 index built from exactly that FASTA;
- a compatible gene annotation release;
- a compatible blacklist if one exists and is appropriate;
- confirmed library-specific adapter sequences; and
- the effective genome size required by the chosen peak-calling setup.

Do not copy the internship's GRCh38 paths or thresholds into another organism
or assembly without review.

## Setup gate

Before analysis, confirm that:

- every required executable resolves inside the intended environment;
- versions and the R session are recorded;
- compute limits are known and explicitly set;
- raw and reference inputs are read-only;
- outputs cannot overwrite inputs or previous attempts;
- available disk space is suitable for FASTQ, BAM, and temporary files;
- all configuration placeholders have been replaced; and
- secrets, controlled data, and private absolute paths will not enter Git.

Continue to [samples, metadata, and statistical design](02_samples_and_metadata.md)
only after this gate passes.

---

[← Previous: Overview](00_overview.md) · [Tutorial contents](README.md) · [Next: Samples and metadata →](02_samples_and_metadata.md)
