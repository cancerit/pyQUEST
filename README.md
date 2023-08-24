# pyQUEST

Count unique reads and optionally match them to a given library (exact matching only).

Input files:

- SAM/BAM/CRAM/FASTQ file
- [library file](#library) (library-dependent mode only)

Output files:

- **library-independent count**:
  - [counts](#library-independent-counts)
  - [statistics](#library-independent-stats-file)
- **library-dependent count**:
  - [counts](#library-dependent-counts)
  - [statistics](#library-dependent-stats-file)

Notes:

- only supports single-sample input files
- reads with ambiguous nucleotides are discarded
- masked reads are discarded

## Setup

Using a Python virtual environment:

```sh
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install .
```

The Docker image can be built as follows:

```sh
docker build -t pyquest .
```

## Usage

```
Usage: pyquest [OPTIONS] QUERIES

  Count reads and optionally map them to a library.

  QUERIES: Query sequence file (fastq[.gz], sam, bam, cram)

Options:
  -o, --output PATH               Final output to this filename prefix
                                  [required]
  --min-length INTEGER RANGE      Minimum read length  [default: 1; x>=1]
  --most-common INTEGER RANGE     Output top X most common unique read
                                  sequences in FASTA format  [1<=x<=50]

Input sample metadata:         Options adding information to the input
    -s, --sample TEXT             Sample name to apply to count column,
                                  required for fastq, interrogate header for
                                  others when not defined.
    -r, --reference FILE          Required for CRAM

Library-dependent:             Options specific to library-dependent
                                  counting
    -l, --library FILE            Expanded library definition TSV file with
                                  optional headers (common format for
                                  single/dual/other)
    --low-count INTEGER RANGE     *.stats.json includes
                                  low_count_guides_lt_{15,30}, this option
                                  allow specification of an additional cut-
                                  off.  [x>=0]

Performance:                   Options to tune the performance
    -c, --cpus INTEGER RANGE      CPUs to use (0 to detect)  [default: 1;
                                  x>=0]

Debug:                         Options specific to troubleshooting, testing
                                  and debugging
    --loglevel [WARNING|INFO|DEBUG]
                                  Set logging verbosity  [default: INFO]
    --no-compression              Disable output compression
  --version                       Show the version and exit.
  --help                          Show this message and exit.
```

With Docker:

```sh
# Output in the current directory
mkdir -p output
docker run \
    -v "$PWD/test.queries.bam":/tmp/x.bam:ro \
    -v "$PWD/output":/output \
    pyquest \
        pyquest \
            -o /output/something \
            --sample XYZ \
            --no-compression \
            /tmp/x.bam
```

## File header formats

TSV headers may contain metadata in the form of key-value pairs thus formatted:

```
##<KEY>: <VALUE>
```

The column headers, separated by tabs, immediately follow the metadata lines and are preceded by a single `#` character, *e.g.*:

```
#<FIELD 1>	<FIELD 2>	<FIELD 3>
```

### Count header

|Field|Format|Description|
|-|-|-|
|`Command`|string|Full command|
|`Version`|`x.y.z`|Tool version|

### Library header

Currently, ignored.

## File formats

### Library

**Format**: TSV with [library header](#library-header)

The headers are ignored, and therefore the relevant fields are identified by their position. Here we indicate the field positions as one-based, with their corresponding field names in the [library-dependent counts](#library-dependent-counts).

|Position|Counts field|Format|Description|
|-|-|-|-|
|1|`ID`|string|Library sequence identifier|
|2|`NAME`|string|Library sequence name|
|3|`SEQUENCE`|`[ACGT]+`|DNA sequence|

*E.g.*:

```
## ...
# ...
1	some-name-1	AAAAAAAAATCCAGAACCT
2	some-name-2	AAAAAAATATGCCCGTGGA
3	some-name-3	AAAAAAGCATTTAGGCAGG
4	some-name-4	AAAAAAGCTTGCATTAGAC
5	some-name-5	AAAAAATATCGTGTCAAGT
6	some-name-6	AAAAAATCAGCCACGCGAC
```

### Library-independent counts

**Format**: TSV with [count header](#count-header) (gzip'ed by default)

|Field|Format|Description|
|-|-|-|
|`SEQUENCE`|`[ACGT]+`|Unique DNA sequence|
|`LENGTH`|integer|Length of the sequence|
|`COUNT`|integer|Number of reads|

*E.g.*:

```
##Command: pyquest -o output --min-length 0 --low-count 2 -l guides.tsv --sample XYZ --no-compression test.queries.bam
##Version: 1.0.0
#SEQUENCE	LENGTH	COUNT
AAAAAAGCTTGCATTAGAC	19	25
AAAAAATATCGTGTCAAGT	19	26
AAAAAATGTCAGTCGAGTG	19	34
AAAAACAAGCGCACCACCG	19	1
AAAAACACTTCCATGCAAA	19	25
AAAAACGTATTTAGCCGAA	19	23
```

### Library-dependent counts

**Format**: TSV with [count header](#count-header) (gzip'ed by default)

|Field|Format|Description|
|-|-|-|
|`ID`|string|Library sequence identifier|
|`NAME`|string|Library sequence name|
|`SEQUENCE`|`[ACGT]+`|DNA sequence|
|`LENGTH`|integer|Length of the DNA sequence|
|`COUNT`|integer|Number of reads|
|`UNIQUE`|0\|1|Whether the sequence is unique in the library|
|`SAMPLE`|string|Name of the sample of origin of the reads|

*E.g.*:

```
##Command: pyquest -o output --min-length 0 --low-count 2 -l guides.tsv --sample XYZ --no-compression test.queries.bam
##Version: 1.0.0
#ID	NAME	SEQUENCE	COUNT	UNIQUE	SAMPLE
1	some-name-1	AAAAAAAAATCCAGAACCT	0	1	XYZ
2	some-name-2	AAAAAAATATGCCCGTGGA	0	1	XYZ
3	some-name-3	AAAAAAGCATTTAGGCAGG	0	1	XYZ
4	some-name-4	AAAAAAGCTTGCATTAGAC	25	1	XYZ
5	some-name-5	AAAAAATATCGTGTCAAGT	26	1	XYZ
6	some-name-6	AAAAAATCAGCCACGCGAC	0	1	XYZ
```

### Library-independent stats file

**Format**: JSON

|Field|Format|Description|
|-|-|-|
|`sample_name`|string|Name of the sample|
|`input_reads`|integer|Total input reads|
|`total_reads`|integer|Total reads passed on to counting|
|`discarded_reads`|integer|Total reads discarded before counting|
|`vendor_failed_reads`|integer|Total reads with the `QCFAIL` flag|
|`length_excluded_reads`|integer|Total reads discarded because shorter than a user-defined threshold|
|`ambiguous_nt_reads`|integer|Total reads with ambiguous nucleotides|
|`masked_reads`|integer|Total soft-masked reads|
|`zero_length_reads`|integer|Total zero-length reads|

*E.g.*:

```json
{
    "version": "1.0.0",
    "command": "pyquest -o output --min-length 0 --low-count 2 -l guides.tsv --sample XYZ --no-compression test.queries.bam",
    "sample_name": "XYZ",
    "total_reads": 1020769,
    "vendor_failed_reads": 0,
    "length_excluded_reads": 0,
    "ambiguous_nt_reads": 0,
    "masked_reads": 0
}
```

### Library-dependent stats file

**Format**: JSON

The library-dependent count statistics include the [library-dependent count statistics](#library-independent-stats-file).

All statistics are computed on the read counts of unique targets, excluding those discarded based on their length. The number of low count templates (`zero_count_templates` and `low_count_templates_*`) also excludes the targets with short sequences.

|Field|Format|Description|
|-|-|-|
|`mapped_to_template_reads`|integer|Total reads mapping to the library|
|`mean_count_per_template`|decimal|Mean reads per template|
|`median_count_per_template`|decimal|Median reads per template|
|`multimap_reads`|integer|Total reads mapping to more than one template|
|`unmapped_reads`|integer|Total reads mapping to no template|
|`total_templates`|integer|Total number of templates|
|`total_unique_templates`|integer|Total number of unique templates|
|`length_excluded_templates`|integer|Total number of unique templates excluded by length|
|`zero_count_templates`|integer|Total number of unique templates with no reads mapping to them|
|`low_count_templates_lt_15`|integer|Total number of unique templates with less than 15 reads mapping to them|
|`low_count_templates_lt_30`|integer|Total number of unique templates with less than 30 reads mapping to them|
|`low_count_templates_user`|object\|`null`|Total number of unique templates with less than a user-defined number of reads mapping to them (optional)|
|`gini_coefficient`|decimal|Gini coefficient of the mapping read counts|

*E.g.*:

```json
{
    "version": "1.0.0",
    "command": "pyquest -o output --min-length 3 --low-count 2 -l guides.tsv --sample XYZ --no-compression test.queries.sam",
    "sample_name": "XYZ",
    "input_reads": 1020770,
    "total_reads": 1020766,
    "discarded_reads": 4,
    "vendor_failed_reads": 0,
    "length_excluded_reads": 1,
    "ambiguous_nt_reads": 2,
    "masked_reads": 2,
    "mapped_to_template_reads": 1020766,
    "mean_count_per_template": 10.1,
    "median_count_per_template": 0,
    "multimap_reads": 0,
    "unmapped_reads": 0,
    "total_templates": 101064,
    "total_unique_templates": 101064,
    "length_excluded_templates": 0,
    "zero_count_templates": 60927,
    "low_count_templates_lt_15": 72265,
    "low_count_templates_lt_30": 84339,
    "low_count_templates_user": {
      "lt": 2,
      "count": 61744
    },
    "gini_coefficient": 0.73
}
```
