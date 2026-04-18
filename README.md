# new_gwes

Tree-weighted conditional mutual information tools for GWES.

This repository currently contains:

- a Python command-line tool to compute:
  - tree-weighted pooled MI
  - tree-weighted conditional MI
  - delta MI
  - monostate / informative weight mass
- an R plotting script for visualizing the output TSV

## Requirements

### Python

* Python >= 3.10
* `pip`
* recommended: a virtual environment or conda environment

### Python package dependencies

Installed automatically by `pip`:

* `numpy`
* `biopython`

### R

The plotting script requires:

* `Rscript`
* R package:

  * `data.table`

Optional R packages for improved plotting:

* `ggplot2`
* `hexbin`
* `ggrastr`
* `ggthemes`

If these optional packages are not available, the plotting script falls back to base R plotting.

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/jeju2486/cluster-GWES.git
cd cluster-GWES
```

### 2. Install in editable mode

Editable install is recommended during development.

```bash
pip install -e .
```

This installs the Python package and exposes the command:

```bash
tree-weighted-conditional-mi
```

### 3. Or install as a regular package

```bash
pip install .
```

---

## Verify installation

Check that the command is available:

```bash
tree-weighted-conditional-mi --help
```

You should see the CLI help for the MI calculation script.

---

## Python command-line tool

### Purpose

`tree-weighted-conditional-mi` computes, for each candidate locus pair:

* `mi_tree_pool`
* `mi_tree_cond`
* `delta_mi`
* `mono_weight_mass`
* `informative_weight_mass`

using:

* a Newick tree with branch lengths
* a whitespace-delimited pair file
* a fake FASTA with:

  * `A = absence = 0`
  * `C = presence = 1`

### Main inputs

#### Tree file

* Newick format
* branch lengths required

#### Pair file

Whitespace-delimited text file with at least:

* column 1: `u`
* column 2: `v`
* column 3: `distance`

Additional columns may be carried through to output.

#### Fake FASTA

Requirements:

* FASTA headers must match tree tip names exactly
* all sequences must have equal length
* only `A` and `C` are allowed

  * `A = 0`
  * `C = 1`

---

## Example usage

### Basic serial run

```bash
tree-weighted-conditional-mi \
  --tree /path/to/tree.nwk \
  --pairs /path/to/pairs.tsv \
  --fasta /path/to/fake.fasta \
  --n-min 50 \
  --L-min 0.02 \
  --out /path/to/results.tsv
```

### Parallel run

```bash
tree-weighted-conditional-mi \
  --tree /path/to/tree.nwk \
  --pairs /path/to/pairs.tsv \
  --fasta /path/to/fake.fasta \
  --n-min 50 \
  --L-min 0.02 \
  --chunk-size 5000 \
  --jobs 8 \
  --out /path/to/results.tsv
```

### Without midpoint rooting

```bash
tree-weighted-conditional-mi \
  --tree /path/to/tree.nwk \
  --pairs /path/to/pairs.tsv \
  --fasta /path/to/fake.fasta \
  --n-min 50 \
  --L-min 0.02 \
  --no-midpoint-root \
  --out /path/to/results.tsv
```

### Drop extra pair columns from output

```bash
tree-weighted-conditional-mi \
  --tree /path/to/tree.nwk \
  --pairs /path/to/pairs.tsv \
  --fasta /path/to/fake.fasta \
  --n-min 50 \
  --L-min 0.02 \
  --drop-extra-cols \
  --out /path/to/results.tsv
```

---

## Command-line options

```text
--tree              Newick tree file with branch lengths
--pairs             Whitespace-delimited pair file
--fasta             Fake FASTA with A=0 and C=1
--n-min             Minimum cluster size
--L-min             Minimum relative subtree length in [0,1]
--out               Output TSV path
--alpha             Pseudocount alpha (default: 0.5)
--chunk-size        Number of pair rows per batch (default: 2000)
--jobs              Number of worker processes (default: 1)
--no-midpoint-root  Do not midpoint-root unrooted tree
--drop-extra-cols   Do not carry through pair columns after distance
```

---

## Output files

### Main output

The main TSV contains one row per pair.

Minimum output columns:

```text
u
v
distance
mi_tree_pool
mi_tree_cond
delta_mi
mono_weight_mass
informative_weight_mass
```

If extra pair columns are preserved, they appear between `distance` and the MI columns.

### Sidecar files

In addition to the main output, the script writes:

#### Tip weights

```text
<out>.tip_weights.tsv
```

Columns:

```text
tip_index
tip_name
tip_weight
cluster_id
```

#### Cluster summary

```text
<out>.cluster_summary.tsv
```

Columns:

```text
cluster_id
n_tips
effective_mass
cluster_weight
```

#### Cluster members

```text
<out>.cluster_members.tsv
```

Columns:

```text
cluster_id
tip_index
tip_name
```

---

## Plotting

The repository also includes an R plotting script:

```text
new_gwes/tree_mi_plotting.r
```

This script is not currently installed as a separate shell command by `pip`, so run it directly with `Rscript` from the repository checkout.

### Example

```bash
Rscript new_gwes/tree_mi_plotting.r \
  -i /path/to/results.tsv \
  -o /path/to/out.png
```

### Example with custom y column

```bash
Rscript new_gwes/tree_mi_plotting.r \
  -i /path/to/results.tsv \
  -o /path/to/out.png \
  -y delta_mi
```

### Example with distance filters and downsampling

```bash
Rscript new_gwes/tree_mi_plotting.r \
  -i /path/to/results.tsv \
  -o /path/to/out.png \
  --min-dist 1000 \
  --max-points 500000 \
  --seed 1
```

### Plotting script options

```text
-i, --input            Output TSV from tree_weighted_conditional_mi.py
-o, --output           Output PNG
-y, --y-col            y column (default: delta_mi)
--n-rows               Read first N rows
-l, --ld-dist          LD distance line
--ld-dist-alt          Second LD distance line
--min-dist             Filter distance >= min-dist
--max-points           Downsample background only
--seed                 RNG seed
--inblock              Prefix/dir for per-block unitig files
--mask-cache           RDS cache path
--no-deps              Force base plotting
--true-block-ids       Comma-separated block ids to highlight pairwise
```

---

## Notes on performance

* The Python script uses bit-packing and vectorized lookup tables for speed.
* Best parallel behavior is usually on Linux.
* Speedup is typically sublinear because runtime may become memory-bandwidth limited.

---

## Development notes

Reinstall after packaging changes if needed:

```bash
pip install -e . --upgrade
```

Remove old build artifacts if you later add packaging outputs:

```bash
rm -rf build dist *.egg-info
```

---

## License

See `LICENSE`.


