# BioinfUitls

**bioinfutils** is a command-line utility toolset designed for various bioinformatics tasks. It provides functionalities to manipulate genetic sequence data stored in FASTA and GenBank formats, perform sequence filtering, concatenation, modification, and analysis.

## Features

- **Modular Structure**: The toolset is modularized into individual scripts, each targeting specific bioinformatics tasks.
- **Command-Line Interface**: Supports a CLI interface for ease of use and integration into automated workflows.
- **Extensible**: Additional functionalities can be easily added by including new scripts and updating the CLI setup.
- **Dependencies**: Uses essential bioinformatics libraries such as Biopython for sequence manipulation and Pandas for data handling.

## Installation

Install the **bioinfutils** package using the following steps:

1. **From PIP**
   ```
   UNDER DEVELOPMENT pip install .
   ```

2. **From Source**

   ```
   git clone https://github.com/your/repository.git
   cd repository
   ```

## How to run
After installation, you can use bioinfutils command followed by one of the available subcommands:
   ```
   bioinfutils <command> [<args>]
   ```


Available commands include:

- `change_gbk_origin`: Modify GenBank files.
- `filter_fasta`: Filter sequences in FASTA files.
- `filter_genbank_contigs`: Filter GenBank files based on contigs.
- `find_gbk`: Search for specific entries in GenBank files.
- `filter_var_sites`: Filter variable sites in genetic sequences.
- `gbk_reverse_complement`: Reverse complement sequences in GenBank format.
- `gbk2fasta`: Convert GenBank files to FASTA format.
- `gene_pa_matrix`: Generate a matrix of gene presence/absence.
- `get_fasta_stat`: Compute statistics for sequences in FASTA files.
- `multi_fasta_concatenator`: Concatenate sequences from multiple FASTA files.
- `remove_from_fasta`: Remove sequences from a FASTA file based on IDs.
- `split_fasta_desc`: Split FASTA sequences based on description.
- `split_fasta_seq`: Split FASTA sequences into individual files.
- `update_fasta_seqid`: Update FASTA sequence IDs.

Replace <command> with any of the listed commands and provide necessary <args> as described in each command's usage.


## Getting Help
For a particular command, help can be elicited from command line using `--help` or `-h` in the following manner:
```
bioinfutils <command> --help
```

## What it can do

### Fasta manipulation
- Filter fasta
- Remove from fasta
- Split fasta desc
- split fasta seq
- update fasta seqid

### Fasta Utilities
-  Get fasta stat
-  Multi fasta concatenator
-  Filter variable sites
-  Gene presence absense matrix

### Genbank file manipulation
- Filter genbank
- Find genbank 
- Convert to fasta
- Change genbank origin
- Reverse complement Genbank file
 
# Detail

## Change genbank originn
This script modifies the origin coordinate of a circular DNA sequence in a GenBank file. It allows users to rotate the origin of the sequence and updates the positions of associated features (such as genes) accordingly, ensuring that the annotations remain consistent with the new origin.

### Features

- **Change sequence origin**: Repositions the origin of a circular DNA sequence in a GenBank file.
- **Update feature locations**: Adjusts the genomic features (e.g., genes) to match the new sequence origin.
- **Custom output naming**: Automatically generates an output filename based on the new origin or allows users to specify a custom output file.

### Usage
The script can be run from the command line and requires three arguments:
- `--origin` or `-n` (required): The new origin's 1-offset location.
- `input` or `-i` (required): The path to the input GenBank file.
- `output` or `-o` (optional): The path to the output GenBank file. If not provided, the script will generate a default filename based on the new origin.

### Examples

**Change origin and specify output:**

   ```
   bioinfutils change_gbk_origin --origin 23304 --input test/dataset/pTi_GV3101.gbk --output test/dataset/pTi_GV3101_repC.gbk
   ```



## Filter Fasta
This script filters sequences from a FASTA file based on different criteria such as sequence IDs or sequence length. It uses Biopython's SeqIO module to parse and manipulate FASTA files. The user can specify the filtering option through command-line arguments.

### Features
- **Filter by ID list** Selects sequences based on a list of sequence IDs.
- **Filter by size** Retains only the longest sequences when duplicates exist.
- **Filter by sequence** IDs from a file: Filters sequences by matching IDs from a separate file.

### Usage

The program accepts several arguments:

- `-f`, `--fasta` (required): Path to the input FASTA file.
- `-i`, `--ids` (optional): Path to the text file containing sequence IDs (used for filtering by ID).
- `-o`, `--output` (required): Path to the output FASTA file where filtered sequences will be saved.
- `-n`, `--option` (required): Type of filter to apply. Use:
      `1` to filter by sequence IDs from a list.
      `2` to filter by sequence length (keep longest sequences).
      `3` to filter by sequence IDs from a separate file.
- `-l`, `--list` (optional): List of sequence IDs for filtering (used with option 1).

### Example Commands

Filter by ID file:

```bioinfutils filter_fasta --fasta path/to/fna --output path/to/output -option 1 ```

This will filter sequences from input.fasta based on the sequence IDs in ids.txt and write the result to output.fasta.

Filter by size:

```bioinfutils filter_fasta -f input.fasta -o output.fasta -n 2```

This will filter duplicate sequences in input.fasta, keeping only the longest ones.

Filter by list of IDs:

```bioinfutils filter_fasta  input.fasta -l id_list.txt -o output.fasta -n 1```

This will filter sequences in input.fasta based on the IDs listed in id_list.txt.


## Remove from fasta

This script removes specific sequences from a FASTA file based on a provided list of sequence IDs. It reads the input FASTA file, compares each sequence ID against a list of IDs to remove, and writes a new FASTA file excluding the specified sequences.

### Features

- **Sequence Removal**: Removes sequences from a FASTA file based on a list of sequence IDs.
- **Input Flexibility**: Accepts paths to the input FASTA file, a file containing sequence IDs to remove, and specifies the output file path.
- **Output Creation**: Generates a new FASTA file with sequences removed, preserving the original FASTA file.


### Arguments
- `--fasta` or `-f` Path to the input FASTA file from which sequences will be removed.
- `--remove_id` or `-r`  Path to a text file containing sequence IDs (one ID per line) that should be removed from the input FASTA file.
- `--output` or `-o` Path to the output FASTA file where the filtered sequences will be saved.

### Usage

Run the script from the command line with the following arguments:

```bash
bioinfutils sequence_to_remove --fasta input.fasta --remove_id sequences_to_remove.txt --output output.fasta
```

## Split fasta desc

This script processes multiple FASTA files in a directory, updating the sequence IDs based on a specified split symbol in the description. It retains only the portion of the ID up to a specified position after splitting, ensuring unique IDs for each sequence. The updated sequences are then saved to new FASTA files with a user-defined suffix.

### Features

- **Directory Processing**: Processes all FASTA files in a specified directory.
- **ID Update**: Updates sequence IDs based on a split symbol in the description.
- **Longest Sequence Retention**: Ensures that only the longest sequence is kept for each unique ID if duplicates exist.
- **Output File Creation**: Saves updated sequences to new FASTA files with a specified suffix.


### Arguments
- `--input_dir`, `-i`: Path to the directory containing the input FASTA files.
- `--fasta_suffix`, `-t`: Suffix of the FASTA files to be processed (e.g., .fasta).
- `--split_symbol`, `-s`: Symbol used to split the FASTA descriptions to update the IDs.
- `--position_to_keep`, `-p`: Position to keep after splitting the description to form the updated ID.
- `--updated_fasta_suffix`, `-u`: Suffix for the updated FASTA output files.

### Usage
```bash
bioinfutils split_fasta_desc --input_dir path/to/fasta_directory --fasta_suffix .fasta --split_symbol "_" --position_to_keep 2 --updated_fasta_suffix _updated.fasta
```


##  Multi-Seq Fasta Splitter `split_fasta_seq`

This script splits a single FASTA file containing multiple sequences into individual FASTA files. Each sequence from the input file is saved as a separate FASTA file in the specified output directory, optionally with a user-defined prefix added to each filename.

### Features

- **FASTA Parsing**: Reads sequences from a single input FASTA file.
- **File Output**: Writes each sequence to a separate FASTA file in the specified output directory.
- **Prefix Option**: Allows adding a prefix to the filenames of the output FASTA files (optional).

### Arguments
- `-i`,` --input`: Path to the input FASTA file containing multiple sequences.
- `-o`, `--output`: Path to the output directory where individual FASTA files will be saved.
- `-p`, `--prefix` (optional): Prefix to prepend to the filenames of the output FASTA files. If not specified, filenames will start with the sequence IDs.

### Usage

```
bioinfutils split_fasta_seq -i input_file.fasta -o output_directory [-p prefix]
```

## Get Fasta Stat

This script calculates the total genome size for each FASTA file in a specified directory. It reads all FASTA files with a given suffix (e.g., `.fas`), computes the cumulative size of sequences within each file, and outputs the results to a CSV file (`fasta_stat.csv`).

### Features

- **Directory Parsing**: Iterates through all files in a specified directory.
- **FASTA Processing**: Parses each FASTA file to compute the cumulative size of sequences.
- **Output**: Writes the computed genome sizes to a CSV file (`fasta_stat.csv`).

### Arguments
- `--input` or `-i` Path to directory containing fasta files
- `--suffix` or `-s` Suffix of the fasta filename, i.e. fasta, fna etc.

### Usage
```
bioinfutils get_fasta_stat --input path/to/input/directory --suffix .fasta
```

## Update Fasta ID
This script updates FASTA sequence IDs by appending the filename prefix to each ID. It takes an input FASTA file, reads each sequence record, modifies the ID by adding the filename prefix, and writes the updated sequences to an output FASTA file.

### Features

- **ID Update**: Modifies each sequence ID in the input FASTA file by prefixing it with the filename (without extension).
- **Input and Output**: Handles paths to input and output FASTA files specified via command line arguments.
- **Output Confirmation**: Prints a confirmation message displaying the path to the updated output FASTA file.

### Arguments
- `--input` or `-i` Path to input Fasta files
- `--output` or `-o` Output fasta file path


### Usage

```
bioinfutils update_fasta_seqid -i input.fasta -o output.fasta
```

## Multi Fasta Concatenator

This script concatenates sequences from multiple FASTA files, each representing genes from individual strains, into a single FASTA file per strain. It ensures that the order of genes is maintained across all strains, filling gaps for missing genes to maintain alignment integrity.

### Features

- **Input Handling**: Parses multiple FASTA files from a specified input directory.
- **Sequence Concatenation**: Concatenates gene sequences for each strain, ensuring alignment consistency.
- **Output**: Writes concatenated sequences for each strain into a single output FASTA file.
- **Gap Handling**: Inserts gaps for missing genes to maintain alignment across strains.

### Arguments
- `--input` or `-i`: Input directory containing fasta file
- `--output` or `-o`: Output fasta file name with or without path

### Usage

```
bioinfutils multi_fasta_concatenator input_directory output_file.fasta
```

## Filter Variable Sites

This script detects and removes non-synonymous sites from nucleotide codon alignments based on amino acid alignments. It reads input files, identifies variable sites in the amino acid alignment, maps these to codon positions, and produces modified nucleotide alignments with non-synonymous sites removed.

### Features

- **Variable site detection**: Identifies positions in amino acid alignments where there are non-synonymous changes.
- **Conversion of positions**: Converts positions of amino acid changes to nucleotide codon positions for removal.
- **Alignment modification**: Produces new nucleotide codon alignments with non-synonymous sites removed.
- **Batch processing**: Supports processing of multiple strains or genes, matching amino acid and nucleotide alignment files by filename prefix.

### Arguments

The script can be run from the command line with various options:
- `--input_aa_fasta`: Path to amino acid alignment file in FASTA format.
- `--input_na_fasta`: Path to nucleotide codon alignment file in FASTA format.
- `--output_fasta`: Path to save modified codon alignment output file.
- `--input_aa_dir`: Directory path containing multiple amino acid alignment files.
- `--input_na_dir`: Directory path containing multiple nucleotide codon alignment files.
- `--input_aa_suffix`: Common suffix of amino acid alignment files.
- `--input_na_suffix`: Common suffix of nucleotide codon alignment files.

## Gene Presence Absence Matrix 

This script processes multiple FASTA files containing gene sequences from isolates and generates a matrix of presence or size information for each gene across isolates. It outputs this matrix as a CSV file.

### Features

- **Directory scanning**: Automatically identifies and processes all FASTA files in a specified directory.
- **Isolate presence information**: Generates a matrix indicating presence (+) or absence (-) of genes across isolates when option 'a' is used.
- **Gene size information**: Generates a matrix indicating gene lengths across isolates when option 'b' is used.
- **Output**: Saves the resulting matrix as a CSV file named `pa_mat.csv`.

### Arguments

The script is run from the command line and accepts the following arguments:
- `--option` or `-o`: Option 'a' to generate presence matrix or 'b' to generate size matrix.
- `--input` or `-i`: Directory path containing the FASTA files.
- `--suffix` or `-s`: Suffix of the FASTA files (e.g., `.fasta`, `.fa`).

## Filter GenBank

This script filters a GenBank file based on specified contig IDs. It reads a GenBank file, extracts records that match the provided contig IDs, and writes them to a new GenBank file.

### Features

- **Contig filtering**: Filters GenBank records by a list of contig IDs.
- **Input and output customization**: Accepts paths for both input and output GenBank files.
- **Flexible contig ID input**: Allows multiple contig IDs to be specified for filtering.

### Usage

The script is intended to be run from the command line with the following arguments:
- `-i`, `--input`: Path to the input GenBank file.
- `-o`, `--output`: Path to the output GenBank file where filtered records will be saved.
- `-c`, `--contigs`: List of contig IDs (space-separated) to filter by.

## Find in genbank

This script checks for the presence of specified genes in a GenBank file. It reads a list of gene names from a text file and searches through the features of the GenBank file to identify matches.

### Features

- **Gene presence check**: Identifies genes present in a GenBank file based on a provided list.
- **Flexible input**: Accepts any GenBank file and a text file listing gene names.
- **Output**: Prints the GenBank file name and the names of genes found within it.

### Usage

The script is designed to be run from the command line, requiring two arguments:
- `-g`, `--genbank`: Path to the input GenBank file.
- `-l`, `--gene_list`: Path to a text file containing a list of gene names, with each gene name on a separate line.


## gbk2fasta

This script converts GenBank files to FASTA format, offering three output options: whole contigs, CDS nucleotide sequences, or translated CDS amino acid sequences. It can process a single GenBank file or multiple files within a directory. The script is designed for command-line use, enabling efficient batch processing of genomic data.

### Features

- **Convert whole contigs**: Extract entire nucleotide sequences from a GenBank file and convert them to FASTA (`.fna`).
- **Extract CDS nucleotide sequences**: Retrieve CDS features and output them as nucleotide sequences in FASTA format (`.ffn`).
- **Translate CDS to amino acids**: Extract CDS features and translate them into amino acid sequences in FASTA format (`.faa`).
- **Batch processing**: Process a single file or all GenBank files in a directory.

### Usage

The script requires an input GenBank file or directory, an output format, and a directory to save the resulting files. The output format can be one of the following:
- `fna`: for whole contig nucleotide sequences
- `ffn`: for CDS nucleotide sequences
- `faa`: for translated CDS amino acid sequences

### Command-line arguments:
- `--input` or `-i`: Path to a GenBank file or directory containing GenBank files.
- `--output_type` or `-t`: Type of output (`fna`, `ffn`, or `faa`).
- `--output_dir` or `-o`: Directory to save the output files.




## Reverse complement Genbank file

This script processes GenBank files by generating the reverse complement of both the sequences and their associated features (e.g., genes, regulatory elements). It takes an input GenBank file, reverses the sequence and its annotations, and writes the reversed GenBank file to an output location.

### Features

- **Reverse complement sequences**: Reverse complements the nucleotide sequence from a GenBank file.
- **Reverse complement features**: Adjusts the genomic feature annotations (such as genes) to match the reversed sequence.
- **Supports annotations**: Maintains essential annotations like the molecule type in the output.

### Usage

The script can be run from the command line and requires two arguments:
- `--input` or `-i`: The path to the input GenBank file.
- `--output` or `-o`: The path to the output GenBank file where the reverse complemented sequences and features will be saved.