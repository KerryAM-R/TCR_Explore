# Quality control process of T cell receptor from nested PCR experiments.

## Paired Sanger sequencing QC process

Step 1. Convert the .seq to .fasta files (QC -> SEQ to FASTA file merger)
- add in the required file naming conversion e.g. Individual.groupChain-initialwell
- 
- If more identifiers are required for your analysis, they can be added in later. 

Step 2. Upload the .fasta file to [IMGT Vquest](https://www.imgt.org/IMGT_vquest/input)
- select the species (e.g. Homo Sapiens) and receptor type or locus (e.g. TR)
- Upload .fasta file to chose file
- go to C. Excel file -> unchecked all -> check 1. (summary) and 6. (Junction)
  + If you download all 12 tab, move the Junction tab to position 2
- Download vquest.xls file by selecting the start button
- Repeat this process for every .fasta file

Step 3. Downloading TCR_Explore QC file
- Upload the vquest.xls file to QC -> IMGT (Sanger sequencing) -> tab 1. Create QC file. 
- Download the .csv file called "IMGT_onlyQC.csv"
- The program extracts the required columns needed either for the filtering QC process, TCR_Explore or TCRdist

***NOTE:*** *If steps 1 and 2 were completed prior to this process and the Junction sheet is missing, the following process can be completed with __'Summary'__ sheet only. However, without this sheet, one cannot download the TCRdist3 file in the 'TCR analysis -> overview of TCR pairing -> summary table', as the JUNCTION nucleotide sequences will be missing.* 

Step 4. Fill in the QC file
- Copy all sequences into the one .csv file, if more than one plate is present add a number infront of the intial well
- The program adds three columns to the end of the file
  - 'V.sequencing.quality.check' 
    - Will list the following issues: No alignment, Unproductive issue, V identity issue, J identity issue or No issue flagged by IMGT
  - 'clone_quality'
    - If V.sequencing.quality.check reported 'No issue flagged by IMGT' the program pre-fills in as a pass, as the test-data highlighted all 151 sequences had high quality chromatograms
    - The remaining rows in will be filled in as `NA`
    - The `NA` will need to be filled in based on the chromatogram quality which can be check in QC -> Check .ab1 files
    - We recommend checking the quality of the chromatogram from the pre-filled in 'pass' sequences
  - 'comments'
    - Fill in based on `NA` e.g. High quality sequence (pass), two sequence (fail), cannot resolve frameshift/stop sequence (fail), messy sequences (fail)


Step 5. Creating the file needed for paired TCR_Explore
- Upload the filled in IMGT_onlyQC.csv to QC -> IMGT (Sanger sequencing) -> tab 2. Paired chain file -> Completed QC file (.csv)
- Select if the data is alpha-beta (ab) or gamma-delta (gd)
- Select if the data was created from either "summary+JUNCTION" or "Summary"
- This will render two tables in this section:
  - Table 1. summary of QC process
  - Table 2. paired TCR file
  
- Download the paired chain file. 
- If needed also download the TCRdist output file. 

## Other data types

#### Sanger sequencing (single chain)
- Complete steps 1-4

Step 5. Downloading filtered file
- Upload the filled in IMGT_onlyQC.csv to QC -> IMGT (Sanger sequencing) -> tab 2. Paired chain file -> Completed QC file (.csv) 
- Go to the Tab 4. Single chain file, which does not require specifyign if the data was alpha-beta or gamma-delta
- Download single chain file

#### Non-Sanger sequencing data
- Data that has undergone alignment with other processes (e.g. ImmunoSEQ, MiXCR etc.), will require conversion to the TCR_Explore format. 
- Go to QC -> Convert to TCR_Explore file format
- There is two input types ImmunoSEQ and Other. 
- The user can upload either a .tsv, .csv or .txt file. Headers must be in row 1. 
- Select the column with the counts, variable, Diversity, Junction and amino acid column
- Remove all unnecessary columns
- This process re-orders the datasheet so the countColumn is in Column A of the .csv file, as well as adding in TRV, TRJ, TRVJ and TRVJ_CDR3 to the end of the file to aid producing the graphs
- A video explain this process is in available. 

- For the ImmunoSEQ processed data, TCR_Explore will remove rows with missing information (e.g. NA in both V and J genes)
***NOTE:*** *For other sources, the user will have to manually remove non-functional sequence. Contact* __kerry.mullan@monash.edu__ *if you have a specific filtering requirement.*

