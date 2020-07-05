# GENDULF_CODES


This file describes a python script that implements the method GENetic
moDULators identiFication (GENDULF) as applied to data for spinal
muscular atrophy (SMA) and Cystic Fibrosis (CF) [1].
 
![Figure1](https://user-images.githubusercontent.com/18428559/68185201-fd242280-ff6e-11e9-9866-17c97408d431.png)


The purpose of GENDULF is to find possible modifier genes for
autosomal recessive diseases in which the primary causative mutations
are loss-of-function mutations that diminish or abrogate the
expression of the mutated genes. In the initial study, we focused on
two of the most common recessive diseases cystic fibrosis (CF), which
is caused by mutations in CFTR, and spinal muscular atrophy (SMA)
which is caused by mutations in SMN1. As implemented, GENDULF searches
for modifier genes for which decreased expression of the modifier is
associated with milder disease (decreased severity).


# Code to run (for all results and some basic plots):
GENDULF_CODE.py

The code should be able to run in any python3 environment. The code
has been tested on python3.6 and python3.7

To run the SMA code from start to finish enter the command:
python3 -c "import GENDULF_CODE; GENDULF_CODE.main_SMA()”



Alternatively, start by entering the python interpreter with
python3

and then in response to the interpreter prompt (shown here as >>>) enter
>>> import GENDULF_CODE
>>> GENDULF_CODE.main_SMA()

The code relies on the availability of seven widely-used python
packages: scipy, pandas, numpy, random, os ,glob, matplotlib and seaboarn

GENDULF_CODE.py code has 2 classes

1. GTEx - builds the GTEx data structure for a selected tissue (tissue
name given as input) and has functions for the GENDULF steps that use
GTEx data - step 1 and step 3 (step 3 is the SMA-specific step)

For GTEx data, the valid gene symbols with respect to this version of GTEx are in
Data/GTEx/GN.csv
and the valid tissues are the prefixes of file names in
Data/GTEx

2. CaseCTRL - builds the case-control gene expression data structure, and
has a function for GENDULF step 2 which uses the case-control data.

The inputs data are in the Data/ directory

# Input files
1. GTEx/ - .csv of each GTEx tissue and gene names
2. CaseCtrl_SMA/ - case control data for SMA
3. CaseCtrl_CF - case control data for CF


# The main_SMA script does seven steps:

1. Generates GTEx data structures for muscle and spinal cord

2. Runs GENDULF step1 for muscle and then for spinal cord, and
identifies potential modifiers with corrected P-value (the P-value
correction for spinal cord is applied to the filtered targets passing
muscle tissue)
None of the p-values for any specific gene depend on  how many or which genes are tested.
The order in which the genes in a set are listed does not matter.

3. Generates SMA case-control data structures for muscle and for spinal cord

4. Runs GENDULF step 2 for muscle and for spinal cord for the PMs
(potential modifiers) passing step 1 for both tissues, and identifies
PMs that pass step 2 in both tissues

5. Runs GENDULF step 3 for muscle, for the PM passing GENDULF step 2 

6. Runs GENDULF step 3 for spinal cord, for the PM passing GENDULF
step 2 (not used to identify targets, not enough samples for this
step, explained in manuscript)

7. Saves the results of the steps, and plots for U2AF1, SF1 and SRSF4

By default, all output files are saved in the subdirectory Outputs.
In the code distribution from GitHub, this output directory does not exist a priori. Therefore, the user should
create the directory using the UNIX command:
mkdir Outputs
or more generally

mkdir <Name of output directory>
before running GENDULF.

In the GENDULF code, the default name of the output directory can be changed by changing the line of code
OutputPath = os.getcwd()+'/Outputs/'
As presently implemented one should do only a single run of GENDULF at a time, so as not to overwrite the files in
the Outputs subdirectory.

# Output files for SMA (in Outputs/):
1. GENDULF_res_SMA_step1.txt gene name and GENDULF step1 P-values for PMs identified in GENDULF step1 for both muscle and spinal cord
2. GENDULF_res_SMA_step23.txt gene name and GENDULF steps 2 and 3 P-values for PMs identified in GENDULF step2 for both muscle and spinal cord
3. GENDULF_SMA_step1M.png/GENDULF_step1SP.png setter plots for GENDULF step 1 for U2AF1, SF1 and SRSF4, for muscle and spinal cord, respectively 
4. GENDULF_SMA_step2M.png/GENDULF_step2SP.png box plots for GENDULF step 1 for U2AF1, SF1 and SRSF4, for muscle and spinal cord, respectively 

In GENDULF_SMA_res_step1.txt, the P-value threshold to pass the muscle
tissue test is 0.05/(56330 - #transcripts in GTEx) = 8.8e-7. For the
spinal cord tissue test, the P-value threshold is
(0.05/#remaining_candidates) = 4.52e-5

In GENDULF_res_step23.txt the P-value threshold to pass step 2 is
0.05/#remaining_candidates = 0.00014 for muscle and 0.00011 for spinal
cord.

In GENDULF_res_step23.txt P = 0.0 really means P < 1/1000; if one
wishes to increase the stringency, one can increase the number of
permutations simply by changing the line of code that reads

    rep = 1000 

In all outputs, the order of the rows from top to bottom is arbitrary.

# How to run GENDULF for a list of genes such as candidates derived from a prior from GWAS in one or more loci:

When running GENDULF step 1, define the parameter LocGenes to be the list of genes.
For example: LocGenes = [‘gene1’,’gene2’,’gene3’]. To get all genes in locus between 2 given bounding 
genes on the same chromosome, use the GetGenesInLoci function. In general, the genes in the list do not need to be adjacnt or even on the
same chromosome. When selecting specific genes by gene symbol, make sure all selected gene symbols are in the file
Data/GTEx/GN.csv

# How to run for other diseases:
The alternative main function main_CF that runs GENDULF for cystic fibrosis exemplifies how to run for different dieases.
To use that run:
>>> GENDULF_CODE.main_CF()

Unlike the SMA example, the CF example does not have a step 3 and thus, may be more typical of usage for other diseases.

# Sensitivity analysis:
One internal parameter defines the range of expression level ranks considered as low. The baseline setting for
this is 0.1 (i.e, 10th percentile), which works well for SMA and CF, but may be adjusted for other diseases.To change this, change the value of Qexp in GetStep1Pval.
The user may use the function SensitivityAnalysis to evaluate how the results change if the threshold 0.1
To apply a sensitivity analysis, use the SensativityAnalysis function.
The SensitivityAnalysis performs: (1) measure overlap in predicted modifier with these predicted with 0.1 threshold, and (2) measure overlap and P-value with actual modifiers reported in the literature (thats the list of genes).
It is applied for CF lung tissues, within the main_CF, and saves a figure with the results to Outputs/SensitivityRes.png
The X-axis in SensitivityRes.png are threshold used, and the Y-axis are the (a) Blue - overlap of predicted modifiers with those used with 0.1 threshold (0 is none and 1 is all), and (b) Orange - the percentage of true modifiers (given as mods variable) that are predicted with each threshold along the X-axis. 
Currently the alternative values tested are 0.05, 0.1, 0.15, 0.2, 0.25, 0.3
and this is controlled by the list assignment
local_increment = 0.05
local_num_steps = 7
thrs = [local_increment * i for i in range(1, local_num_steps)]
which the user can change. 

It is deliberate that the loop for sensitivity analysis reruns the
analysis including the baseline setting (default 0.1) as one of the
values, so that in the plot one can compare the results for the
baseline setting to other settings; the plot uses the same colors for
every value on the x-axis and does not distinguish the baseline
setting from the alternatives by color.  The intent is that the user
will run once with the baseline value, then modify the array currently
called MOD_CF_LUNG to list the modifier genes found at the baseline
setting and call SensitivityAnalysis with the array of genes as the
last argument, and rerun including a call to SensitivityAnalysis.  In
the code as distributed, the value of MOD_CF_LUNG is pre-set.

# Power Analysis:
To estimate the number of case and control samples required to accept
DPM with confidence for a specific tissue and GCD, run the MinSampleToCollect
Function. 
For example, to run it for Lung tissue with CFTR as the GCD, use the following:
> smpCollect = GENDULF_CODE.MinSampleToCollect('CFTR', ['Lung'])

The required first argument is the disease-causing gene in single quotes.
The required second argument is a comma-separated non-empty list of GTEx tissues.
The optional parameter ConfTHR is a floating point threshold for the  mean power inferred for all PMs (default is 0.8). This value should be between 0 and 1.
The optional parameter MaxIter is an integer and is the maximal number of iterations (or maximal number of case and control samples one is willing to collect).
The optional parameter PrintProgress, by default set to True, prints an update each time the number of samples is incremented and each time the
power has been estimated for the newly incremented number of samples
 
The procedure GENDULF_CODE.MinSampleToCollect  estimates power iteratively for the number of samples being 2, 3, 4, ..., MaxIter until either the number of samples
equals MaxIter or the power exceeds the specified threshold of ConfTHR. If the power never exceeds ConfTHR, then None is returned. If the power does exceed
ConfTHR, then the minimum number of samples achieving the sufficient power is returned as an integer value. To see the actual power estimates, one must
keep PrintProgress set to True and save the contents stdout, where the updates are printed.

Importantly, GENDULF_CODE.MinSampleToCollect internally runs step 1, which uses GTEx, but not step 2 because the intent of the power calculation is to estimate the number 
of samples needed at step 2. One does not need to run step 1 separately for GENDULF_CODE.MinSampleToCollect with default parameter settings. In practice,
step 2 involves collecting both cases and controls, but varying both numbers independently is time consuming. Therefore, we made the time-saving
assumption that the two numbers are equal. An alternative time-saving assumption would be that one number is fixed and the other number varies. To 
implement that assumption, one would modify the code in the helper pocedure PowerAnalysis as follows. At the two pairs of lines of code
        posS = random.sample(posD, SmpNum)
        negS = random.sample(posD, SmpNum)
and
            RSM = random.sample(posD, SmpNum)
            RPM = random.sample(negD, SmpNum)
one would change one occurrence of SmpNum to a constant and keep the other as a variable; the substring pos represents cases and the substring neg
represents controls. If for example, the number of controls was fixed at 50, then one should add the line
NUM_CONTROLS = 50
and change the above four lines to:
        posS = random.sample(posD, SmpNum)
        negS = random.sample(posD, NUM_CONTROLS)
and
            RSM = random.sample(posD, SmpNum)
            RPM = random.sample(negD, NUM_CONTROLS)



Reference:
[1] Auslander N, Ramos DM, Zelaya I, Karathia H, Crawford
TO, Schaffer AA, Sumner CJ, Ruppin E.  The GENDULF algorithm: mining
transcriptomics to uncover modifier genes for monogenic diseases
caused by loss-off-function mutations, submitted.


Contact:
Noam Auslander noam.auslander@nih.gov
