# Imports --------------------------------------------------------------------------------------------------------------

from scipy import io
import os
import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.stats import hypergeom
import random
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob

# Constants ------------------------------------------------------------------------------------------------------------

PVALUE = 0.05
DataPathGTEx = os.getcwd() + '/Data/GTEx/'
DataPathCaseCtrl_SMA = os.getcwd() + '/Data/CaseCtrl_SMA/'
DataPathCaseCtrl_CF = os.getcwd() + '/Data/CaseCtrl_CF/'
OutputPath = os.getcwd() + '/Outputs/'
IDT = 56317

# Classes ------------------------------------------------------------------------------------------------------------

class GTEx:
    def __init__(self, LocGenes=None, Tissues=None):
        ''' Constructor for GTEx class
         reate data structure of GTEx RPKM values (analysis V6p data)  '''

        print('Generating GTEx Data structure')

        fns = (glob(DataPathGTEx + '/*.csv'))

        if Tissues is None:
            expf = [i for i in fns if 'GN.csv' not in i]
        else:
            expf = flatten([[i for i in fns if 'GN.csv' not in i and j in i] for j in Tissues])

        GE = [];
        tis = []
        for ef in expf:
            print('loading: ' + ef.split('/')[-1].replace('.csv', ''))
            x = pd.read_table(ef, delimiter=',')
            GE.append(x)
            tis.append([str(ef.split('/')[-1].replace('.csv', '')) for i in range(x.shape[1])])

        self.tissue_spec = flatten(tis)
        if LocGenes is None:
            self.gene = pd.read_csv(DataPathGTEx + 'GN.csv', header=None)[0]
            self.exp = pd.concat(GE, axis=1, sort=False)
        else:
            gene = pd.read_csv(DataPathGTEx + 'GN.csv', header=None)[0]
            ids = [i for i, e in enumerate(gene) if e in LocGenes]
            self.gene = [gene[i] for i in ids]
            exp = pd.concat(GE, axis=1, sort=False)
            self.exp = exp.loc[ids, :]

    def GetStep1Pval(self, GCD, Tissues=None, Qexp=0.1):
        ''' Run GENDULF step 1 - get PM that whose low expression is associated
            with low GCD expression in relevant healthy tissue; Qexp is the
            quantile threshold for low expression; this threshold can be
            vaied in the function SensitivityAnalysis'''

        print('Running GENDULF step1')


        if Tissues == None:
            self.T1 = self.exp.loc[0:IDT, :]
        else:
            TI = [i for i, e in enumerate(self.tissue_spec) if e in Tissues]
            self.T1 = self.exp.loc[0:IDT, TI]

        self.GCDi = [i for i, e in enumerate(self.gene) if e == GCD]

        quants = self.T1.quantile(q=Qexp, axis=1)

        self.T2 = self.T1.lt(quants, axis=0)

        mascG = self.T2.loc[self.GCDi[0]]
        pv = []
        for i in range(0, IDT):
            masci = self.T2.loc[i]
            pvi = hypergeom.sf(sum(mascG & masci), len(mascG), sum(mascG), sum(masci)) if sum(masci) > 0 else 1
            pv.append(float(pvi))

        self.pv1 = pv

    def GetStep3Pval(self, gn, Tissues=None):
        ''' Run GENDULF step 3 - SMN specific step: Get step 3 P-value -
            PM whose low expression is associated with higher FL_SMN2/d7_SMN2 '''

        print('Running GENDULF step3')

        if Tissues == None:
            T1 = self.exp.loc[:, :]
        else:
            TI = [i for i, e in enumerate(self.tissue_spec) if e in Tissues]
            T1 = self.exp.loc[:, TI]

        RT = (T1.loc[56320]) / (T1.loc[56318])  ##Ratio of FL_SMN2/d7_SMN2

        pRATIO = []
        for i in range(len(gn)):
            x = self.gene.index[self.gene == gn[i]]
            vi = T1.loc[x[0]]

            (s, p) = stats.mannwhitneyu(vi[RT > 1], vi[RT < 1], alternative='less')
            try:
                (s, p) = stats.mannwhitneyu(vi[RT > 1], vi[RT < 1], alternative='less')
            except:
                p = 1;
            pRATIO.append(p)

        self.pRATIO = pRATIO

    def sctterD(self, PMn, nameF):
        ''' plot GENDULF step 1 scatters - SMN1 vs. PM expression in healthy tissue '''

        from matplotlib.axes._axes import _log as matplotlib_axes_logger
        matplotlib_axes_logger.setLevel('ERROR')

        PMi = [i for i, e in enumerate(self.gene) if e in PMn]

        fig1, ax1 = plt.subplots(len(PMi), 1, sharex=False, sharey=False)
        for i in range(len(PMi)):
            g2 = self.T1.loc[PMi[i]]
            g1 = self.T1.loc[self.GCDi]

            colors = np.array([0.5, 0.5, 0.5])
            ax1[i].scatter(g1, g2, c=colors, alpha=0.5, edgecolors=colors, label='P = ' + str(self.pv1[PMi[i]]))
            ax1[i].set_xlabel(self.gene[self.GCDi[0]])
            ax1[i].set_ylabel(self.gene[PMi[i]])
            ax1[i].legend(loc=4)
        fig1.savefig(nameF + '.png')


class CaseCTRL:
    def __init__(self, filepath, filen, gnfile):
        ''' Constructor for Case-control data class '''
        ## Create data structure of SMA case-control studies FPKM values
        print('Generating Case-Control Data structure:' + filen)

        x = pd.read_table(filepath + filen + '.csv', delimiter=',')
        gene = pd.read_csv(filepath + gnfile + '.csv', header=None)
        smp = pd.read_csv(filepath + filen + '_SMP.csv', header=None)

        self.exp = x
        self.gene = gene[0]
        self.smp = smp[0]

    def GetCaseControlSMPS(self, diseasen):
        ''' This function returns the indices of case and those of controls in the case-control data '''

        sms = list(filter(lambda x: diseasen in x, self.smp))

        posi = [i for i, e in enumerate(self.smp) if e in sms]
        # print(posi)
        ns = list(filter(lambda x: diseasen not in x and 'CNTL' in x, self.smp))

        negi = [i for i, e in enumerate(self.smp) if e in ns]
        # print(negi)
        self.categ = [diseasen if i in posi else 'CTRL' for i in range(len(self.smp))]
        return posi, negi

    def GetStep2Pval(self, gn, Data, diseasen, DegF=0, rep=10000):
        ''' This function calculates the permutation P-value, examining H0 that the GCD-DPM relation is similar in healthy and disease,
        vs. the one-sided H1 that these are low-low only in healthy tissues '''

        print('Running GENDULF step2')

        gi = [i for i, e in enumerate(self.gene) if e in gn]
        gn2 = [self.gene[i] for i in gi]

        (self.posi, self.negi) = self.GetCaseControlSMPS(diseasen)
        permP = []
        for i in range(len(gi)):
            x = Data.gene.index[Data.gene == gn2[i]]

            GCDe = list(Data.T2.loc[Data.GCDi[0]])
            PM =list(Data.T1.loc[x[0]])
            PMp = [PM[v] for v, z in enumerate(GCDe) if z==True]
            PMn = [PM[v] for v, z in enumerate(GCDe) if z==False]

            GCDSM = self.exp.loc[gi[i]]

            ## True case/ctrl PairScore
            PSC = sum(flatten([[p < n for p in GCDSM[self.posi]] for n in GCDSM[self.negi]]))

            # print(PSC)

            PSCd = []
            for j in range(rep):
                RSM = random.sample(list(PMp), len(self.posi))
                RPM = random.sample(list(PMn), len(self.negi))
                PSCd.append(sum(flatten([[p < n for p in RSM] for n in RPM])))

            permP.append(len([i for i in PSCd if i + DegF <= PSC]) / rep)

        self.gn = gn2
        self.permP = permP

    def boxpltD(self, PMn, nameF):
        ''' plot GENDULF step 2 boxplot of PM in case-control '''
        PMi = [i for i, e in enumerate(self.gene) if e in PMn]

        plt.clf()
        D1 = {'Group': self.categ}
        D2 = {self.gene[PMi[i]]: self.exp.loc[PMi[i]].T for i in range(len(PMi))}
        D1.update(D2)
        df = pd.DataFrame(D1)
        dd = pd.melt(df, id_vars=['Group'], value_vars=[self.gene[PMi[i]] for i in range(len(PMn))], var_name='gene')
        sns.boxplot(x='Group', y='value', data=dd, hue='gene')
        plt.savefig(nameF + '.png')


# Functions  general   -----------------------------------------------------------------------------------------------

flatten = lambda l: [item for sublist in l for item in sublist]

def GetGenesInLoci(geneStart, geneEnd):
    ''' returns all genes in loci from gene geneStart to gene  geneEnd'''
    f = open(DataPathGTEx + 'GN.csv')
    l = [i.replace('\n', '').replace('\ufeff', '') for i in f.readlines()]

    assert geneStart in l and geneEnd in l, "One of the genes not in GTEx"

    ind1 = [i for i, j in enumerate(l) if j == geneStart]
    ind2 = [i for i, j in enumerate(l) if j == geneEnd]

    assert ind1[0] <= ind2[0], "geneStart is located after geneEnd"

    return [l[i] for i in range(ind1[0], ind2[0])]


def makeBarPlot(labels, M1, M2, pv):
    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width / 2, M1, width, label='overall percentage')
    rects2 = ax.bar(x + width / 2, M2, width, label='modifiers percentage')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Scores')
    ax.set_xlabel('Threshold')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()

    def autolabel(rects):
        """Attach a text label above each bar in *rects*, displaying its height."""
        cnt = 0
        for rect in rects:
            height = rect.get_height()
            ax.annotate('{}'.format(pv[cnt]),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom')
            cnt += 1

    autolabel(rects2)

    fig.tight_layout()

    plt.show()
    return plt


#
def SensitivityAnalysis(tis, gene, mods, SavePlot=True):
    """sensitivity analysis for GENDULF step1
    measuring overlap of predicted moodifiers with varying thresholds and plots the results"""
    local_increment = 0.05
    local_num_steps = 7
    thrs = [local_increment * i for i in range(1, local_num_steps)]
    DataL = GTEx(Tissues=[tis])

    preds = []
    for thr in thrs:
        DataL.GetStep1Pval(gene, Qexp=thr)
        preds.append([DataL.gene[i] for i, e in enumerate(DataL.pv1) if e < PVALUE / len(DataL.gene)])

    overlaps = []
    for pr in preds:
        overlaps.append(len(list(set(pr) & set(preds[1]))) / len(preds[1]))

    mod_overlaps = [];
    pmod = []
    for pr in preds:
        mod_overlaps.append(len(list(set(pr) & set(mods))) / len(mods))
        pmod.append(hypergeom.sf(len(list(set(pr) & set(mods))), len(DataL.gene), len(list(set(pr))), len(set(mods))))

    plt = makeBarPlot([str(i)[0:4] for i in thrs], [round(i, 3) for i in overlaps], [round(i, 3) for i in mod_overlaps],
                      [round(i, 20) for i in pmod])

    if SavePlot:
        plt.savefig(OutputPath + 'SensitivityRes' + '.png')

    return overlaps, mod_overlaps, pmod

def PowerAnalysis(GCD, tissues, PM, Pdata = None, SmpNum = 5, rep=10000, iterCnt = 10, PrintRes=False):
    ''' Power analysis to estimate of SmpNum positive and negative samples
    are sufficient for GENDULF step 2, given specific tissues, GCD and PM from step1. '''

    if Pdata is None:
        Pdata = GTEx(Tissues=tissues)
        ##Get the PMs from step1
        Pdata.GetStep1Pval(GCD)

    mascG = list(Pdata.T2.loc[Pdata.GCDi[0]])

    PMi = [i for i, j in enumerate(Pdata.gene) if j == PM]
    expP = list(Pdata.T1.loc[PMi[0]])

    posD = [expP[i] for i, j in enumerate(mascG) if j==True]
    negD = [expP[i] for i, j in enumerate(mascG) if j == False]
    PosPerm = []

    for k in range(iterCnt):
        ## Simulated true modifier (different distributions)
        posS = random.sample(posD, SmpNum)
        negS = random.sample(posD, SmpNum)

        PSC = sum(flatten([[p < n for p in posS] for n in negS]))
        PSCd = []
        for r in range(rep):
            RSM = random.sample(posD, SmpNum)
            RPM = random.sample(negD, SmpNum)
            PSCd.append(sum(flatten([[p < n for p in RSM] for n in RPM])))

        PosPerm.append(len([i for i in PSCd if i <= PSC]) / rep)

    ConfPos = len([i for i in PosPerm if i < PVALUE])/iterCnt

    if PrintRes:
        print("Confidence identifying DPM:: %s \n" % (ConfPos))

    return ConfPos,PosPerm

def MinSampleToCollect(GCD, tissues, ConfTHR = 0.8, MaxIter = 100, PrintProgress=True):
    '''Apply power analysis to PM with a tissues and a GCD, and return the minimal number of samples to collect'''
    DataP = GTEx(Tissues=tissues)
    DataP.GetStep1Pval(GCD)
    g0 = [DataP.gene[i] for i, e in enumerate(DataP.pv1) if e < PVALUE / len(DataP.gene)]
    smpC = 2

    while smpC <= MaxIter:
        MeanConf = []
        if PrintProgress:
            print("Started testing power for %d cases and controls each" % smpC)
        for j in g0:
            ConfPos, PosPerm = PowerAnalysis(GCD, tissues, j, SmpNum=smpC, Pdata=DataP)
            MeanConf.append(ConfPos)

        meanPower = np.mean(MeanConf)    
        if PrintProgress:
            print("Power for %d cases and controls each is %.3f" % (smpC, meanPower))
        if (meanPower>ConfTHR):
            return smpC

        smpC+=1

    ## Maximum iterations reached, no number of sampels
    # below MaxIter that would allow inference with confidence
    return None

# Main Functions -----------------------------------------------------------------------------------------------------


def main_SMA():
    ''' run the code for SMA from start to finish '''

    DataM = GTEx(Tissues=['Muscle - Skeletal'])
    DataM.GetStep1Pval('SMN1')
    DataSP = GTEx(Tissues=['Brain - Spinal cord (cervical c-1)'])
    DataSP.GetStep1Pval('SMN1')

    x0 = [i for i, e in enumerate(DataM.pv1) if e < PVALUE / len(DataM.gene)]  # Bonferroni corrected P<0.05
    x1 = [i for i, e in enumerate(DataSP.pv1) if e < PVALUE / len(x0)]  # Bonferroni corrected for what's left

    v = set(x0).intersection(x1)
    y = [DataM.gene[i] for i in v]

    with open(OutputPath + 'GENDULF_res_SMA_step1.txt', 'w') as f:
        for item in v:
            f.write("%s Muscle P = %s Spinal cord  P = %s \n" % (
            DataM.gene[item], str(DataM.pv1[item]), str(DataSP.pv1[item])))

    SMA_SP = CaseCTRL(DataPathCaseCtrl_SMA, 'SP', 'GN2')
    SMA_SP.GetStep2Pval(y, DataSP, 'SMA')
    SMA_M = CaseCTRL(DataPathCaseCtrl_SMA, 'MS', 'GN2')
    SMA_M.GetStep2Pval(y, DataM, 'SMA')

    xx0 = [SMA_M.gn[i] for i, e in enumerate(SMA_M.permP) if e < PVALUE]
    xx1 = [SMA_SP.gn[i] for i, e in enumerate(SMA_SP.permP) if e < PVALUE]
    v2 = list(set(xx0).intersection(xx1))

    DataM.GetStep3Pval(v2)
    DataSP.GetStep3Pval(v2)  ## This step doesnt pass for any gene (explained in the manuscript and shown in supp.)

    with open(OutputPath + 'GENDULF_res_SMA_step23.txt', 'w') as f:
        for item in v2:
            f.write("%s STEP2 Muscle P = %s Spinal cord P = %s STEP3 Muscle P = %s Spinal cord P = %s \n" % (
            item, str(SMA_M.permP[SMA_M.gn.index(item)]), str(SMA_SP.permP[SMA_SP.gn.index(item)]),
            str(DataM.pRATIO[v2.index(item)]), str(DataSP.pRATIO[v2.index(item)])))

    mod = ['U2AF1', 'SF1', 'SRSF4']

    ####Plot U2AF1,SF1,SRSF4
    DataM.sctterD(mod, OutputPath + 'GENDULF_SMA_step1M')
    DataSP.sctterD(mod, OutputPath + 'GENDULF_SMA_step1SP')
    SMA_M.boxpltD(mod, OutputPath + 'GENDULF_SMA_step2M')
    SMA_SP.boxpltD(mod, OutputPath + 'GENDULF_SMA_step2SP')


def main_CF():
    ''' run the code for CF from start to finish '''

    ##Run Lung tissue
    DataL = GTEx(Tissues=['Lung'])
    DataL.GetStep1Pval('CFTR')

    x0 = [i for i, e in enumerate(DataL.pv1) if e < PVALUE / len(DataL.gene)]  # Bonferroni corrected P<0.05

    y = [DataL.gene[i] for i in x0]

    with open(OutputPath + 'GENDULF_res_CF_Lung_step1.txt', 'w') as f:
        for item in x0:
            f.write("%s Lung P = %s \n" % (DataL.gene[item], str(DataL.pv1[item])))

    CF_LUNG = CaseCTRL(DataPathCaseCtrl_CF, 'LUNG', 'GN1')
    CF_LUNG.GetStep2Pval(y, DataL, 'CF')
    xx0 = [CF_LUNG.gn[i] for i, e in enumerate(CF_LUNG.permP) if e < PVALUE]

    with open(OutputPath + 'GENDULF_res_CF_LUNG_step2.txt', 'w') as f:
        for item in xx0:
            f.write("%s STEP2 P = %s \n" % (
                item, str(CF_LUNG.permP[CF_LUNG.gn.index(item)])))

    ##Run colon tissue
    DataC = GTEx(Tissues=['Colon - Transverse', 'Colon - Sigmoid'])
    DataC.GetStep1Pval('CFTR')

    x0 = [i for i, e in enumerate(DataC.pv1) if e < PVALUE / len(DataC.gene)]  # Bonferroni corrected P<0.05

    y = [DataC.gene[i] for i in x0]

    with open(OutputPath + 'GENDULF_res_CF_Colon_step1.txt', 'w') as f:
        for item in x0:
            f.write("%s Colon P = %s \n" % (DataC.gene[item], str(DataC.pv1[item])))

    CF_COLON = CaseCTRL(DataPathCaseCtrl_CF, 'COLON', 'GN2')
    CF_COLON.GetStep2Pval(y, DataC, 'CF')
    xx0 = [CF_COLON.gn[i] for i, e in enumerate(CF_COLON.permP) if e < PVALUE]

    with open(OutputPath + 'GENDULF_res_CF_COLON_step2.txt', 'w') as f:
        for item in xx0:
            f.write("%s STEP2 P = %s \n" % (
                item, str(CF_COLON.permP[CF_COLON.gn.index(item)])))

    ##sensitivity analysis:
    MOD_CF_LUNG = ['TGFB1', 'MBL2', 'TNFA', 'IL10', 'GSTM1', 'GSTP1', 'ADRB2', 'EDNRA', 'IFRD1', 'IL8', 'ACE',
                   'SLC26A9', 'NOS3', 'NOS1', 'AAT', 'HLA2', 'B2AR', 'CLC2', 'MUC4', 'MUC20', 'SLC9A3', 'AGTR2',
                   'SLC6A14', 'EHF', 'APIP', 'SFTPA1', 'SFTPA2', 'SFTPB', 'SFTPC', 'SFTPD', 'KRT8']

    overlaps, mod_overlaps,pmod = SensitivityAnalysis('Lung', 'CFTR', MOD_CF_LUNG)


