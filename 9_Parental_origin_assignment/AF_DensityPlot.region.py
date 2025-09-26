import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests

parser = argparse.ArgumentParser(description="Large-scale mutation density plot.")
parser.add_argument("--input", "-i", help="phased_input", required=True)
parser.add_argument("--output", "-o", required=True)
parser.add_argument("--region", "-r", required=True)

args = parser.parse_args()

AF_file = args.input
output_path = args.output
region = args.region

if not os.path.exists(output_path):  
    os.makedirs(output_path)

def parse_region(region_str):
    try:
        chrom, coords = region_str.split(":")
        start, end = map(int, coords.split("-"))
        if start >= end:
            raise ValueError("Start position must be less than end position")
        return chrom, start, end
    except Exception as e:
        raise ValueError(f"Invalid region format. Expected 'chr:start-end'. Error: {e}")

def deplot_AF_region(i,sample,data):
    F_AF = data
    M_AF = 1-F_AF
    sns.kdeplot(x=F_AF,ax=ax_sample,color='blue',label='F',warn_singular=False)
    sns.kdeplot(x=M_AF,ax=ax_sample,color='red',label='M',warn_singular=False)

    _, p_value = wilcoxon(F_AF, M_AF)
    ax_sample.text(0.5,0.1,'Wilcoxon, p='+format(p_value,'.4f'),ha='center', va='top',transform=ax_sample.transAxes)
    ax_sample.set_title(sample)
    ax_sample.set_xlabel("AF")

chrom, start, end = parse_region(args.region)

sample_number = list(range(1,8))
sample_list = ['B','P1','P2','P3','P4','P5','E']

fig = plt.figure(figsize=(11,5))
for i in sample_number:
    sample = sample_list[i-1]
    AF_data = pd.read_csv(input,sep='\t',header=None,low_memory=False,na_values="")
    AF_data = AF_data.dropna(axis=1, how='all')
    AF_data = AF_data[(AF_data.iloc[:,1] >= int(start)) & (AF_data.iloc[:,1] <= int(end))]
    AF_data.columns = ['chr','pos']+sample_list

    ax_sample = fig.add_subplot(2, 4, i)
    deplot_AF_region(i,sample,AF_data[sample])
plt.suptitle(region, x=0.5, y=0.97, fontsize=15)
lines, labels = fig.axes[-1].get_legend_handles_labels()
fig.legend(lines, labels, bbox_to_anchor=(0.82,0.22), prop={'size': 10})
plt.subplots_adjust(wspace=0.5,hspace=0.5)

fig.savefig(f'{output_path}/{region}_density.jpg')
plt.close()
