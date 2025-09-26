import argparse
import pandas as pd
from scipy.stats import wilcoxon

parser = argparse.ArgumentParser(description="LOH boundary refinement.")
parser.add_argument("--input", "-i", help="phased_input", required=True)
parser.add_argument("--region", "-r", required=True)

args = parser.parse_args()

AF_file = args.input
AF_data = pd.read_csv(AF_file,sep='\t',header=None,low_memory=False)
AF_data = AF_data.dropna(axis=1, how='all')
region = args.region

def redefine_LOH_region(AF_data,origin_region):
    origin_start = origin_region.split(":")[1].split("-")[0]
    origin_end = origin_region.split(":")[1].split("-")[1]
    start_all = []
    end_all = []
    if AF_data.shape[0] > 500:
        return origin_start,origin_end
    else:
        for i in range(2,AF_data.shape[1]):
            d_AF = pd.DataFrame(list(AF_data.iloc[:,1]),columns=['pos'],dtype=int)
            d_AF['dAF'] = AF_data.iloc[:,i]*2-1
            d_AF['win'] = d_AF.dAF.rolling(window=5,center=True).mean()
            if d_AF.dAF.mean() > 0:
                d_AF['is_LOH'] = (d_AF.win > d_AF.dAF.mean()/2) | (d_AF.win > 3*d_AF.dAF.var()) | (d_AF.win > 0.1)
            else:
                d_AF['is_LOH'] = (d_AF.win < -d_AF.dAF.mean()/2) | (d_AF.win < -3*d_AF.dAF.var()) | (d_AF.win < -0.1)
            d_AF['crossing'] = (d_AF.is_LOH != d_AF.is_LOH.shift()).cumsum()
            d_AF['counts'] = d_AF.groupby(['is_LOH', 'crossing']).cumcount(ascending=True) + 1
            d_AF.loc[d_AF.is_LOH == False, 'counts'] = 0
            start = d_AF.shift(d_AF.counts.max()-1).loc[d_AF.counts == d_AF.counts.max()]['pos'].values[0]
            end = d_AF.loc[d_AF.counts == d_AF.counts.max()]['pos'].values[0]
            F_AF = AF_data[(d_AF['pos'] >= int(start)) & (d_AF['pos'] <= int(end))].iloc[:,i]
            M_AF = 1-F_AF
            origin_F_AF = AF_data[(d_AF['pos'] >= int(origin_start)) & (d_AF['pos'] <= int(origin_end))].iloc[:,i]
            origin_M_AF = 1-origin_F_AF
            _, p_value = wilcoxon(F_AF,M_AF)
            _, origin_p_value = wilcoxon(origin_F_AF,origin_M_AF)
            if p_value < 0.01 and origin_p_value < 0.01:
                start_all.append(int(start))
                end_all.append(int(end))
        if len(start_all)==0 or len(end_all)==0:
            return 0,0
        else:
            re_start = max(start_all)
            re_end = min(end_all)
            return re_start,re_end
    

re_start,re_end = redefine_LOH_region(AF_data,region)
print(re_start,re_end)

