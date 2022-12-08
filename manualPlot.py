from pathlib import Path
from os import listdir, path
from natsort import natsorted

import matplotlib.pyplot as plt
import seaborn as sns
import colorcet as cc
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from scipy.signal import medfilt
#  from sklearn import datasets, linear_model
#  from sklearn.metrics import mean_squared_error, r2_score

def boxPlot(data, title_list, color_list, segment, folderName):
    #  matplotlib.rc('font', family='SimHei')
    f = plt.figure(figsize=(4, 6))
    medianprops = dict( color='#323232')
    bp = plt.boxplot(data.transpose(),
                     vert=True,
                     patch_artist=True,
                     labels=title_list,
                     medianprops=medianprops
    )
    for patch, color in zip(bp['boxes'], color_list):
        patch.set(facecolor = color, alpha =0.8)

    nn = len(data)
    # Add it to the plot
    row_means_ori = data.mean(axis=1)
    row_means = [None] * nn
    ratio_means = [None] * nn

    max2d = max(map(max, data)) # max in 2d, used for text
    for i in range(nn):
        row_means[i] = round(row_means_ori[i], 2)

    segmentCount = 1
    for i in range(len(data)):
        if i >= segment[segmentCount - 1] and i < segment[segmentCount]:
            ratio_means[i] = row_means[i] / row_means[segment[segmentCount] - 1]
            if i == segment[segmentCount] - 1:
                segmentCount = segmentCount + 1
        ratio_means[i] = round(ratio_means[i], 2)


    pos = range(len(data))
    #  for tick in range(nn):
        #  plt.text(pos[tick]+1, 302, ratio_means[tick],\
        #  horizontalalignment='center', size='x-small', color='C0', weight='semibold')
#
        #  plt.text(pos[tick]+1, 315, row_means[tick],\
        #  horizontalalignment='center', size='x-small', color='C2', weight='semibold')
#
    #  plt.suptitle('boxPlot_' + label + '_' + 'P' +  str(i))
    plt.suptitle(folderName)
    #  plt.ylim([0, 300])
    plt.ylabel('Bead number per 100 X 100 $\mu$$m^2$')
    plt.xticks(rotation=45)
    plt.xticks(fontsize = 10)
    plt.yticks(fontsize = 12)
    #  plt.yscale('log')
    plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.48)
    plt.savefig(f"{folderName}_rm.png", dpi=300)
    #  plt.show()
    return f


def adjust_lightness(color, amount=0.5):
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])

def colorGroup(color_list, segment):
    segmentCount = 0
    for i in range(len(color_list)):
        if i >= segment[segmentCount] and i < segment[segmentCount + 1]:
            color_list[i] = adjust_lightness('C'+str(segmentCount), 1)
        if i == segment[segmentCount + 1] - 1:
            segmentCount = segmentCount + 1

    return color_list

def loadData(dataFolder):
    # load txt data and its conc condition
    rawData_ls = []
    basepath = path.abspath('')
    #  basepath = Path(__file__).parent
    #  basepath = basepath / "expData_all"
    basepath = path.join(basepath, dataFolder)
    #  basepath = basepath / dataFolder

    for file in listdir(basepath):
        if file.endswith(".txt"):
            rawData_ls.append(file)

    rawData_ls.sort()

    return rawData_ls, basepath

def selectByMagets_Z(df):
    # select the largest two Magnets_Z postition
    df['freq'] = df.groupby('Magnets_Z')['Magnets_Z'].transform('count')
    selected_Z = df['freq'][-2:]
    df.loc[df['freq'].isin(selected_Z)]
    df.reset_index(drop=True, inplace=True)
    return df

def removeNan(df):
    # ----- remove one section of nan rows and connect time series -----
    #  df.dropna( axis=0, inplace=True) # simply but didn't connect time
    nanRow = df.loc[df.isna().any(axis=1)]
    temp1 = df.iloc[:nanRow.index[0]].copy()
    temp2 = df.iloc[nanRow.index[-1]+1:].copy()
    skipTime = nanRow['Time'].iloc[-1]- nanRow['Time'].iloc[0]
    temp2['Time'] = temp2['Time'] - skipTime
    df = pd.concat([temp1, temp2])
    return df

def loadDF(rawData_ls):
    DF = pd.DataFrame()
    medfiltFlag =True

    for i, file in enumerate(rawData_ls):
        df = pd.read_csv(path.join(basepath, file), header=None, sep='\t')
        new_header = df.iloc[0]
        df = df.iloc[2:]
        df.columns = new_header
        df.reset_index(drop=True, inplace=True)

        #  breakpoint()

        # convert to float32 before substraction calc
        df = df.astype('float')
        df['Time'] = df['Time'] - df['Time'][0]
        df['Extension'] = df['Extension'] - df['Extension'][0] + (500 * i)
        if medfiltFlag:
            df['Extension'] = medfilt(df['Extension'], kernel_size=3)

        # label DNA type
        if '2K' in file:
            df['DNA'] = '2K'
        else:
            df['DNA'] = '7K'

        # apply 1D GaussianFilter
        sigma = df.shape[0] / df['Time'].tail(1).values

        df['Time_gaus'] = gaussian_filter(df['Time'], sigma)
        df['Ext_gaus'] = gaussian_filter(df['Extension'], sigma)
        df['File'] = i
        #  df['File'] = file.split('-')[-1]

        df = df.loc[df['Time'] > 120]


        DF = pd.concat([DF, df])
    return DF

def GaussianFilterPlot(sigma):
    z = gaussian_filter(df['Extension'], sigma=sigma)
    t = gaussian_filter(df['Time'], sigma=sigma)
    sns.lineplot( x=t, y=z)

#  ====== main =======
folderName = 'temp'
#  folderName = 'Figure 2_LATRX_highConc'
#  folderName = 'Figure 1_2K_7K'
#  if False:
if True:
    rawData_ls, basepath = loadData(folderName)
    print(rawData_ls)
    #  dataFrame = np.empty([1, 15])
    DF = loadDF(rawData_ls)
    DF.to_csv(f"{folderName}.csv", index=False)

DF = pd.read_csv(f"{folderName}.csv")

#  breakpoint()

DF.reset_index(drop=True, inplace=True)
df_t = DF[['Time', 'Extension', 'DNA', 'File']].copy()
df_t.loc[:, 'DataType'] = 'raw'

df_t_gaus = DF[['Time_gaus', 'Ext_gaus', 'DNA', 'File']].copy()

df_t_gaus.rename(columns={'Time_gaus': 'Time', 'Ext_gaus': 'Extension'}, inplace=True)
df_t_gaus['DataType'] = 'gaus'
#  df = pd.concat([df_t_gaus, df_t])
df = pd.concat([df_t, df_t_gaus])

df['DataType_File'] = df['File'].astype(str) + df['DataType']
#  df['DataType_File'] = df['DataType'] + df['File'].astype(str)

#  DF['DNA_File'] = DF['DNA'] + DF['File'].astype(str)

#  print('done')
#  breakpoint()
#  sns.set_theme(style='whitegrid')
sns.set_theme(style='darkgrid')
fig, ax = plt.subplots(figsize=(12, 8))
#  sns.lineplot(data=DF, x='Time', y='Extension', hue='DNA',
             #  units='File', estimator=None,
             #  palette='flare')

hue_order = df['DataType_File'].unique()
print(hue_order)
#  natsorted(hue_order)
hue_order.sort()
#  hue_order = hue_order.reshape(-1, 2)[::-1, ::-1].flatten().tolist()
hue_order = hue_order.reshape(-1, 2)[:, ::-1].flatten().tolist()
print(hue_order)

#  breakpoint()
sns.lineplot(data=df, x='Time', y='Extension', hue='DataType_File',
             hue_order=hue_order,
             #  units='DataType_File', estimator=None,
             #  palette= cc.fire,
             palette='Paired',
             ax=ax)
#  sns.lineplot(data=DF, x='Time', y='Extension', hue='DNA')
#  GaussianFilterPlot(86)

#  sns.lineplot(data=df, x='t', y='z')
#  plt.legend('')
#  plt.savefig(f"{folderName}.svg", format='svg', dpi=75)
# change the fontsize
fontsize=20
plt.xlabel('Time (s)', fontsize=fontsize);
ax.tick_params(axis='x', labelsize=fontsize)
plt.ylabel('Bead height (nm)', fontsize=fontsize);
ax.tick_params(axis='y', labelsize=fontsize)
plt.legend([],[], frameon=False)
plt.savefig(f"{folderName}.png", format='png', dpi=150)
plt.show()

breakpoint()


# ==================
#
# the end
#
# =================

ratioCircle = np.loadtxt('2022-08-30_CR3022 spike-in experiment_All.txt')
#  print(ratioCircle)
#  ratioCircle = np.flip(ratioCircle, axis=0)
#  print(ratioCircle)

conc = np.loadtxt('2022-08-30_POCT_Conc.txt')
#  print(conc)
#  conc = np.flip(conc)
#  print(conc)
#  breakpoint()
#  concData = np.array(conc).transpose()
concData = np.array([ conc ])
concData = np.tile(concData.T, (1, 4))
concData = concData.flatten()
concData = [f"CR3022 { str(float(x)) } pM"  for x in concData]
plt.style.use('seaborn')
N = len(concData)
color_list = [None]*N
segment = list(range(0, 25, 4))
#  segment = [0, 6, 7, 10, 13,  16, 24]
color_list = colorGroup(color_list, segment)
#  boxPlot(ratioCircle, concData, color_list, segment, 'individual testing well')
#  plt.show()
#  breakpoint()
#  np.savetxt(f"ratioCircle after reshape(6, 16).txt", ratioCircle)
#  print(ratioCircle)
# -------------- POCT =================
# combine 4wells as one box for each POCT kit
#  ratioCircle = np.loadtxt('data_2022-08-30_Blood dilution.txt')
#  cc = np.loadtxt('Blood dilution.txt')
cc = np.loadtxt('2022-08-30_POCT_Conc.txt')
conc = [np.nan] * len(cc)
for i, c in enumerate(cc):
    if c == 0.8:
        conc[i] = f"CR3022 { str(float(c)) } pM"
        #  conc = [f"CR3022 { str(float(x)) } pM"  for x in conc]
    else:
        conc[i] = f"CR3022 { str(int(c)) } pM"
        #  conc = [f"CR3022 { str(int(x)) } pM"  for x in conc]
#  conc = [f"CR3022 { str(float(x)) } pM"  for x in conc]
#  breakpoint()
ratioCircle = ratioCircle.reshape(-1, 16)
#  np.savetxt(f"ratioCircle after reshape(6, 16).txt", ratioCircle)

N = len(conc)
color_list = [None]*N
segment = [0, N]
#  segment = list(range(N+1))
color_list = colorGroup(color_list, segment)
#  color_list = [None] * N
plt.style.use('seaborn')
boxPlot(ratioCircle, conc, color_list, segment, 'POCT')
plt.show()
#  boxPlot(ratioCircle, title_list, color_list, segment, folder_ls[k])
#  print(conc)
#  print(SMS)
#  import pdb; pdb.set_trace()
breakpoint()

concData = np.tile(np.array([conc]).transpose(), (1, 16))
#  print(f"transpose(1,4){ concData }")
SMS = SMS.flatten().reshape(-1, 1)
#  SMS = SMS.flatten()
#  print(f"flatten and reshape{SMS}")
#  breakpoint()
concData = concData.flatten().reshape(-1, 1)
#  import pdb; pdb.set_trace()

# Create linear regression object
regr = linear_model.LinearRegression()

# Train the model using the training sets
#  regr.fit(diabetes_X_train, diabetes_y_train)
regr.fit(concData, SMS)
#  import pdb; pdb.set_trace()

# Make predictions using the testing set
SMS_y_pred = regr.predict(concData)

dev = SMS - SMS_y_pred
SD = np.mean(dev ** 2)
LOD = 3.3 * SD / regr.coef_
LOD = LOD[0][0]

# the LOD 
print('LOD: \n', LOD)
# The coefficients
print("Coefficients: \n", regr.coef_)
# The mean squared error
print("Mean squared error: %.2f" % mean_squared_error(concData, SMS))
# The coefficient of determination: 1 is perfect prediction
print("Coefficient of determination: %.2f" % r2_score(SMS, SMS_y_pred))

# Plot outputs
plt.style.use('seaborn')

plt.scatter(concData, SMS, color="black")
plt.plot(concData, SMS_y_pred, color="blue", linewidth=3)
plt.xlabel('Spike Concentration (pM)')
plt.text(np.mean(conc), max(SMS), 'LOD: 3.3x$\u03C3$/slop = '+ ' '+ str(round(LOD, 2)),\
          horizontalalignment='center', size='x-large', color='k', weight='semibold')

plt.ylabel('Bead number per 100 X 100 $\mu$$m^2$')
#  plt.xticks(())
#  plt.yticks(())
#  import pdb; pdb.set_trace()

plt.show()

