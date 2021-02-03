'''
Created on Aug 20, 2020

@author: mina
'''
import numpy as np

seq_dic = {}
var_list = ['wt', 'wtAro', 'allY', 'allF', 'aro1', 'f15y',
                         'f2y', 'f3y', 'k2g', '2r', '7r', '-10r',
                         '12d', '8d', '4d', '-4d', '7k12d', 'rk1',
                         'rk2', '7r12d', '-2k-4r5d']

temp_list = [4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 25, 28, 30, 32, 34, 35, 36, 38, 40]
fmt_list = ['o', 's', 'D', '^', 'v', '>', '<', 'p']

for var in var_list:
    seq_dic[var] = {}
    seq_dic[var]['r'] = 10
    seq_dic[var]['k'] = 2
    seq_dic[var]['d'] = 4
    seq_dic[var]['e'] = 0
    seq_dic[var]['y'] = 8
    seq_dic[var]['f'] = 12
    seq_dic[var]['aro'] = 20
    seq_dic[var]['csat'] = {}
    for temp in temp_list:
        seq_dic[var]['csat'][str(temp)] = 0

#Names
seq_dic['wt']['name'] = 'WT+NLS'
seq_dic['wtAro']['name'] = 'WT-NLS'
seq_dic['k2g']['name'] = '-2K'
seq_dic['aro1']['name'] = '-4F-2Y'
seq_dic['allY']['name'] = '-12F+12Y'
seq_dic['allF']['name'] = '+7F-7Y'
seq_dic['f3y']['name'] = '-9F+3Y'
seq_dic['f2y']['name'] = '-8F+4Y'
seq_dic['f15y']['name'] = '-9F+6Y'
seq_dic['2r']['name'] = '+2R'
seq_dic['7r']['name'] = '+7R'
seq_dic['-10r']['name'] = '-10R'
seq_dic['12d']['name'] = '+12D'
seq_dic['8d']['name'] = '+8D' 
seq_dic['4d']['name'] = '+4D'
seq_dic['-4d']['name'] = '-4D'
seq_dic['7k12d']['name'] = '+7K+12D'
seq_dic['rk1']['name'] = '-3R+3K'
seq_dic['rk2']['name'] = '-6R+6K'
seq_dic['7r12d']['name'] = '+7R+12D'
seq_dic['-2k-4r5d']['name'] = '-4R-2K+5D'

#Color Palette 1
seq_dic['wt']['color'] = 'black'
seq_dic['wtAro']['color'] = 'black'

seq_dic['allY']['color'] = (.8, .8, 0)
seq_dic['f15y']['color'] = 'gold'
seq_dic['f2y']['color'] = 'goldenrod'
seq_dic['f3y']['color'] = 'darkorange'
seq_dic['allF']['color'] = 'chocolate'
seq_dic['aro1']['color'] = 'silver'

seq_dic['7r']['color'] = 'navy'
seq_dic['2r']['color'] = 'blue'
seq_dic['-10r']['color'] = 'lightblue'
seq_dic['7r12d']['color'] = 'rebeccapurple'

seq_dic['12d']['color'] = 'darkred'
seq_dic['8d']['color'] = 'red'
seq_dic['4d']['color'] = 'salmon'
seq_dic['-4d']['color'] = 'pink'

seq_dic['rk2']['color'] = 'olivedrab'
seq_dic['rk1']['color'] = 'limegreen'
seq_dic['k2g']['color'] = 'springgreen'
seq_dic['7k12d']['color'] = 'teal'

seq_dic['-2k-4r5d']['color'] = 'turquoise'

#Color Palette Fig 2
seq_dic['wt']['colorFig2'] = 'black'
seq_dic['wtAro']['colorFig2'] = 'black'

seq_dic['allY']['colorFig2'] = 'black'
seq_dic['f15y']['colorFig2'] = 'black'
seq_dic['f2y']['colorFig2'] = 'black'
seq_dic['f3y']['colorFig2'] = 'black'
seq_dic['allF']['colorFig2'] = 'black'
seq_dic['aro1']['colorFig2'] = 'black'

seq_dic['7r']['colorFig2'] = 'royalblue'
seq_dic['2r']['colorFig2'] = 'royalblue'
seq_dic['-10r']['colorFig2'] = 'royalblue'
seq_dic['7r12d']['colorFig2'] = 'blue'

seq_dic['12d']['colorFig2'] = 'red'
seq_dic['8d']['colorFig2'] = 'red'
seq_dic['4d']['colorFig2'] = 'red'
seq_dic['-4d']['colorFig2'] = 'red'

seq_dic['rk2']['colorFig2'] = 'limegreen'
seq_dic['rk1']['colorFig2'] = 'limegreen'
seq_dic['k2g']['colorFig2'] = 'limegreen'
seq_dic['7k12d']['colorFig2'] = 'limegreen'

seq_dic['-2k-4r5d']['colorFig2'] = 'limegreen'

#NCPR values
seq_dic['wt']['ncpr'] = 0.0579
seq_dic['wtAro']['ncpr'] = 0.0579
seq_dic['k2g']['ncpr'] = 0.0433
seq_dic['aro1']['ncpr'] = 0.0579
seq_dic['allY']['ncpr'] = 0.0578
seq_dic['allF']['ncpr'] = 0.0580
seq_dic['f3y']['ncpr'] = 0.0579
seq_dic['f2y']['ncpr'] = 0.0579
seq_dic['f15y']['ncpr'] = 0.0578
seq_dic['2r']['ncpr'] = 0.0725
seq_dic['7r']['ncpr'] = 0.1090
seq_dic['-10r']['ncpr'] = -0.0151
seq_dic['12d']['ncpr'] = -0.0296
seq_dic['8d']['ncpr'] = 0.000
seq_dic['4d']['ncpr'] = 0.0287
seq_dic['-4d']['ncpr'] = 0.0871
seq_dic['7k12d']['ncpr'] = 0.0214
seq_dic['rk1']['ncpr'] = 0.0579
seq_dic['rk2']['ncpr'] = 0.0579
seq_dic['7r12d']['ncpr'] = 0.0215
seq_dic['-2k-4r5d']['ncpr'] = -0.0151

#Numbers of different residues
seq_dic['2r']['r'] = 12
seq_dic['7r']['r'] = 17
seq_dic['-10r']['r'] = 0
seq_dic['rk1']['r'] = 7
seq_dic['rk2']['r'] = 4
seq_dic['7r12d']['r'] = 17
seq_dic['-2k-4r5d']['r'] = 6

seq_dic['k2g']['k'] = 0
seq_dic['7k12d']['k'] = 9
seq_dic['rk1']['k'] = 5
seq_dic['rk2']['k'] = 8
seq_dic['-2k-4r5d']['k'] = 0

seq_dic['12d']['d'] = 16
seq_dic['8d']['d'] = 12
seq_dic['4d']['d'] = 8
seq_dic['-4d']['d'] = 0
seq_dic['7k12d']['d'] = 16
seq_dic['7r12d']['d'] = 16
seq_dic['-2k-4r5d']['d'] = 9

seq_dic['wtAro']['y'] = 7
seq_dic['aro1']['y'] = 5
seq_dic['allY']['y'] = 19
seq_dic['allF']['y'] = 0
seq_dic['f3y']['y'] = 10
seq_dic['f2y']['y'] = 11
seq_dic['f15y']['y'] = 13
seq_dic['rk1']['y'] = 7
seq_dic['rk2']['y'] = 7

seq_dic['aro1']['f'] = 8
seq_dic['allY']['f'] = 0
seq_dic['allF']['f'] = 19
seq_dic['f3y']['f'] = 3
seq_dic['f2y']['f'] = 4
seq_dic['f15y']['f'] = 3

#Numbers of aromatics
seq_dic['wtAro']['aro'] = 19
seq_dic['aro1']['aro'] = 13
seq_dic['allY']['aro'] = 19
seq_dic['allF']['aro'] = 19
seq_dic['f3y']['aro'] = 13
seq_dic['f2y']['aro'] = 15
seq_dic['f15y']['aro'] = 16
seq_dic['rk1']['aro'] = 19
seq_dic['rk2']['aro'] = 19

#Dilute arm data at different temperatures
seq_dic['wt']['csat']['4'] = 1.10738 * 10 ** -5
seq_dic['wtAro']['csat']['4'] = 1.29 * 10 ** -5
seq_dic['aro1']['csat']['4'] = 0.000339
seq_dic['allY']['csat']['4'] = 0.0000027
seq_dic['allF']['csat']['4'] = np.mean([5.27742 * 10 ** -5,
    4.62913 * 10 ** -5, 4.6695 * 10 ** -5, 5.18412 * 10 ** -5])
seq_dic['f3y']['csat']['4'] = np.mean([0.000113221,
    0.000135638, 0.000111782])
seq_dic['f2y']['csat']['4'] = 6.31483 * 10 ** -5
seq_dic['f15y']['csat']['4'] = 0.00002842
seq_dic['2r']['csat']['4'] = 0.000018
seq_dic['7r']['csat']['4'] = np.mean([0.000160067,
    2.116 * 10 ** -4, 0.000163616])
seq_dic['-10r']['csat']['4'] = np.mean([0.000140478, 0.00013578])
seq_dic['12d']['csat']['4'] = np.mean([0.000108035, 9.91051 * 10 ** -5])
seq_dic['8d']['csat']['4'] = 1.87826 * 10 ** -5
seq_dic['4d']['csat']['4'] = 4.45 * 10 ** -6
seq_dic['-4d']['csat']['4'] = 8.85533 * 10 ** -5
seq_dic['7k12d']['csat']['4'] = np.mean([4.2697 * 10 ** -5,
    4.87509 * 10 ** -5, 0.0000473])
seq_dic['rk1']['csat']['4'] = np.mean([8.87291 * 10 ** -5,
    9.4535 * 10 ** -5, 7.63822 * 10 ** -5, 7.25791 * 10 ** -5])
seq_dic['rk2']['csat']['4'] = np.mean([0.0005, 0.000506871, 0.000472124])
seq_dic['7r12d']['csat']['4'] = 7.74 * 10 ** -7
seq_dic['-2k-4r5d']['csat']['4'] = 0.0001073

seq_dic['wtAro']['csat']['6'] = 1.58 * 10 ** -5
seq_dic['aro1']['csat']['6'] = 0.000421
seq_dic['allY']['csat']['6'] = np.mean([0.0000043, 3.90 * 10 ** -6])
seq_dic['allF']['csat']['6'] = 6.85558 * 10 ** -5
seq_dic['7r']['csat']['6'] = np.mean([0.0002809, 3.360 * 10 ** -4])
seq_dic['12d']['csat']['6'] = 0.000119232
seq_dic['rk2']['csat']['6'] = 0.000663577

seq_dic['wt']['csat']['8'] = 1.31991 * 10 ** -5
seq_dic['wtAro']['csat']['8'] = 2.32 * 10 ** -5
seq_dic['aro1']['csat']['8'] = 0.000512
seq_dic['allY']['csat']['8'] = np.mean([0.0000062,
    5.49 * 10 ** -6, 6.53 * 10 ** -6])
seq_dic['allF']['csat']['8'] = np.mean([6.40081 * 10 ** -5,
    5.62282 * 10 ** -5])
seq_dic['f3y']['csat']['8'] = 0.000211931
seq_dic['7r']['csat']['8'] = np.mean([4.020 * 10 ** -4,
    8.023 * 10 ** -4, 0.000811708])
seq_dic['-10r']['csat']['8'] = 0.000195004
seq_dic['-4d']['csat']['8'] = 0.000194072
seq_dic['7k12d']['csat']['8'] = 8.13758 * 10 ** -5
seq_dic['rk1']['csat']['8'] = np.mean([0.000152051, 0.000157574])
seq_dic['rk2']['csat']['8'] = 0.000922233
seq_dic['7r12d']['csat']['8'] = 1.96682 * 10 ** -6
seq_dic['-2k-4r5d']['csat']['8'] = 0.0001457

seq_dic['wt']['csat']['10'] = np.mean([2.15045 * 10 ** -5,
    2.17511 * 10 ** -5])
seq_dic['wtAro']['csat']['10'] = 2.76606 * 10 ** -5
seq_dic['aro1']['csat']['10'] = 0.000591
seq_dic['allY']['csat']['10'] = np.mean([9.26753 * 10 ** -6,
    1.23455 * 10 ** -5, 8.26 * 10 ** -6])
seq_dic['allF']['csat']['10'] = 0.000100491
seq_dic['f2y']['csat']['10'] = 0.000198444
seq_dic['12d']['csat']['10'] = 0.00020563
seq_dic['7k12d']['csat']['10'] = 1.02629 * 10 ** -4

seq_dic['wt']['csat']['12'] = 2.53169 * 10 ** -5
seq_dic['wtAro']['csat']['12'] = 4.68 * 10 ** -5
seq_dic['aro1']['csat']['12'] = 0.000803
seq_dic['allY']['csat']['12'] = np.mean([1.50791 * 10 ** -5,
    1.37 * 10 ** -5])
seq_dic['allF']['csat']['12'] = np.mean([9.94979 * 10 ** -5, 0.000104539])
seq_dic['f3y']['csat']['12'] = 0.000418548
seq_dic['f15y']['csat']['12'] = 0.000126897
seq_dic['2r']['csat']['12'] = 0.0000482
seq_dic['8d']['csat']['12'] = 0.000079013
seq_dic['4d']['csat']['12'] = 0.00000936
seq_dic['-4d']['csat']['12'] = 0.000316695
seq_dic['7k12d']['csat']['12'] = np.mean([0.000123683, 1.29847 * 10 ** -4])
seq_dic['rk1']['csat']['12'] = np.mean([0.000225759, 0.000308322])
seq_dic['7r12d']['csat']['12'] = 3.62603 * 10 ** -6
seq_dic['-2k-4r5d']['csat']['12'] = 0.0002174

seq_dic['wtAro']['csat']['14'] = 5.29 * 10 ** -5
seq_dic['aro1']['csat']['14'] = 0.000891
seq_dic['allY']['csat']['14'] = np.mean([2.03727 * 10 ** -5,
    1.85 * 10 ** -5, 2.35959 * 10 ** -5])
seq_dic['allF']['csat']['14'] = 0.000164624
seq_dic['rk1']['csat']['14'] = 0.000338969

seq_dic['wt']['csat']['16'] = 4.69146 * 10 ** -5
seq_dic['wtAro']['csat']['16'] = 5.79738 * 10 ** -5
seq_dic['allY']['csat']['16'] = np.mean([2.75851 * 10 ** -5,
    2.73 * 10 ** -5])
seq_dic['allF']['csat']['16'] = 0.000152364
seq_dic['2r']['csat']['16'] = 0.0000955
seq_dic['-10r']['csat']['16'] = 0.000382047
seq_dic['12d']['csat']['16'] = 0.000377364
seq_dic['4d']['csat']['16'] = 0.0000179
seq_dic['7r12d']['csat']['16'] = 5.73266 * 10 ** -6
seq_dic['-2k-4r5d']['csat']['16'] = 0.0002742

seq_dic['wtAro']['csat']['18'] = 7.70 * 10 ** -5
seq_dic['allY']['csat']['18'] = 3.36 * 10 ** -5
seq_dic['7r12d']['csat']['18'] = 6.32 * 10 ** -6

seq_dic['wt']['csat']['20'] = 7.3406 * 10 ** -5
seq_dic['wtAro']['csat']['20'] = 1.02 * 10 ** -4
seq_dic['k2g']['csat']['20'] = np.mean([1.31898 * 10 ** -5,
    1.49119 * 10 ** -5])
seq_dic['allY']['csat']['20'] = np.mean([5.22 * 10 ** -5,
    5.48717 * 10 ** -5])
seq_dic['allF']['csat']['20'] = 0.000209056
seq_dic['f2y']['csat']['20'] = 0.000695682
seq_dic['f15y']['csat']['20'] = 0.000341813
seq_dic['12d']['csat']['20'] = 0.000529269
seq_dic['8d']['csat']['20'] = 0.00013873
seq_dic['7k12d']['csat']['20'] = 0.000270442
seq_dic['rk1']['csat']['20'] = 0.000595078
seq_dic['7r12d']['csat']['20'] = np.mean([9.34 * 10 ** -6,
    9.64765 * 10 ** -6])

seq_dic['wtAro']['csat']['22'] = 0.000124321
seq_dic['7r12d']['csat']['22'] = 1.14 * 10 ** -5

seq_dic['wtAro']['csat']['24'] = 0.000146096
seq_dic['k2g']['csat']['24'] = 2.46749 * 10 ** -5

seq_dic['allF']['csat']['25'] = 0.000452827
seq_dic['7r12d']['csat']['25'] = 1.57 * 10 ** -5

seq_dic['k2g']['csat']['28'] = 3.45973 * 10 ** -5
seq_dic['7r12d']['csat']['28'] = np.mean([1.97 * 10 ** -5,
    2.51 * 10 ** -5])

seq_dic['allY']['csat']['30'] = np.mean([0.000232034,
    0.000257222, 0.0002156])
seq_dic['7r12d']['csat']['30'] = np.mean([3.37 * 10 ** -5,
    3.38 * 10 ** -5])

seq_dic['k2g']['csat']['32'] = 5.4656 * 10 ** -5
seq_dic['7r12d']['csat']['32'] = 4.61 * 10 ** -5

seq_dic['7r12d']['csat']['34'] = 5.74 * 10 ** -5

seq_dic['allY']['csat']['35'] = 0.000472027

seq_dic['7r12d']['csat']['36'] = 6.65 * 10 ** -5

seq_dic['7r12d']['csat']['38'] = 7.12 * 10 ** -5

seq_dic['allY']['csat']['40'] = 0.000889713
seq_dic['7r12d']['csat']['40'] = 0.000104782
