#!/home/dell/miniconda3/envs/TB_pip2/bin/python
import pandas as pd 
import os 
import subprocess
from sys import argv 
import sys
import argparse
from itertools import combinations
def generate_combinations(lst):
    return list(combinations(lst, 2))
def merge_lists_with_common_elements(lst):
    merged_list = []
    for sublist in lst:
        if not any(set(sublist) & set(item) for item in merged_list):
            merged_list.append(sublist)
        else:
            merged = False
            for i, existing in enumerate(merged_list):
                if set(sublist) & set(existing):
                    merged_list[i] = list(set(sublist) | set(existing))
                    merged = True
            if not merged:
                merged_list.append(sublist)
    return merged_list
def Pair_dis(inf,ofn,nsnp,ref='/data/deploy/TB_soft/ref/TB/ref.fa',binmod=0,dism='TN93'):
	nsnp = int(nsnp)
	#1.core.aln 2.snp-dist 3.fastaANI 3.柱状图矩阵 4.热图矩阵
	if not os.path.isdir(ofn):
		os.makedirs(ofn)
	os.chdir(ofn)
	tmpdb = pd.read_table(inf,names=['Sam','Sam_path'])
	tmpdb[['Sam']].to_csv('Samlist.txt',sep='\t',header=False,index=False)
	Samlist = tmpdb['Sam'].tolist()
	#--- snippy 生成core.aln
	subprocess.run(f'/data/deploy/TB_soft/other_soft/snippy/bin/snippy-multi {inf} --ref {ref} --cpu 10 > runme.sh',shell=True)
	subprocess.run(f'sed -i \'s@^@/data/deploy/TB_soft/other_soft/snippy/bin/@g\' runme.sh',shell=True)
	with open('snp.log','w') as f:
		subprocess.run(f'sh runme.sh',shell=True,stdout=f,stderr=f)
	subprocess.run(f'seqkit grep -v -p Reference core.aln > core_dref.aln',shell=True)
	subprocess.run(f'seqkit grep -v -p Reference core.full.aln > core_dref.full.aln',shell=True)
	subprocess.run(f'rm split*.aln',shell=True)
	subprocess.run(f'''seqkit split -i core_dref.full.aln -O ./ --by-id-prefix 'split_full' ''',shell=True)
	subprocess.run(f'''seqkit split -i core_dref.aln -O ./ --by-id-prefix 'split_core' ''',shell=True)
	if os.path.isfile('genomelist.txt'):
		subprocess.run(f'rm genomelist.txt',shell=True)
	if os.path.isfile('genome_corelist.txt'):
		subprocess.run(f'rm genome_corelist.txt',shell=True)
	for i in os.listdir():
		if i.startswith('split_full') and i.endswith('aln'):
			open('genomelist.txt','a').write(f'{os.getcwd()}/{i}\n')
		elif i.startswith('split_core') and i.endswith('aln'):
			open('genome_corelist.txt','a').write(f'{os.getcwd()}/{i}\n')
	subprocess.run(f'snp-dists core_dref.aln > dis.mat.txt',shell=True)
	subprocess.run(f'snp-dists -m core_dref.aln > t_dis.tsv',shell=True)
	df = pd.read_table('t_dis.tsv',names=['A','B','D'])
	df = df[df['A']!=df['B']]
	df['AB_tuple'] = df[['A', 'B']].apply(lambda x:sorted(list(x)), axis=1)
	df = df[~df['AB_tuple'].duplicated()]
	if not binmod:
		#---0;1-10;10-100;100-1000;1000+分档计数
		disdict = {'0':0,'1-10':0,'10-100':0,'100-1000':0,'1000+':0}
		disdict['0'] = len([i for i in df['D'].tolist() if i == 0])
		disdict['1-10'] = len([i for i in df['D'].tolist() if i > 0 and i <= 10])
		disdict['10-100'] = len([i for i in df['D'].tolist() if i > 10 and i <= 100])
		disdict['100-1000'] = len([i for i in df['D'].tolist() if i > 100 and i <= 1000])
		disdict['1000+'] = len([i for i in df['D'].tolist() if i > 1000])
		
	else:
		disdict = {}
		binlist =  sorted([int(i) for i in binmod.split(',')])
		if len(binlist) >= 2:
			for modnum in range(len(binlist)):
				mod = int(binlist[modnum])
				if modnum == 0 :
					disdict[f'<={mod}'] = len([i for i in df['D'].tolist() if i <= mod])
				elif modnum == len(binlist)-1:
					oldmod = int(binlist[modnum-1])
					disdict[f'{oldmod}-{mod}'] = len([i for i in df['D'].tolist() if i > oldmod and i <= mod])
					disdict[f'{mod}+'] = len([i for i in df['D'].tolist() if i > mod])
				else:
					oldmod = int(binlist[modnum-1])
					disdict[f'{oldmod}-{mod}'] = len([i for i in df['D'].tolist() if i > oldmod and i <= mod])
		else:
			mod = binlist[0]
			disdict[f'<={mod}'] = len([i for i in df['D'].tolist() if i <= mod])
			disdict[f'>{mod}'] = len([i for i in df['D'].tolist() if i > mod])

	disdf = pd.DataFrame(disdict,index=['0'])
	disdf.to_csv('dis_bin.tsv',sep='\t',index=False)
	# 根据差异位点数量和差异位点阈值进行聚类
	Clusterlist = []
	for tmpi in df.index:
		tmpdf = df.loc[df.index==tmpi,:]
		if tmpdf['D'].tolist()[0] <= nsnp:
			Clusterlist.append([tmpdf['A'].tolist()[0],tmpdf['B'].tolist()[0]])
	merlist = merge_lists_with_common_elements(Clusterlist)
	#print(merlist)
	Clusternum = 1
	if os.path.isfile('Cluster.tsv'):
		subprocess.run(f'rm Cluster.tsv',shell=True)
	open('Cluster.tsv','w').write('聚类名称\t样本数量\t聚类样本\t最大snp差异\t最小snp差异\t平均snp差异\n')
	for Cluster in merlist:
		snplist = []
		All_clu = generate_combinations(Cluster)
		for Sam1,Sam2 in All_clu:	
			tmpdisdf = df.loc[((df['A']==Sam1) & (df['B']==Sam2))|((df['A']==Sam2) & (df['B']==Sam1)),:]
			snplist.append(tmpdisdf['D'].tolist()[0])
		open('Cluster.tsv','a').write(f'''Cluster{Clusternum}\t{len(Cluster)}\t{','.join(Cluster)}\t{max(snplist)}\t{min(snplist)}\t{sum(snplist)/len(snplist)}\n''')
		Clusternum+=1
	# fastANI mat
	'''
	'''
	#--- 计算距离矩阵 提供snp 和 TN93 ANI三种方式
	if dism in ['TN93','SNP']:
		subprocess.run(f'java -jar /data/deploy/TB_soft/script/SeqRuler.jar -i core_dref.aln -d {dism} -o t_Gdis.txt -a average -c 10',shell=True)
		misdf = pd.read_table('t_Gdis.txt',sep=',')
		misdf.to_csv('Gdis.txt',sep='\t',index=False)
		rmisdf = misdf[['Target','Source','Distance']]
		rmisdf.columns =  ['Source','Target','Distance']
		misdf = pd.concat([misdf,rmisdf]).reset_index(drop=True)
		for Samnum in range(len(Samlist)):
			tSam = Samlist[Samnum]
			misdf.loc[len(misdf)] = [tSam,tSam,0]
		misdf['xlab'] = '-'
		misdf['ylab'] = '-'
		for tindex in misdf.index:
			xSam = misdf.loc[misdf.index==tindex,'Source'].tolist()[0]
			ySam = misdf.loc[misdf.index==tindex,'Target'].tolist()[0]
			xSamloc = Samlist.index(xSam)
			ySamloc = Samlist.index(ySam) 
			misdf.loc[misdf.index==tindex,'xlab'] = xSamloc 
			misdf.loc[misdf.index==tindex,'ylab'] = ySamloc 

		misdf[['Source','Target','Distance','xlab','ylab']].to_csv('Gdis_core.txt',sep='\t',index=False)


	elif dism == 'ANI':
		if int(os.popen('seqkit stat core.full.aln -T |cut -f5|tail -n1').read().strip()) > 1000000:
			subprocess.run(f'/data/deploy/TB_soft/other_soft/fastANI --rl genomelist.txt  --ql genomelist.txt -o Full_ANI.txt --matrix -t 10',shell=True)
		else:
			subprocess.run(f'/data/deploy/TB_soft/other_soft/fastANI --rl genomelist.txt  --ql genomelist.txt -o Full_ANI.txt --matrix -t 10 --fragLen 100 ',shell=True)
		misdf = pd.read_table('Full_ANI.txt',names=['Source','Target','Distance','A','B'])
		misdf['Source'] = misdf['Source'].str.split('/').str[-1].str.replace('split_full','').str.replace('.aln','')
		misdf['Target'] = misdf['Target'].str.split('/').str[-1].str.replace('split_full','').str.replace('.aln','')
		misdf = misdf[misdf['Source']!=misdf['Target']]
		misdf['AB_tuple'] = misdf[['Source', 'Target']].apply(lambda x:sorted(list(x)), axis=1)
		misdf = misdf[~misdf['AB_tuple'].duplicated()]
		misdf = misdf[['Source','Target','Distance']]
		misdf.to_csv('Gdis.txt',sep='\t',index=False)
		rmisdf = misdf[['Target','Source','Distance']]
		rmisdf.columns =  ['Source','Target','Distance']
		misdf = pd.concat([misdf,rmisdf]).reset_index(drop=True)
		
		for Samnum in range(len(Samlist)):
			tSam = Samlist[Samnum]
			misdf.loc[len(misdf)] = [tSam,tSam,100]
		misdf['xlab'] = '-'
		misdf['ylab'] = '-'
		for tindex in misdf.index:
			xSam = misdf.loc[misdf.index==tindex,'Source'].tolist()[0]
			ySam = misdf.loc[misdf.index==tindex,'Target'].tolist()[0]
			xSamloc = Samlist.index(xSam)
			ySamloc = Samlist.index(ySam) 
			misdf.loc[misdf.index==tindex,'xlab'] = xSamloc 
			misdf.loc[misdf.index==tindex,'ylab'] = ySamloc 

		misdf[['Source','Target','Distance','xlab','ylab']].to_csv('Gdis_full.txt',sep='\t',index=False)

		#misdf[['Source','Target','Distance']].to_csv('Gdis_full.txt',sep='\t',index=False)


		if int(os.popen('seqkit stat core.aln -T |cut -f5|tail -n1').read().strip()) > 1000000:
			subprocess.run(f'/data/deploy/TB_soft/other_soft/fastANI --rl genome_corelist.txt  --ql genome_corelist.txt -o Core_ANI.txt --matrix -t 10',shell=True)
		else:
			subprocess.run(f'/data/deploy/TB_soft/other_soft/fastANI --rl genome_corelist.txt  --ql genome_corelist.txt -o Core_ANI.txt --matrix -t 10  --fragLen 100',shell=True)
		misdf = pd.read_table('Core_ANI.txt',names=['Source','Target','Distance','A','B'])
		misdf['Source'] = misdf['Source'].str.split('/').str[-1].str.replace('split_core','').str.replace('.aln','')
		misdf['Target'] = misdf['Target'].str.split('/').str[-1].str.replace('split_core','').str.replace('.aln','')
		misdf = misdf[misdf['Source']!=misdf['Target']]
		misdf['AB_tuple'] = misdf[['Source', 'Target']].apply(lambda x:sorted(list(x)), axis=1)
		misdf = misdf[~misdf['AB_tuple'].duplicated()]
		misdf = misdf[['Source','Target','Distance']]
		misdf.to_csv('Gdis.txt',sep='\t',index=False)

		rmisdf = misdf[['Target','Source','Distance']]
		rmisdf.columns =  ['Source','Target','Distance']
		misdf = pd.concat([misdf,rmisdf]).reset_index(drop=True)
		for Samnum in range(len(Samlist)):
			tSam = Samlist[Samnum]
			misdf.loc[len(misdf)] = [tSam,tSam,100]
		misdf['xlab'] = '-'
		misdf['ylab'] = '-'
		for tindex in misdf.index:
			xSam = misdf.loc[misdf.index==tindex,'Source'].tolist()[0]
			ySam = misdf.loc[misdf.index==tindex,'Target'].tolist()[0]
			xSamloc = Samlist.index(xSam)
			ySamloc = Samlist.index(ySam) 
			misdf.loc[misdf.index==tindex,'xlab'] = xSamloc 
			misdf.loc[misdf.index==tindex,'ylab'] = ySamloc
		misdf[['Source','Target','Distance','xlab','ylab']].to_csv('Gdis_core.txt',sep='\t',index=False)

		#misdf[['Source','Target','Distance']].to_csv('Gdis_core.txt',sep='\t',index=False)














if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Snp_distpip')
	parser.add_argument('--inputlist','-l',type=str,default=False,help='文件列表')
	parser.add_argument('--output','-o',type=str,default=False,help='输出文件')
	parser.add_argument('--transcutoff','-u',type=str,default='5',help='传播分析设置的snp数量阈值')
	parser.add_argument('--snpdis','-s',type=str,default='10,20,50,100,200,500,1000',help='snp区段分布')
	parser.add_argument('--method','-m',type=str,default='TN93',help='计算相关性方法')
	parser.add_argument('--ref','-r',type=str,default='/data/deploy/TB_soft/ref/TB/ref.fa',help='计算相关性方法')
	argv = parser.parse_args()
	infile = argv.inputlist
	outfile = argv.output
	ClusterNum = argv.transcutoff
	Dislist = argv.snpdis
	Dismeth = argv.method
	Ref = argv.ref
	Pair_dis(infile,outfile,ClusterNum,Ref,Dislist,Dismeth)
