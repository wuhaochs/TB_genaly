#!/home/dell/miniconda3/envs/TB_pip2/bin/python
import os
import subprocess
import pandas as pd 
import argparse
import sys
from Bio import Phylo
import re
import pymysql
pt = re.compile('( {1,})<')
def Treebuild(Sinf,inf,infn,ofn,exbed,csn='goalSample',lin=0,method='ml',boots=1000):
    falist = ['fa','fas','fna','fasta']
    if os.path.isfile(Sinf):
        Sinf = os.path.abspath(Sinf)
    if os.path.isfile(inf):
        inf = os.path.abspath(inf)
    ofn = os.path.abspath(ofn)
    if inf != 'null':
        infn = os.path.abspath(infn)
    if not os.path.isdir(ofn):
        os.makedirs(ofn)
    os.chdir(ofn)
    subprocess.run(f'rm -r *',shell=True)
    if infn != 'null':
        print(infn)
        os.makedirs('upload')
        if os.path.isfile(infn):
            subprocess.run(f'''seqkit seq {infn} > upload/{infn.split('/')[-1]}''',shell=True)
    #判断文件夹内哪些为fa
        elif os.path.isdir(infn):
            fafilelist = []
            for i in os.listdir(infn):
                for tfa in falist:
                    if i.endswith(tfa):
                        Pre=i.replace(f'.{tfa}','')
                        subprocess.run(f'seqkit seq {infn}/{i} > ./upload/{i}',shell=True)
        else:
            print('输入路径错误')
            sys.exit()
    if os.path.isfile(Sinf):
        subprocess.run(f'sort -u {Sinf} > tmplist',shell=True)
    else:
        subprocess.run(f'touch tmplist',shell=True)
    if os.path.isfile(inf):
        open('tmplist','a').write(f'{csn}\t{inf}\n')

    if infn != 'null':
        wkdir = os.getcwd()
        for i in os.listdir('upload'):
            Pre = i.split('.')[0]
            open('tmplist','a').write(f'{Pre}\t{wkdir}/upload/{i}\n')
    query = "SELECT * FROM tcu_database"
    connection = pymysql.connect(host='127.0.0.1',user='baiyi',password='baiyi123@+1s',database='baiyi')
    metadb =  pd.read_sql(query, connection)
    tmetadb = metadb[['name','drtype','lineage','address','age','sex']].fillna('-')
    tmetadb.to_csv('tbprofiler.txt',sep='\t',index=False)
    subprocess.run(f'/data/deploy/TB_soft/other_soft/snippy/bin/snippy-multi tmplist --ref /data/deploy/TB_soft/ref/TB/ref.fa --cpu 10 > runme.sh',shell=True)
    with open('task.log','a') as f:
        subprocess.run(f'sed -i \'s@^@/data/deploy/TB_soft/other_soft/snippy/bin/@g\' runme.sh',shell=True)
        if exbed:
            subprocess.run(f'''sed -i 's@snippy-core@snippy-core --mask {exbed}@g' runme.sh''',shell=True)
        subprocess.run(f'sh runme.sh',shell=True,stdout=f,stderr=f)
        #subprocess.run(f'',shell=True)
        numsam = int(os.popen(f'cat tmplist|wc -l').read().strip())
        if method=='ml':
            if numsam > 4:
                if boots > 0:
                    subprocess.run(f'iqtree2 -s core.aln -m MFP --quiet -b {boots} -T 10',shell=True,stdout=f,stderr=f)
                    if boots ==1:
                        subprocess.run(f'mv core.aln.treefile core.aln.contree',shell=True)
                else:
                    subprocess.run(f'iqtree2 -s core.aln -m MFP --quiet -T 10',shell=True,stdout=f,stderr=f)
                    subprocess.run(f'mv core.aln.treefile core.aln.contree',shell=True)
            else:
                subprocess.run(f'iqtree2 -s core.aln -m MFP --quiet -T 10',shell=True,stdout=f,stderr=f)
                subprocess.run(f'mv core.aln.treefile core.aln.contree',shell=True)

        if method=='nj':
            subprocess.run(f'grapetree -p core.aln --n_proc 10 -m NJ > core.aln.contree',shell=True,stdout=f,stderr=f)
        if method=='grapetree':
            subprocess.run(f'grapetree -p core.aln --n_proc 10  > grapetree.nwk',shell=True,stdout=f,stderr=f)


def annoxml(inf,ofn,meta,meta1=0):
    Phylo.convert(inf,'newick',ofn,'phyloxml')
    mdb = pd.read_table(meta)
    if meta1:
        mdb1 = pd.read_table(meta1)
    with open(ofn) as f:
        for line in f:
            tmpline = line.strip()
            open(f'Anno_{ofn}','a').write(line)
            if tmpline.startswith('<name'):
                node = 1
                snum = pt.findall(line)[0].count('')-1
                name =  tmpline.replace('<','').replace('>','').replace('name','').replace('/','').strip()
                if len(mdb.loc[mdb['毒株名']==name,'家系']):
                    main_lin = mdb.loc[mdb['毒株名']==name,'家系'].tolist()[0]
                    drtype = mdb.loc[mdb['毒株名']==name,'耐药类型'].tolist()[0]
                else:
                    main_lin = 'unknown'
                    drtype = 'unknown'

            if tmpline.startswith('<branch') and node == 1:
                open(f'Anno_{ofn}','a').write(f'{" "*snum}<property ref="ird:main_lin" datatype="xsd:string" applies_to="node">{main_lin}</property>\n')
                open(f'Anno_{ofn}','a').write(f'{" "*snum}<property ref="ird:Dr_type" datatype="xsd:string" applies_to="node">{drtype}</property>\n')
                node = 0 
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='TBpipline')
    parser.add_argument('--inf1','-i1',type=str,default=False,help='选择文件')
    parser.add_argument('--inf2','-i2',type=str,default=False,help='当前样本')
    parser.add_argument('--output','-o',type=str,default='Treebuild',help='输出文件')
    parser.add_argument('--boots','-b',type=int,default=100,help='bootstrap')
    parser.add_argument('--Samname','-s',type=str,default='focus_Sam',help='当前样本名称')
    parser.add_argument('--Lin','-l',type=int,default=0,help='家系名称')
    parser.add_argument('--method','-mt',type=str,default='ml',help='家系名称')
    parser.add_argument('--inf3','-i3',type=str,default=0,help='自上传文件')
    parser.add_argument('--exbed','-eb',type=str,default='/data/deploy/TB_soft/ref/TB/black_list.bed',help='排除位点')
    argv = parser.parse_args()
    inf1 = argv.inf1
    inf2 = argv.inf2
    inf3 = argv.inf3
    ofn = argv.output
    Sam = argv.Samname
    lineage = argv.Lin
    Md = argv.method
    Bs = argv.boots
    eb = argv.exbed
    Treebuild(inf1,inf2,inf3,ofn,eb,Sam,lineage,Md,Bs)
    #if Bs > 1:
    #    annoxml(f'core.aln.contree',f'core.aln.xml','metadb.tsv')
