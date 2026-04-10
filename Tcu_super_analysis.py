#!/home/dell/miniconda3/envs/TB_pip2/bin/python
# -*- coding: utf-8 -*-
#--------env Tb_pip------------
import sys
import os
import argparse
import datetime
import subprocess
import json
import pandas as pd
import re
import time
from Bio.Seq import Seq
import math
import pysam
import pymysql
__author__='wsh'
__version__='1.6.6'
pd.options.mode.chained_assignment = None
#------------------
#1.3 不用snippy 来获得一致性序列和vcf文件，写脚本来完成
#1.4 更新mst算法
#1.4 fast 去掉两步fastqc 增加分支杆菌单独的库，增加vcf毒力注释文件 增加iqtree 和 fasttree本地绘图部分
#1.5 增加三代数据处理通过传入包含多fq的文件夹或者多fast5的文件夹
#1.5 增加fasta接口
#1.6 增加clockwork算法
#1.6.1 计算物种成分由kraken2 更换成kraken2+braken
#1.6.2 更新同一密码子内多个突变的问题
#1.6.3 密码子矫正bug处理，耐药矫正
#1.6.5 可重复运行
#1.6.6 修改耐药注释逻辑 增加 vcf和bam作为输入
#20241112 遇见有问题的fastq可以滤过运行后面
#20241125 遇见有问题的tbprofilerbug
#20241126 vcf bam fasta在用药指导的bug
#------------------
parser = argparse.ArgumentParser(description='TBpipline')
parser.add_argument('--list','-l',type=str,default=False,help='文件列表')
parser.add_argument('--thread','-t',type=int,default=10,help='线程数量')
parser.add_argument('--output','-o',type=str,default=False,help='输出文件')
parser.add_argument('--transcutoff','-u',type=str,default='5,12',help='传播分析设置的snp数量阈值')
parser.add_argument('--vcf_dpt','-v',type=int,default=5,help='测序深度')
parser.add_argument('--db','-d',type=str,default='TB_library',help='Standard.标准库;TB_library.TB库;GT_DB.微生物全库,分析时间较慢')
parser.add_argument('--mode','-m',type=str,default='hac',help='fast5文件 guppy转换config选择 hac or fast')
parser.add_argument('--maf','-mf',type=float,default=0.25,help='次等位基因频率')
parser.add_argument('--clockwork','-ck',type=str,default=False,help='是否利用clockwork进行snp callling')
argv = parser.parse_args()
Tims = datetime.datetime.now().strftime('%m%d_%H%M')
cutoff = argv.transcutoff
Mode = argv.mode
maf_cutoff = argv.maf
vcfdep = argv.vcf_dpt
stime=time.time()
if argv.clockwork:
    if argv.clockwork.upper() in ['FALSE','0','F'] :
        ck = False
    else:
        ck = True
print(f'clockwork 是否启动:{ck}')
def Pro_fun(a,b):
    print(f'{b} analysis START')
    ttime1=time.time()
    newrun = exec(a)
    ttime2=time.time()
    runtime = ttime2 - ttime1
    print(f'{b} analysis END;run time {runtime:0.2f}s')

#-----chech_file_dir-------
def check_file(nfile):
    if not os.path.isfile(nfile):
        raise Exception(f'{nfile}不是列表文件,或者该文件不存在请仔细核对')
    else:
        return(os.path.abspath(nfile))

def check_dir(ndir):
    if not os.path.isdir(ndir):
        raise Exception(f'{ndir}不是目录，或者该目录不存在请仔细核对')
    else:
        return(os.path.abspath(ndir))

def make_dir(ndir):
    if not os.path.exists(ndir):
        os.mkdir(ndir)
        os.chdir(ndir)
    else:
        print('路径已存在，请核实')
        os.chdir(ndir)

def format_seconds(seconds):
    # 将秒数转换为整数
    total_seconds = int(seconds) 
    # 计算小时、分钟和秒
    hours = total_seconds // 3600
    minutes = (total_seconds % 3600) // 60
    seconds = total_seconds % 60
    # 格式化输出
    return f"{hours}小时{minutes}分钟{seconds}秒"
#---------路径设置------
Listn = check_file(argv.list)
make_dir(argv.output)
wkdir=os.getcwd()
nT = argv.thread
optstd = f'''输入文件:\t{Listn}
输出文件:\t{argv.output}
线程数:\t{nT}
'''
print(optstd)
sys.stdout.flush()
ref='/data/deploy/TB_soft/ref/TB/ref.fa'
drtab = '/data/deploy/TB_soft/ref/TB/newDr2.tsv'
DrugGfile = '/data/deploy/TB_soft/ref/TB/Drug_guide.tsv'
#ref2='/data/mlf/miniconda3/envs/TB/share/tbprofiler/tbdb.fasta'
snippydir='/data/deploy/TB_soft/other_soft/snippy/'
if argv.db=="TB_library":
    Krdb='/data/deploy/TB_soft/Database/TB_base'
else:
    Krdb='/data/deploy/TB_soft/Database/All_base'
sc1='/data/deploy/TB_soft/other_soft/3_kreport2krona.py'
listt='/data/deploy/TB_soft/other_soft/taxa_info_0717.txt'
gff = '/data/deploy/TB_soft/ref/TB/ref.gff'
config = '/data/deploy/TB_soft/ref/TB/snpeff.config'
black_list = "/data/deploy/TB_soft/ref/TB/black_list.bed"
MICfile='/data/deploy/TB_soft/ref/TB/DR_MIC_new.tsv'
#----------定义所有工作函数-----------
#-----参考基因组准备------

def ref_process(ref,gff,rwkdir):
    #---------准备ref文件------
    os.chdir(rwkdir)
    if not os.path.isfile('./ref.fa'):
        subprocess.run(f'ln -s {ref} ./ref.fa',shell=True)
    subprocess.run(f'samtools faidx ref.fa -o ./ref.fa.fai',shell=True)
    subprocess.run(f'bwa index ref.fa > {wkdir}/task.log 2>&1',shell=True)
    new_ref = os.path.abspath('ref.fa')
    subprocess.run(f'{snippydir}/binaries/noarch/fasta_generate_regions.py ref.fa.fai 200000 > ./ref.txt',shell=True)
    make_dir('reference')
    if not os.path.exists('ref') and not os.path.exists('genomes'):
        os.mkdir('ref')
        os.mkdir('genomes')
    if not os.path.isfile(f'./ref/genes.gff'):
        subprocess.run(f'ln -s {gff} ./ref/genes.gff ',shell=True)
    if not os.path.isfile(f'./genomes/ref.fa'):
        subprocess.run(f'ln -s {new_ref} ./genomes/ref.fa ',shell=True)
    subprocess.run(f'cp {config} ./',shell=True)
        #-------准备gff文件------
    subprocess.run(f'/home/dell/biosoft/snpEff/scripts/snpEff build -c snpeff.config -dataDir . ref -gff3 -q ',shell=True)
#----------fastqc--------
def qc_funt(Out,Pre,nT,read1,read2=0):
    code1=f'''fastqc {read1} -t {nT} -o {Out}
    unzip   {Out}/{Pre}_1_fastqc.zip -d {Out}
    convert -append \\
        {Out}/{Pre}_1_fastqc/Images/per_base_quality.png \\
        {Out}/{Pre}_1_fastqc/Images/per_sequence_gc_content.png \\
        {Out}/{Pre}_1A.png

    convert -append \\
        {Out}/{Pre}_1_fastqc/Images/per_base_sequence_content.png \\
        {Out}/{Pre}_1_fastqc/Images/adapter_content.png \\
        {Out}/{Pre}_1B.png

    convert -append \\
        {Out}/{Pre}_1_fastqc/Images/sequence_length_distribution.png \\
        {Out}/{Pre}_1_fastqc/Images/duplication_levels.png \\
        {Out}/{Pre}_1C.png

    convert +append {Out}/{Pre}_1A.png {Out}/{Pre}_1B.png {Out}/{Pre}_1D.png
    convert +append {Out}/{Pre}_1C.png {Out}/{Pre}_1D.png {Out}/{Pre}_1E.png

    convert {Out}/{Pre}_1E.png       \\
        -background White -pointsize 96        \\
        -gravity center label:'{Pre}_1'     \\
        +swap                                  \\
        -append                                \\
        {Out}/FastQC_{Pre}_1.png
    '''
    code2=f'''fastqc {read2} -t {nT} -o {Out}
    echo "what happpened"
    unzip   {Out}/{Pre}_2_fastqc.zip -d {Out}
    convert -append \\
        {Out}/{Pre}_2_fastqc/Images/per_base_quality.png \\
        {Out}/{Pre}_2_fastqc/Images/per_sequence_gc_content.png \\
        {Out}/{Pre}_2A.png

    convert -append \\
        {Out}/{Pre}_2_fastqc/Images/per_base_sequence_content.png \\
        {Out}/{Pre}_2_fastqc/Images/adapter_content.png \\
        {Out}/{Pre}_2B.png

    convert -append \\
        {Out}/{Pre}_2_fastqc/Images/sequence_length_distribution.png \\
        {Out}/{Pre}_2_fastqc/Images/duplication_levels.png \\
        {Out}/{Pre}_2C.png

    convert +append {Out}/{Pre}_2A.png {Out}/{Pre}_2B.png {Out}/{Pre}_2D.png
    convert +append {Out}/{Pre}_2C.png {Out}/{Pre}_2D.png {Out}/{Pre}_2E.png

    convert {Out}/{Pre}_2E.png       \\
        -background White -pointsize 96        \\
        -gravity center label:'{Pre}_2'     \\
        +swap                                  \\
        -append                                \\
        {Out}/FastQC_{Pre}_2.png
    '''
    
    code3=f'''convert +append {Out}/FastQC_{Pre}_1.png {Out}/FastQC_{Pre}_2.png {Out}/FastQC_{Pre}.png
    
    rm -r {Out}/{Pre}_*_fastqc/
    rm {Out}/FastQC_{Pre}_*.png
    '''
    code4=f'''fastqc {read1} -t {nT} -o {Out}
    unzip   {Out}/{Pre}_fastqc.zip -d {Out}
    convert -append \\
        {Out}/{Pre}_fastqc/Images/per_base_quality.png \\
        {Out}/{Pre}_fastqc/Images/per_sequence_gc_content.png \\
        {Out}/{Pre}_A.png

    convert -append \\
        {Out}/{Pre}_fastqc/Images/per_base_sequence_content.png \\
        {Out}/{Pre}_fastqc/Images/adapter_content.png \\
        {Out}/{Pre}_B.png

    convert -append \\
        {Out}/{Pre}_fastqc/Images/sequence_length_distribution.png \\
        {Out}/{Pre}_fastqc/Images/duplication_levels.png \\
        {Out}/{Pre}_C.png

    convert +append {Out}/{Pre}_A.png {Out}/{Pre}_B.png {Out}/{Pre}_D.png
    convert +append {Out}/{Pre}_C.png {Out}/{Pre}_D.png {Out}/{Pre}_E.png

    convert {Out}/{Pre}_E.png       \\
        -background White -pointsize 96        \\
        -gravity center label:'{Pre}'     \\
        +swap                                  \\
        -append                                \\
        {Out}/FastQC_{Pre}.png
    '''
    if read2:
        subprocess.run(code1,shell=True)
        subprocess.run(code2,shell=True)
        subprocess.run(code3,shell=True)
    else:
        subprocess.run(code4,shell=True)
    print(f'{Pre}样本质控完成')

#------trim-reads--------
def fa_process(fa,outdir,samplename):
    #----生成30x左右数据-s step -W windowsize
    subprocess.run(f'seqkit sliding -s 5 -W 150 {fa} -w 0 > {outdir}/{samplename}.fa',shell=True)
    ofn = f'{outdir}/{samplename}.fq'
    with open(f'{outdir}/{samplename}.fa') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                readsID=line.split('>')[1]
                header=f'@{readsID}'
                open(ofn,'a').write(f'{header}\n')
            else:
                reads=line
                readslen=len(reads)
                qreads='J'*readslen
                newline=f'''{reads}
+
{qreads}
'''
                open(ofn,'a').write(f'{newline}')


def trim_fun(outdir,threads,samplename,gs,read1,read2=0):
    clean_read1 = f'{samplename}_clean_1.gz'
    clean_read2 = f'{samplename}_clean_2.gz'
    if read2:
        subprocess.run(f'''fastp --in1 {read1} \\
                --out1 {outdir}/{clean_read1} \\
                --in2 {read2} \\
                --out2 {outdir}/{clean_read2} \\
                --thread {threads} \\
                --length_required=50 \\
                --n_base_limit=6 \\
                --compression=6 \\
                --detect_adapter_for_pe \\
                --json {outdir}/{samplename}.fastp.json \\
                2> {outdir}/{samplename}.fastp.log
                
                ''',shell=True )
        if not os.popen(f'head -n1 {outdir}/{samplename}.fastp.log').read().strip().startswith('ERROR'):
            with open(f'{wkdir}/trim_fqlist','a') as f:
                f.write(f'{Pre}_{gs}\t{outdir}/{clean_read1}\t{outdir}/{clean_read2}\tfastq\n')
            if not os.path.isfile(f'{wkdir}/fq_file/{samplename}_OK'):
                open(f'{wkdir}/Samplelist.txt','a').write(f'{Pre}\n')
                open(f'{wkdir}/fq_file/{samplename}_OK','w').write('OK')

    else:
        clean_read = f'{samplename}_clean.gz'
        subprocess.run(f'''fastp --in1 {read1} \\
                --out1 {outdir}/{clean_read} \\
                --thread {threads} \\
                --compression=6 \\
                -Q \\
                -A \\
                --json {outdir}/{samplename}.fastp.json \\
                2>{outdir}/{samplename}.fastp.log

                ''',shell=True)
        if not os.popen(f'head -n1 {outdir}/{samplename}.fastp.log').read().strip().startswith('ERROR'):
            with open(f'{wkdir}/trim_fqlist','a') as f:
                f.write(f'{Pre}_{gs}\t{outdir}/{clean_read}\t-\tfastq\n')
            if not os.path.isfile(f'{wkdir}/fq_file/{samplename}_OK'):
                open(f'{wkdir}/Samplelist.txt','a').write(f'{Pre}\n')
                open(f'{wkdir}/fq_file/{samplename}_OK','w').write('OK')



        
#-----trans vcf or bam Chrom
def modify_bam_reference(input_bam, output_bam, new_reference_name):
    # 打开输入 BAM 文件
    with pysam.AlignmentFile(input_bam, "rb") as in_bam:
        # 创建一个新的 BAM header
        new_header = in_bam.header.to_dict()

        # 修改参考基因组名称
        for seq in new_header['SQ']:
            seq['SN'] = new_reference_name

        # 打开输出 BAM 文件
        with pysam.AlignmentFile(output_bam, "wb", header=new_header) as out_bam:
            # 遍历输入 BAM 文件中的所有记录
            for read in in_bam:
                # 将记录写入新的 BAM 文件
                out_bam.write(read)


def modify_vcf_reference(input_vcf, output_vcf, new_reference_name):
    # 打开输入 VCF 文件
    with pysam.VariantFile(input_vcf, 'r') as in_vcf:
        # 创建一个新的 VCF header
        new_header = in_vcf.header.copy()

        # 修改参考基因组名称
        open('rename.txt','w').write('')
        for contig in new_header.contigs:
            open('rename.txt','a').write(f'{new_header.contigs[contig].name}\tChromosome\n')
        subprocess.run(f'bcftools annotate --rename-chrs rename.txt {input_vcf} > {output_vcf}',shell=True)

        # 打开输出 VCF 文

#-----call_snp and Drug guide--------
drug_dict = {'AMK':'阿米卡星:amikacin','amikacin':'阿米卡星:amikacin','BDQ':'贝达喹啉:bedaquiline','bedaquiline':'贝达喹啉:bedaquiline','CAP':'卷曲霉素:capreomycin','capreomycin':'卷曲霉素:capreomycin','CFZ':'氯法齐明:Clofazimine','Clofazimine':'氯法齐明:Clofazimine','DLM':'德拉马尼:delamanid','delamanid':'德拉马尼:delamanid','EMB':'乙胺丁醇:ethambutol','ethambutol':'乙胺丁醇:ethambutol','ETO':'乙硫异烟胺:ethionamide','ethionamide':'乙硫异烟胺:ethionamide','FQ':'氟喹诺酮:fluoroquinolone','fluoroquinolone':'氟喹诺酮:fluoroquinolone','INH':'异烟肼:isoniazid','isoniazid':'异烟肼:isoniazid','KAN':'卡那霉素:kanamycin','kanamycin':'卡那霉素:kanamycin','LFX':'左氧氟沙星:levofloxacin','levofloxacin':'左氧氟沙星:levofloxacin','LZD':'利奈唑胺:linezolid','linezolid':'利奈唑胺:linezolid','MFX':'莫西沙星:moxifloxacin','moxifloxacin':'莫西沙星:moxifloxacin','OFX':'氧氟沙星:ofloxacin','ofloxacin':'氧氟沙星:ofloxacin','PTO':'丙硫异烟胺:prothionamide','prothionamide':'丙硫异烟胺:prothionamide','PZA':'吡嗪酰胺:pyrazinamide','pyrazinamide':'吡嗪酰胺:pyrazinamide','rifampicin':'利福平:rifampicin','RIF':'利福平:rifampicin','STM':'链霉素:streptomycin','streptomycin':'链霉素:streptomycin','LEV':'左氧氟沙星:levofloxacin','levofloxacin':'左氧氟沙星:levofloxacin','AMI':'阿米卡星:amikacin','amikacin':'阿米卡星:amikacin','MXF':'莫西沙星:moxifloxacin','moxifloxacin':'莫西沙星:moxifloxacin','ETH':'乙硫异烟胺:ethionamide','ethionamide':'乙硫异烟胺:ethionamide','cycloserine':'环丝氨酸:cycloserine','CYS':'环丝氨酸:cycloserine'}
Firstline = ['吡嗪酰胺','异烟肼','乙胺丁醇','利福平']
#Secondline = ['乙硫异烟胺','利奈唑胺','氯法齐明','贝达喹啉','氟喹诺酮','左氧氟沙星','莫西沙星','阿米卡星','丙硫异烟胺']
def Drugvcf(annofile):
    micdb = pd.read_table('/data/deploy/TB_soft/ref/TB/MIC_edit.tsv')
    metadb1 = pd.read_table('/data/deploy/TB_soft/ref/TB/TB_GeneName.tsv')
    query = "SELECT * FROM tb_dr_database"
    connection = pymysql.connect(host='127.0.0.1',user='baiyi',password='baiyi123@+1s',database='baiyi')
    df = pd.read_sql(query, connection)
    df = df[['id','drug','gene','mut','mut_two','assow','raw_geneid','gene_name']]
    metadb2 = df.rename(columns={'id':'Index','drug':'Drug','mut':'Mut','mut_two':'Mut2','assow':'Assow','raw_geneid':'RawGeneID','gene_name':'GeneName'})
    metadb2 = metadb2[(~metadb2['Assow'].isna()) & (~metadb2['Drug'].isna())]
    #metadb2 = pd.read_table('/data/deploy/TB_soft/ref/TB/TB_dr_edit.tsv')
    def getMIC(tmpdb):
        if micdb.loc[(micdb['new']==tmpdb['Mut2']) & (micdb['Drug']==tmpdb['Drug'])].shape[0] != 0:
            print(tmpdb['Drug'])
            MIC = micdb[micdb['new']==tmpdb['Mut2']]['log2MIC'].tolist()[0]
        else:
            MIC = '-'
        return MIC
    drug_dict = {'Amikacin':'阿米卡星:amikacin','Bedaquiline':'贝达喹啉:bedaquiline','Capreomycin':'卷曲霉素:capreomycin','CFZ':'氯法齐明:clofazimine','Delamanid':'德拉马尼:delamanid','Ethambutol':'乙胺丁醇:ethambutol','Ethionamide':'乙硫异烟胺:ethionamide','FQ':'氟喹诺酮:fluoroquinolone','Isoniazid':'异烟肼:isoniazid','Kanamycin':'卡那霉素:kanamycin','Levofloxacin':'左氧氟沙星:levofloxacin','Linezolid':'利奈唑胺:linezolid','Moxifloxacin':'莫西沙星:moxifloxacin','OFX':'氧氟沙星:ofloxacin','PTO':'丙硫异烟胺:prothionamide','Pyrazinamide':'吡嗪酰胺:pyrazinamide','Rifampicin':'利福平:rifampicin','Streptomycin':'链霉素:streptomycin','LEV':'左氧氟沙星:levofloxacin','AMI':'阿米卡星:amikacin','MXF':'莫西沙星:moxifloxacin','ETH':'乙硫异烟胺:ethionamide','CYS':'环丝氨酸:cycloserine','Clofazimine':'氯法齐明:Clofazimine'}
    drug_dict2 = {'RIF':'Rifampicin','EMB':'Ethambutol','ETH':'Ethionamide','INH':'Isoniazid','KAN':'Kanamycin','MXF':'Moxifloxacin','SZ':'Streptomycin'}
    subprocess.run(f'''grep -v "^#" {annofile} | cut -f2,8  > annotated_file.txt''',shell=True)
    subprocess.run(f'''python /data/deploy/TB_soft/other_soft/load_annotation_variants.py --file annotated_file.txt''',shell=True)
    vcfsknum = int(os.popen(f'''grep '##' {annofile} |wc -l''').read())
    tmpvcfdb = pd.read_table(annofile,skiprows=vcfsknum)
    SamN = tmpvcfdb.columns[-1]
    tmpvcfdb['altn'] = tmpvcfdb[SamN].str.split(':').str[2].str.split(',').str[1]
    tmpvcfdb['refn'] = tmpvcfdb[SamN].str.split(':').str[2].str.split(',').str[0]
    tmpvcfdb = tmpvcfdb[['POS','REF','ALT','altn','refn']]
    rawvcfdb = pd.read_table('annotation_normalized.tsv',header=None)
    rawvcfdb['GeneID'] = rawvcfdb.apply(lambda x:x[2].replace('gene:','').split('-')[0].replace('GENE_','').replace('EBG00000313325','rrs').replace('EBG00000313339','rrl').replace('EBG00000313349','rrf'),axis=1)
    rawvcfdb = rawvcfdb.merge(metadb1,left_on='GeneID',right_on='RawGeneID')
    rawvcfdb['Mut2'] = rawvcfdb['GeneName'] + '_' + rawvcfdb[5]
    drugdb = rawvcfdb.merge(metadb2,on='Mut2')
    drugdb['染色体'] = 'Chromosome'
    drugdb = drugdb.merge(tmpvcfdb,left_on = 0,right_on = 'POS')
    if drugdb.shape[0]>0:
        drugdb['变异频率'] = drugdb.apply(lambda x :f'''{x['REF']}:{x['refn']} {x['ALT']}:{x['altn']} ({round(float(x['altn'])/(float(x['altn'])+float(x['refn']))*100,2)}%)'''  if (float(x['altn'])+ float(x['refn'])) != 0 else '100%',axis=1)
        #drugdb.to_csv('tmp111.tsv',sep='\t',index=False)
        drugdb['药物名'] =  drugdb.apply(lambda x:drug_dict.get(x['Drug'].capitalize(),x['Drug']),axis=1)
        drugdb.fillna('-',inplace=True)
        drugdb['MIC'] = drugdb.apply(lambda x:getMIC(x),axis=1)
        drugdb = drugdb[['染色体','POS','REF','ALT','变异频率','GeneName_x',1,'药物名','Assow','Mut2','MIC']]
        drugdb.rename(columns={'Mut2':'突变结果','Assow':'相关性','POS':'变异位置','GeneName_x':'基因名',1:'变异效应','REF':'参考碱基','ALT':'突变碱基'},inplace=True)
        drugdb = drugdb[['染色体','变异位置','参考碱基','突变碱基','变异频率','基因名','变异效应','药物名','相关性','突变结果','MIC']]
        drugdb = drugdb.drop_duplicates()
        drugdb.to_csv('Chin_snpdr.tsv',sep='\t',index=False)
        MIC_tab = drugdb.rename(columns={'染色体':'CHROM','变异位置':'POS','参考碱基':'REF','突变碱基':'ALT','变异频率':'EVIDENCE','基因名':'Gene','变异效应':'EFFECT','药物名':'Drug','相关性':'Corr','突变结果':'EFFECT2'})
        MIC_tab = MIC_tab[['CHROM','POS','REF','ALT','EVIDENCE','Gene','EFFECT','Drug','Corr','EFFECT2','MIC']].fillna('-').drop_duplicates().sort_values('POS')
        MIC_tab.to_csv('snpdr.tsv',sep='\t',index=False)
    else:
        open('snpdr.tsv','w').write(f'CHROM\tPOS\tREF\tALT\tEVIDENCE\tGene\tEFFECT\tDrug\tCorr\tEFFECT2\tMIC')
        open('Chin_snpdr.tsv','w').write(f'染色体\t变异位置\t参考碱基\t突变碱基\t变异频率\t基因名\t变异效应\t药物名\t相关性\t突变结果\tMIC')

def get_key(dn,valuen):
        return [k for k,v in dn.items() if v.split(':')[0] == valuen]

def drug_guide(tdgfile,tdrfile,ofn,Sample):
    INH_D,RIF_D,PZA_D,EMB_D,LFX_D,BDQ_D,CFZ_D,LZD_D,CYS_D,MXF_D=['敏感']*10
    dgfile = pd.read_table(tdgfile)
    drfile = pd.read_table(tdrfile)
    drfile = drfile.fillna('-')
    #drfile = drfile.loc[~(drfile['Corr']=="Not assoc w R")]
    drfile = drfile.loc[(drfile['Corr'].str.contains("Assoc w R"))]
    drlist = drfile['Drug'].apply(lambda x:x.split(':')[0]).unique().tolist()
    loc=locals()
    if drlist:
        for i in drlist:
            if len(get_key(drug_dict,i)) > 0:
                DEN = get_key(drug_dict,i)[0]
            else:
                DEN = i
            loc[f'{DEN}_D']='耐药'
            drfile.loc[drfile['Drug'].str.contains(i)][['Drug','EFFECT2','Corr','MIC']].rename(columns={'Drug':'药物名称','EFFECT':'突变效应','Corr':'关联性'}).to_csv(f'{i}_snpMIC.tsv',sep='\t',index=0,mode='w',header=False)
            if i in Firstline:
                open('firstline.tsv','a').write(f'{i}\t耐药\n')
            else:
                open('Secondline.tsv','a').write(f'{i}\t耐药\n')
        '''
        if len(drlist) < 2:
            drtype='DR'
        else:
            if '利福平' in drlist and '异烟肼' in drlist:
                if '莫西沙星' in drlist or '左氧氟沙星' in drlist:
                    if '阿米卡星' in drlist or '卡那霉素' in drlist or '卷曲霉素' in drlist:
                        drtype='XDR'
                    else:
                        drtype='Pre-XDR'
                else:
                    drtype='MDR'
            else:
                if '利福平' in drlist or '异烟肼' in drlist:
                    drtype='Pre-MDR'
                else:
                    drtype='DR'
        '''
        flq = ['氧氟沙星','左氧氟沙星','莫西沙星']
        gpa = ['贝达喹啉','利奈唑胺']
        if '利福平' in drlist and not '异烟肼' in drlist and not any([i for i in drlist if i in flq]):
            drtype = 'RR-TB'
        elif not '利福平' in drlist and '异烟肼' in drlist:
            drtype = 'HR-TB'
        elif ('利福平' in drlist and '异烟肼' in drlist) and not any([i for i in drlist if i in flq]):
            drtype = "MDR-TB"
        elif '利福平' in drlist and (any([i for i in drlist if i in flq]) and not any([i for i in drlist if i in gpa])):
            drtype = "Pre-XDR-TB"
        elif  '利福平' in drlist and (any([i for i in drlist if i in flq]) and any([i for i in drlist if i in gpa])):
            drtype = "XDR-TB"
        else:
            drtype = "Other"
    else:
        drtype='Sensitive'
    drline = '\t'.join(drlist)
    open('All_snpMIC.tsv','w').write(f'抗生素名称\t突变效应\t关联性\tMIC\n')
    subprocess.run('cat *snpMIC.tsv >> All_snpMIC.tsv',shell=True)
    if drlist:
        open('sam_dr.tsv','w').write(f'{Sample}\t{drtype}\t{drline}\n')
    else:
        open('sam_dr.tsv','w').write(f'{Sample}\t{drtype}\n')
    tmpdg = dgfile.columns.tolist()
    tmpdg = tmpdg[1:-2]
    for drug in tmpdg:
        if drug in drlist:
            dgfile = dgfile.loc[~dgfile[drug].str.contains("S")]
    dgfile.to_csv(ofn,sep='\t',index=False)
    if dgfile.shape[0] > 0:
        dgfile[['推荐方案','备注']].iloc[0,:].fillna('-').to_csv('DG.tsv',sep='\t',index=False,header=False)
    else:
        open('DG.tsv','w').write('没有适配方案\n')
    if '利福平' not in drlist:
           subprocess.run('cp /data/deploy/TB_soft/ref/TB/RifS_Drug_guide.tsv ./',shell=True)
    open('DG1.tsv','w').write('用药方案')
  
    DGmodel=f'''用药方案
抗生素名称\t是否敏感
（A组）未知情况（2HRZE/4HR）
异烟肼  {loc['INH_D']}
利福平  {loc['RIF_D']}
吡嗪酰胺    {loc['PZA_D']}
乙胺丁醇    {loc['EMB_D']}

（B组）异烟肼耐药（6-9RZELfx）
利福平  {loc['RIF_D']}
吡嗪酰胺    {loc['PZA_D']}
乙胺丁醇    {loc['EMB_D']}
左氧氟沙星  {loc['LFX_D']}
（C组）利福平耐药-氟喹诺酮类不耐药（6Lfx（Mfx）BdqLzd（Cs）Cfz /12Lfx（Mfx）Lzd（Cs）Cfz）
左氧氟沙星  {loc['LFX_D']}
贝达喹啉    {loc['BDQ_D']}
氯法齐明    {loc['CFZ_D']}
莫西沙星    {loc['MXF_D']}
利奈唑胺    {loc['LZD_D']}
（D组）氟喹诺酮类耐药（6 BdqLzd Cfz Cs/14 Lzd Cfz Cs）
贝达喹啉    {loc['BDQ_D']}
氯法齐明    {loc['CFZ_D']}
环丝氨酸    {loc['CYS_D']}
利奈唑胺    {loc['LZD_D']}
备注：
MIC:最低抑菌浓度、MDR（Multidrug-Resistant）:利福平、异烟肼耐药和另外两种耐药   XDR(Extensively Drug-resistant)：利福平、异烟肼耐药 包含氟喹诺酮类耐药以及二线注射药物耐药

'''            
    open('DG1.tsv','a').write(DGmodel)
    DG2f = pd.read_table('All_snpMIC.tsv')
    DG2f['tmpn'] =  DG2f['抗生素名称'].apply(lambda x : x.split(':')[0])
    DG2f.rename(columns={'抗生素名称':'抗生素'},inplace=True)
    FDG =  DG2f[DG2f['tmpn'].isin(Firstline)]
    SDG =  DG2f[~(DG2f['tmpn'].isin(Firstline))]
    FDG.drop('tmpn',inplace=True,axis=1)
    FDG['耐药性'] = '耐药'
    FDG = FDG[['抗生素','耐药性','突变效应','MIC']]
    SDG.drop('tmpn',inplace=True,axis=1)
    SDG['耐药性'] = '耐药'
    SDG = SDG[['抗生素','耐药性','突变效应','MIC']]
    FDG.drop_duplicates().to_csv('FirstData.tsv',sep='\t',index=False)
    SDG.drop_duplicates().to_csv('SecondData.tsv',sep='\t',index=False)
def rcdr(inf,ofn):
    three_letter ={'V':'VAL', 'I':'ILE', 'L':'LEU', 'E':'GLU', 'Q':'GLN', \
'D':'ASP', 'N':'ASN', 'H':'HIS', 'W':'TRP', 'F':'PHE', 'Y':'TYR',    \
'R':'ARG', 'K':'LYS', 'S':'SER', 'T':'THR', 'M':'MET', 'A':'ALA',    \
'G':'GLY', 'P':'PRO', 'C':'CYS'}
    gffdb = '/data/deploy/TB_soft/ref/TB/gene1.tsv'
    refl = '/data/deploy/TB_soft/ref/TB/ref.fa'
    gffdf = pd.read_table(gffdb)
    #---get info from gff----
    afile = pd.read_table(inf)
    poslist = afile[afile['TYPE']=='snp'][['CHROM','POS']].drop_duplicates()['POS'].tolist()
    tmpdic = {}
    n=1
    for i in range(len(poslist)):
        if i > 0:
            if poslist[i] - poslist[i-1] < 3:
                tmpdic[n] = [poslist[i-1],poslist[i]]
                n+=1
    for key,value in tmpdic.items():
        snppos1 = value[0]
        snppos2 = value[1]
        ref1 = afile[afile['POS']==snppos1]['REF'].tolist()[0]
        ref2 = afile[afile['POS']==snppos2]['REF'].tolist()[0]
        snp1 = afile[afile['POS']==snppos1]['ALT'].tolist()[0]
        snp2 = afile[afile['POS']==snppos2]['ALT'].tolist()[0]
        gene = afile[afile['POS']==snppos1]['GENE'].tolist()[0]
        if gene != 'nogene':
            try:
                strand = gffdf[gffdf['gene']==gene]['Strand'].tolist()[0]
            except:
                gene = gffdf[(gffdf['GStart']<snppos1) & (gffdf['GEnd']>snppos1)]['gene'].tolist()[0]               
                
            strand = gffdf[gffdf['gene']==gene]['Strand'].tolist()[0]
            geneS = gffdf[gffdf['gene']==gene]['GStart'].tolist()[0]
            geneE = gffdf[gffdf['gene']==gene]['GEnd'].tolist()[0]
            Chr = gffdf[gffdf['gene']==gene]['Contig'].tolist()[0]
            if strand == '+':
                GeneP1 = snppos1-geneS+1
                GeneP2 = snppos2-geneS+1
            else:
                GeneP1 = geneE - snppos1 + 1
                GeneP2 = geneE - snppos2 + 1
                
            if math.ceil((GeneP1/3)) == math.ceil((GeneP2/3)):
                coden = math.ceil((GeneP1/3))
                codenS = geneS+(coden-1)*3+1
                codenE = geneS+(coden-1)*3+3
                codenP1 = GeneP1%3-1
                codenP2 = GeneP2%3-1
                open('tmp.bed','w').write(f'{Chr}\t{codenS-2}\t{codenE-1}\n')
                subprocess.run(f'bedtools getfasta -fi {refl} -bed tmp.bed > tmp.fa',shell=True)
                with open('tmp.fa') as f1:
                    for line in f1:
                        if not line.startswith('>'):
                            Codenl = list(line.strip())
                            CodenRef = Seq(''.join(Codenl))
                            Codenl[codenP1] = snp1
                            Codenl[codenP2] = snp2
                            CodenAlt = Seq(''.join(Codenl))
                            if CodenRef.translate() == '*':
                                ProRef = 'fs'
                            else:
                                ProRef = three_letter.get(str(CodenRef.translate())).capitalize() 
                            if three_letter.get(str(CodenAlt.translate())):
                                ProALT = three_letter.get(str(CodenAlt.translate())).capitalize()
                            else:
                                ProALT = 'fs'
                afile['EFFECT'] = afile.apply(lambda x : x['EFFECT'] if x['POS']!=snppos1 else f'''{x['EFFECT'].split(' ')[0]} c.{GeneP1}{ref1}>{snp1} p.{ProRef}{coden}{ProALT}''',axis=1)
                afile['EFFECT'] = afile.apply(lambda x : x['EFFECT'] if x['POS']!=snppos2 else f'''{x['EFFECT'].split(' ')[0]} c.{GeneP2}{ref2}>{snp2} p.{ProRef}{coden}{ProALT}''',axis=1)
                #afile.iloc[afile['POS']==snppos2,'EFFECT'] = f'{gene};{ProRef}{coden}{ProALT}'
    afile.rename(columns = {'GENE':'Gene','Drug_x':'Drug'},inplace=True)
    afile.to_csv(ofn,sep='\t',index=False)
#----转换在启动子区域的突变为'-' 例如-15C>T-----
genedb = pd.read_table('/data/deploy/TB_soft/ref/TB/gene1.tsv')
Drgene = ['ahpC','inhA','katG','mshA','ndh','Rv1258c','Rv2752c','rpoB','rpoA','rpoC','embA','embB','embC','embR','ubiA','pncA','clpC1','panD','PPE35','Rv3236c','gyrA','gyrB','pepQ','Rv0678','mmpL5','mmpS5','atpE','Rv1979c','rplC','rrl','fgd1','ddn','fbiA','fbiB','fbiC','Rv2983','rrs','eis','whiB7','rpsL','gid','inhA','ethA','tlyA','whiB6','ccsA','fprA','aftB','ethR','mshA','Rv3083','ndh','ccsA','fprA','fabG1']
def transMut(pos,ref,alt,db=genedb,exon=1000):
    tmpdb = db[(db['GStart']<pos+1000) & (db['GEnd']>pos-1000)]
    tmplist = []
    for i in tmpdb['gene'].tolist():
        #i = i.replace('gene:','')
        std = tmpdb.loc[tmpdb['gene']==i,'Strand'].tolist()[0]
        GS = tmpdb.loc[tmpdb['gene']==i,'GStart'].tolist()[0]
        GE = tmpdb.loc[tmpdb['gene']==i,'GEnd'].tolist()[0]
        i = i.replace('gene:','')
        if i in Drgene:
            if std == '+':
                if pos < GS:
                    npos=GS-pos
                    tmplist.append(f'intergenic_region c.-{npos}{ref.upper()}>{alt.upper()}')
            else:
                if pos > GE:
                    npos=pos - GE
                    tmplist.append(f'intergenic_region c.-{npos}{ref.upper()}>{alt.upper()}')
    if not tmplist:
        tmplist = [f'intergenic_region c.-{pos}{ref.upper()}>{alt.upper()}']
    return tmplist[0]

    
def fq_process(gs,sam,intype,Ref,nt,dep,fwkdir,fq1,fq2=0):
    os.chdir(fwkdir)
    make_dir(f'{fwkdir}/{sam}')
    if gs==2:
        bwacode=f'''bwa mem -v 1  -Y -M -R '@RG\\tID:{sam}\\tSM:{sam}' -t {nt} {Ref} {fq1} {fq2} | samtools sort -n -l 0 -T /tmp --threads 1 -m 1600M | samtools fixmate -m --threads {nt} - - | samtools sort -l 0 -T /tmp --threads 1 -m 1600M | samtools markdup -T /tmp --threads {nt} --verbosity 1 -r - - > snps.bam '''
        bwacodes=f'''bwa mem -v 1 -Y -M -R '@RG\\tID:{sam}\\tSM:{sam}' -t {nt} {Ref} {fq1}| samtools sort -n -l 0 -T /tmp --threads {nt} -m 1600M | samtools fixmate -m --threads {nt} - - | samtools sort -l 0 -T /tmp --threads 1 -m 1600M | samtools markdup -T /tmp --threads {nt} --verbosity 1 -r - - > snps.bam'''
    elif gs==3:
        bwacodes=f'''minimap2 -ax map-ont -R '@RG\\tID:{sam}\\tSM:{sam}' -t {nt} {Ref} {fq1}|samtools sort -n -l 0 -T /tmp --threads 1 -m 1600M | samtools fixmate -m --threads {nt} - - | samtools sort -l 0 -T /tmp --threads 1 -m 1600M | samtools markdup -T /tmp --threads {nt} -r -s - - > snps.bam'''
    csnpcode=f'''{snippydir}/binaries/noarch/freebayes-parallel {wkdir}/ref/ref.txt {nt} -p 2 -P 0 -C 2 -F 0.05 --min-coverage {dep} --min-repeat-entropy 1.0 -q 30 -m 30 --strict-vcf -f {Ref} snps.bam > snps.raw.vcf'''
    #-----质控1 q30 Q30---------
    if gs==2:
        fil1snp=f''' bcftools view --include 'QUAL>=20 && FMT/DP>={dep} && (FMT/AO)/(FMT/DP)>=0.1' snps.raw.vcf  | {snippydir}/binaries/linux/vt normalize -r {Ref} - | bcftools annotate --remove '^INFO/TYPE,^INFO/DP,^INFO/RO,^INFO/AO,^INFO/AD,^INFO/AB,^FORMAT/GT,^FORMAT/DP,^FORMAT/AD,^FORMAT/RO,^FORMAT/AO,^FORMAT/QR,^FORMAT/QA,^FORMAT/GL' > snps.filt1.vcf'''
    elif gs==3:
        fil1snp=f''' bcftools view --include 'QUAL>=10 && FMT/DP>={dep} && (FMT/AO)/(FMT/DP)>=0.1' snps.raw.vcf  | {snippydir}/binaries/linux/vt normalize -r {Ref} - | bcftools annotate --remove '^INFO/TYPE,^INFO/DP,^INFO/RO,^INFO/AO,^INFO/AD,^INFO/AB,^FORMAT/GT,^FORMAT/DP,^FORMAT/AD,^FORMAT/RO,^FORMAT/AO,^FORMAT/QR,^FORMAT/QA,^FORMAT/GL' > snps.filt1.vcf'''

    fil2snp=f''' bcftools view --include 'FMT/GT="1/1" && QUAL>=100 && FMT/DP>={dep} && (FMT/AO)/(FMT/DP)>=0.75' snps.raw.vcf  | {snippydir}/binaries/linux/vt normalize -r {Ref} - | bcftools annotate --remove '^INFO/TYPE,^INFO/DP,^INFO/AD,^INFO/RO,^INFO/AO,^INFO/AB,^FORMAT/GT,^FORMAT/DP,^FORMAT/AD,^FORMAT/RO,^FORMAT/AO,^FORMAT/QR,^FORMAT/QA,^FORMAT/GL' > snps.filt.vcf'''    
    annvcf=f'''grep -v '#' snps.filt.vcf |awk -v OFS='\t' '{{print $1,$2,$2+1,$3,$4,$5}}' > vcf.bed
    echo "染色体\t起始位点\t终止位点\tnt\t参考碱基\t突变碱基\t基因染色体\t基因起始\t基因终止\t基因所在链\t注释1\t机制\t基因名称\t机制2\t机制模型\t结果\tRv_num" > {wkdir}/5.virulence_genes/tmp_{sam}_vcf.anno.bed
    bedtools intersect -a vcf.bed -b /data/deploy/TB_soft/ref/TB/VB_anno.bed -wa -wb > {wkdir}/5.virulence_genes/tmp1_{sam}_vcf.anno.bed
    cat  {wkdir}/5.virulence_genes/tmp1_{sam}_vcf.anno.bed >> {wkdir}/5.virulence_genes/tmp_{sam}_vcf.anno.bed
    '''
    #-----质控2 genotype Q > 100 && 变异类型要纯合突变，深度>=5
    fildr1snp=f'''vcftools --vcf snps.filt.vcf --exclude-bed {black_list} --recode --recode-INFO-all --out filtdr'''
    #-----排除耐药基因位点---
    annosnp=f'''/home/dell/biosoft/snpEff/scripts/snpEff ann -noLog -noStats -upDownStreamLen 2000 -spliceSiteSize 0 -spliceRegionExonSize 0 -spliceRegionIntronMax 0 -spliceRegionIntronMin 0 -c {wkdir}/ref/reference/snpeff.config -dataDir . ref snps.subs.vcf > snps.vcf'''
    anno2snp=f'''/home/dell/biosoft/snpEff/scripts/snpEff ann -noLog -noStats -upDownStreamLen 2000 -spliceSiteSize 0 -spliceRegionExonSize 0 -spliceRegionIntronMax 0 -spliceRegionIntronMin 0 -c {wkdir}/ref/reference/snpeff.config -dataDir . ref filtdr.recode.vcf > snps.filtdr.vcf'''
    anno1snp=f'''/home/dell/biosoft/snpEff/scripts/snpEff ann -noLog -noStats -upDownStreamLen 2000 -spliceSiteSize 0 -spliceRegionExonSize 0 -spliceRegionIntronMax 0 -spliceRegionIntronMin 0 -c {wkdir}/ref/reference/snpeff.config -dataDir . ref snps.subs1.vcf > snps1.vcf'''
    #------SNP注释-----
    #cvtab1=f'''{snippydir}/bin/snippy-vcf_to_tab --gff {wkdir}/ref/reference/ref/genes.gff --ref {Ref} --vcf snps.filtdr.vcf > snps.filtdr.tab'''
    #cvtab2=f'''{snippydir}/bin/snippy-vcf_to_tab --gff {wkdir}/ref/reference/ref/genes.gff --ref {Ref} --vcf snps.vcf > snps.tab'''
    cvtab1 = f'''/data/deploy/TB_soft/other_soft/vcf2tab.py snps.filtdr.vcf snps.filtdr.tab'''
    cvtab2 = f'''/data/deploy/TB_soft/other_soft/vcf2tab.py snps.vcf snps.tab'''
    cvtab3 = f'''/data/deploy/TB_soft/other_soft/vcf2tab.py snps1.vcf snps1.tab'''
    #----整理成TAB文件---用于后续分析和统计
    subvcf=f'''{snippydir}/binaries/linux/vcfallelicprimitives -kg snps.filt.vcf > snps.subs.vcf'''
    if ck:
        sub1vcf=f'''cp tt.vcf snps.subs1.vcf'''
        sub2vcf=f'''cp filtdr.recode.vcf snps.drsubs.vcf'''
    else:
        sub1vcf=f'''{snippydir}/binaries/linux/vcfallelicprimitives -kg snps.filt1.vcf > snps.subs1.vcf'''
        sub2vcf=f'''{snippydir}/binaries/linux/vcfallelicprimitives -kg filtdr.recode.vcf > snps.drsubs.vcf'''
    #----将MNP拆散成SNP 去掉INS 和 DEL只保留SNP
    confa=f'''bcftools convert -Oz -o snps.vcf.gz snps.vcf
    bcftools index -f snps.vcf.gz
    awk -v OFS='\t' '{{if(int($4)<{dep})print $1,$2+1,$3+1}}' Sample_depth.regions.bed > mask.bed
    bcftools consensus -f {Ref} -o snps.consensus.fa snps.vcf.gz --mask mask.bed 
    bcftools convert -Oz -o snps.subs.vcf.gz snps.subs.vcf
    bcftools index -f snps.subs.vcf.gz
    bcftools consensus -f {Ref} -o snps.subs.consensus.fa snps.subs.vcf.gz --mask mask.bed
    '''
    con2fa=f'''bcftools convert -Oz -o snps.filtdr.vcf.gz snps.filtdr.vcf
    bcftools index -f snps.filtdr.vcf.gz
    awk -v OFS='\t' '{{if(int($4)<{dep})print $1,$2+1,$3+1}}' Sample_depth_1.regions.bed > mask.bed
    bcftools consensus -f {Ref} -o snps.filtdr.consensus.fa snps.filtdr.vcf.gz --mask mask.bed
    bcftools convert -Oz -o snps.drsubs.vcf.gz snps.drsubs.vcf
    bcftools index -f  snps.drsubs.vcf.gz
    bcftools consensus -f {Ref} -o snps.subs.drsubs.consensus.fa snps.drsubs.vcf.gz --mask mask.bed
    '''
    covcal=f'''mosdepth Sample_depth_1 -t {nt} -b 1 snps.bam
    gunzip -f Sample_depth_1.regions.bed.gz
    mosdepth Sample_depth -t {nt} -b 7500 snps.bam
    gunzip -f Sample_depth.regions.bed.gz
    '''
    #wkdirI="/".join(wkdir.split('/')[3:])
    #igvhtml=f'''python /data/deploy/TB_soft/other_soft/IGV_js/IGV_new.py -r /outputs/{wkdirI}/ref/ref.fa -m  /outputs/{wkdirI}/1.snp_calling/{sam}/snps.bam -o {wkdir}/1.snp_calling/{sam}/'''
    igvhtml=f'''python /data/deploy/TB_soft/other_soft/IGV_js/IGV_new.py -r ../../ref/ref.fa -m snps.bam -o {wkdir}/1.snp_calling/{sam}/ -s {sam}
cp /data/deploy/TB_soft/other_soft/IGV_js/igv.min.js {wkdir}/1.snp_calling/{sam}/
                '''
    cwigvhtml=f'''cp /home/dell/soft/clockwork/Ref.QC_and_map/ref.fa ./ref.fa
    samtools faidx ./ref.fa
    python /data/deploy/TB_soft/other_soft/IGV_js/IGV_new.py -r ../../ref/ref.fa -m snps.bam -o {wkdir}/1.snp_calling/{sam}/ -s {sam}
     cp /data/deploy/TB_soft/other_soft/IGV_js/igv.min.js {wkdir}/1.snp_calling/{sam}/
     sed -i 's/Chromosome/NC_000962.3/g' index_{sam}.html
     '''
    clockcmd1=f'''
        nextflow run /home/dell/soft/clockwork/nextflow/remove_contam.nf -with-singularity /home/dell/soft/clockwork/clockwork_container.img --ref_fasta /home/dell/soft/clockwork/Ref.remove_contam/ref.fa --ref_metadata_tsv /home/dell/soft/clockwork/Ref.remove_contam/remove_contam_metadata.tsv --reads_in1 tmp1.fastq --reads_in2  tmp2.fastq --outprefix remove_contam --mapping_threads 10'''
    clockcmd2=f'''nextflow run /home/dell/soft/clockwork/nextflow/variant_call.nf -with-singularity /home/dell/soft/clockwork/clockwork_container.img  --ref_dir /home/dell/soft/clockwork/Ref.QC_and_map/ --reads_in1 remove_contam.remove_contam.1.fq.gz --reads_in2 remove_contam.remove_contam.2.fq.gz  --output_dir variant_call_out --sample_name {sam}
    mv variant_call_out/minos/final.vcf {wkdir}/1.snp_calling/{sam}/snps.raw.vcf
    mv variant_call_out/samtools/rmdup.bam {wkdir}/1.snp_calling/{sam}/snps.bam
    samtools index {wkdir}/1.snp_calling/{sam}/snps.bam 
    rm -r /home/dell/soft/clockwork/work
    '''
    clockcmd3=f'''/data/deploy/TB_soft/other_soft/vcf2tab1.py snps.filtdr.vcf snps.filtdr.tab'''
    clockcmd4=f'''/data/deploy/TB_soft/other_soft/vcf2tab1.py snps.vcf snps.tab'''
    clockcmd5=f'''/data/deploy/TB_soft/other_soft/vcf2tab1.py snps1.vcf snps1.tab'''
    #-----------------生成一致性序列--------------------
    with open('mapping.task.txt','w') as f:
        if intype == 'fastq':
            if not ck:
                if fq2:
                    subprocess.run(bwacode,shell=True,stdout=f,stderr=f)
                else:
                    subprocess.run(bwacodes,shell=True)
                subprocess.run('samtools index snps.bam',shell=True,stdout=f,stderr=f)
                subprocess.run(csnpcode,shell=True,stdout=f,stderr=f)
                subprocess.run(fil1snp,shell=True,stdout=f,stderr=f)
                subprocess.run(fil2snp,shell=True,stdout=f,stderr=f)
                subprocess.run(fildr1snp,shell=True,stdout=f,stderr=f)
            else:
                subprocess.run(f'gunzip {fq1} -c > /home/dell/soft/clockwork/tmp1.fastq',shell=True,stdout=f,stderr=f)
                subprocess.run(f'gunzip {fq2} -c > /home/dell/soft/clockwork/tmp2.fastq',shell=True,stdout=f,stderr=f)
                os.chdir('/home/dell/soft/clockwork')
                subprocess.run(clockcmd1,shell=True,stdout=f,stderr=f)
                subprocess.run(clockcmd2,shell=True,stdout=f,stderr=f)
                os.chdir(f'{wkdir}/1.snp_calling/{sam}')
                subprocess.run(f''' sed -i 's/NC_000962.3/Chromosome/g' snps.raw.vcf  ''',shell=True,stdout=f,stderr=f)
                with open('snps.raw.vcf') as vcff:
                    for line in vcff:
                        if line.startswith('#'):
                            open('tt.vcf','a').write(line)
                        else:
                            if len(line.split('\t')[3]) <=1000 and len(line.split('\t')[4])<=1000:
                                open('tt.vcf','a').write(line)

            subprocess.run(f'cp tt.vcf ./snps.filt.vcf',shell=True,stdout=f,stderr=f)
            subprocess.run(f'cp tt.vcf ./snps.filt1.vcf',shell=True,stdout=f,stderr=f)
            subprocess.run(f'vcftools --vcf snps.filt.vcf --exclude-bed {black_list} --recode --recode-INFO-all --out filtdr',shell=True,stdout=f,stderr=f)
        elif intype == 'bam':
            modify_bam_reference(fq1, 'snps.bam', 'Chromosome')
            subprocess.run('samtools index snps.bam',shell=True,stdout=f,stderr=f)
            subprocess.run(csnpcode,shell=True,stdout=f,stderr=f)
            subprocess.run(fil1snp,shell=True,stdout=f,stderr=f)
            subprocess.run(fil2snp,shell=True,stdout=f,stderr=f)
            subprocess.run(fildr1snp,shell=True,stdout=f,stderr=f)
        else:
            modify_vcf_reference(fq1, 'snps.raw.vcf', 'Chromosome')
            subprocess.run(fil1snp,shell=True,stdout=f,stderr=f)
            subprocess.run(fil2snp,shell=True,stdout=f,stderr=f)
            subprocess.run(fildr1snp,shell=True,stdout=f,stderr=f)

        subprocess.run(subvcf,shell=True,stdout=f,stderr=f)
        subprocess.run(sub1vcf,shell=True,stdout=f,stderr=f)
        subprocess.run(sub2vcf,shell=True,stdout=f,stderr=f)
        subprocess.run(annvcf,shell=True,stdout=f,stderr=f)
        subprocess.run(annosnp,shell=True,stdout=f,stderr=f)
        subprocess.run(anno1snp,shell=True,stdout=f,stderr=f)
        subprocess.run(anno2snp,shell=True,stdout=f,stderr=f)
        sys.stdout.flush()
        if not ck:
            subprocess.run(cvtab1,shell=True,stdout=f,stderr=f)
            subprocess.run(cvtab2,shell=True,stdout=f,stderr=f)
            subprocess.run(cvtab3,shell=True,stdout=f,stderr=f)
            subprocess.run(igvhtml,shell=True,stdout=f,stderr=f)
        else:
            subprocess.run(clockcmd3,shell=True,stdout=f,stderr=f)
            subprocess.run(clockcmd4,shell=True,stdout=f,stderr=f)
            subprocess.run(clockcmd5,shell=True,stdout=f,stderr=f)
            subprocess.run(cwigvhtml,shell=True,stdout=f,stderr=f)
        subprocess.run(covcal,shell=True,stdout=f,stderr=f)
        subprocess.run(confa,shell=True,stdout=f,stderr=f)
        subprocess.run(con2fa,shell=True,stdout=f,stderr=f)
        #subprocess.run(igvhtml,shell=True,stdout=f,stderr=f)
    #tmpfile = pd.read_table(f'{wkdir}/5.virulence_genes/tmp_{sam}_vcf.anno.bed')
    #tmpfile = tmpfile.drop(['nt','染色体','基因染色体','机制模型','注释1'],axis=1)
    #tmpfile['SN'] = sam
    #tmpfile.to_csv(f'{wkdir}/5.virulence_genes/{sam}_vcf.anno.bed',sep='\t',index=0)
    #subprocess.run(f'rm {wkdir}/5.virulence_genes/tmp_{sam}_vcf.anno.bed',shell=True)
    subprocess.run(f'''sed -i 's/gene:EBG00000313325/rrs/g' snps1.vcf ''',shell=True)
    subprocess.run(f'''sed -i 's/gene:EBG00000313339/rrl/g' snps1.vcf ''',shell=True)
    subprocess.run(f'''sed -i 's/gene:EBG00000313349/rrf/g' snps1.vcf ''',shell=True)
    subprocess.run(f'''sed -i 's/gene:rrf/rrf/g' snps1.tab ''',shell=True)
    subprocess.run(f'''sed -i 's/gene:rrl/rrl/g' snps1.tab ''',shell=True)
    subprocess.run(f'''sed -i 's/gene:rrs/rrs/g' snps1.tab ''',shell=True)
    #20240301 新鉴定
    #rcdr('snps1.tab','snps2.tab')
    #rcdrdb = pd.read_table('snps2.tab')
    #if rcdrdb.shape[0] >0:
    #    rcdrdb.loc[rcdrdb['Gene']=='nogene','EFFECT'] =  rcdrdb[rcdrdb['Gene']=='nogene'].apply(lambda x: transMut(x['POS'],x['REF'],x['ALT']),axis=1)
    #rcdrdb.to_csv('snps3.tab',sep='\t',index=False)
    Drugvcf('snps1.vcf')
    #dr_anno(drtab,'snps3.tab','snpdr.tsv')
    drug_guide(DrugGfile,'snpdr.tsv','DrugGuide.tsv',sam)
    open(f'map_OK','w').write(f'OK')
    
    #transdf = pd.read_table('snpdr.tsv')
    #transdf['Drug'] = transdf['Drug'].str.split(':').str.join(' ')
#------根据snp数量进行成簇分析--------
#-----搜索字典------
def nindic(a,b):
    for value in b.values():
        if a in value:
            return True
            break
    else:
        return False
def findval(a,b):
    for key,value in b.items():
        if a in value:
            return key
            break
def findval2(a,b):
    for key in b.keys():
        if a in b[key]['sample']:
            return key
            break
 
def convcf(filen,pre):
    a = pd.read_table(filen)
    a = a.iloc[:,9:]
    df2 = pd.DataFrame(a.values.T, index=a.columns, columns=a.index)
    df2.replace(1,2,inplace=True)
    df2.replace(0,1,inplace=True)
    df2.to_csv(f'{pre}.vcf',sep='\t',header=True)
           
def trans_cluster(vcf,cutoff,outdir):
    file1 = pd.read_table(vcf)
    Clus=0
    All_cluster = {}
    Cluster_sum = {}
    outsam = []
    samplenames = file1.columns[9:]
    samnum = len(samplenames) 
    tof=f'Transmission_Cluster_cutoff_{cutoff}.txt'
    samplenames = file1.columns[9:]
    samnum = len(samplenames) 
    for i in range(samnum-1):
        n=i+1
        while n < samnum:
            samplename1 = samplenames[i]
            samplename2 = samplenames[n]
            diff_num = sum(map(abs,file1.iloc[:,i+9]-file1.iloc[:,n+9]))
            open('compare_result.txt','a').write(f'{samplename1}\t{samplename2}\t{diff_num}\n')
            if diff_num <=int(cutoff):
                if not nindic(samplename1,All_cluster):
                    if samplename1 not in All_cluster.keys():
                        All_cluster[samplename1] = [samplename2]
                        Clus+=1
                        if Clus not in Cluster_sum.keys():
                            Cluster_sum[Clus]= {}
                            Cluster_sum[Clus]['sample']=[samplename1,samplename2]
                        else:
                            Cluster_sum[Clus]['sample'].append(samplename1)
                            Cluster_sum[Clus]['sample'].append(samplename2)
                            Cluster_sum[Clus]['sample'] = list(set(Cluster_sum[Clus]['sample']))

                    else:
                        All_cluster[samplename1].append(samplename2)
                        Cluster_sum[Clus]['sample'].append(samplename2)
                        #f.write(f'Cluster_{Clus}\t{samplename1}\t{samplename2}\t{diff_num}\n')
                else:
                    All_cluster[findval(samplename1,All_cluster)].append(samplename2)
                    All_cluster[findval(samplename1,All_cluster)] = list(set(All_cluster[findval(samplename1,All_cluster)]))
                    OClus = findval2(samplename1,Cluster_sum)
            else:
                outsam.append(samplename1)
                outsam.append(samplename2)
                outsam = list(set(outsam))
                
            n+=1
    def xx(a,b):
        for key,value in b.items():
            if a in value:
                return 1
                break
            else:
                return 0 

    #new_clus=Cluster_sum.copy()
    for sam in outsam:
        if xx(sam,Cluster_sum):
            pass
        else:
            Clus+=1
            Cluster_sum[Clus]={}
            Cluster_sum[Clus]['sample'] = [sam]
    
    with open(tof,'w') as f:
        f.write(f'Cluster\tsample1\tsample2\tdiffsnps\n')
        for key in Cluster_sum.keys():  
            for i in range(len(Cluster_sum[key]['sample'])-1):
                n=i+1
                while n < len(Cluster_sum[key]['sample']):
                    sample1=Cluster_sum[key]['sample'][i]
                    sample2=Cluster_sum[key]['sample'][n]
                    Diffn = sum(map(abs,file1[sample1]-file1[sample2]))
                    f.write(f'Cluster_{key}\t{sample1}\t{sample2}\t{Diffn}\n')
                    n+=1
    #----统计结果------
    try:
        tof2 = f'Summary_TransCluster_cutoff{cutoff}.txt'
        Clu_file = pd.read_table(tof)
        Clu_sum = Clu_file.pivot_table(values='diffsnps',index='Cluster',aggfunc=['max','min','mean'])
        pd.DataFrame(Clu_sum).to_csv(f'{outdir}/{tof2}',sep='\t',header=True)
        with open(f'cutoff{cutoff}_sample_num_summary','w') as f:
            for key in Cluster_sum.keys():
                Sample_num = len(Cluster_sum[key]['sample'])
                Sample_aa = " ".join(Cluster_sum[key]['sample'])
                f.write(f'Cluster_{key}\tSample_num:\t{Sample_num}\tSample:{Sample_aa}\n')
    except:
        print('所有样本均未成簇')
                          
#-------kraken2-------
def kk2(nT,Pre,fq1,fq2=0):
    with open('kk2.log','w') as kk2l:
        if fq2:
            subprocess.run(f'kraken2 --db {Krdb} --threads {nT} --output {wkdir}/4.kraken_taxonomic/{Pre}.txt --report {wkdir}/4.kraken_taxonomic/{Pre}.kraken2.txt {fq1} {fq2}',shell=True,stdout=kk2l,stderr=kk2l)
        else:
            subprocess.run(f'kraken2 --db {Krdb} --threads {nT} --output {wkdir}/4.kraken_taxonomic/{Pre}.txt --report {wkdir}/4.kraken_taxonomic/{Pre}.kraken2.txt {fq1}',shell=True,stdout=kk2l,stderr=kk2l)
        subprocess.run(f'bracken -d {Krdb} -o {wkdir}/4.kraken_taxonomic/{Pre}.bracken1.txt -w {wkdir}/4.kraken_taxonomic/{Pre}.bracken2.txt -l S -t 10  -i {wkdir}/4.kraken_taxonomic/{Pre}.kraken2.txt',shell=True,stdout=kk2l,stderr=kk2l)
        tmpfile = pd.read_table(f'{wkdir}/4.kraken_taxonomic/{Pre}.bracken1.txt')
        tmpfile = tmpfile[['name','taxonomy_id','taxonomy_lvl','new_est_reads','fraction_total_reads']]
        tmpfile.rename(columns={'name':'物种','taxonomy_id':'taxid','taxonomy_lvl':'水平','new_est_reads':'序列数量','fraction_total_reads':'相对丰度'},inplace=True)
        sanofile = pd.read_table('/data/deploy/TB_soft/other_soft/taxa_info_0717.txt')
        tbfile = pd.read_table('/data/deploy/TB_soft/ref/TB/TB_TAXA_info.tsv')
        tbfile['taxid'] = tbfile['taxid'].apply(lambda x:str(x))
        antmpfile = pd.merge(tmpfile,sanofile,on='taxid')
        antmpfile.fillna('-',inplace=True)
        antmpfile=antmpfile[['物种','taxid','序列数量','相对丰度','中文名','致病源性','可能引起的疾病','基因组长度']]
        antmpfile.to_csv(f'{wkdir}/4.kraken_taxonomic/{Pre}.tb.report1.txt',index=False,sep='\t')
        antmpfile['taxid'] = antmpfile['taxid'].apply(lambda x:str(x))
        tbfile2 = pd.merge(tbfile,antmpfile,on='taxid')
        tbfile2 = tbfile2[['中文名_x','序列数量','相对丰度','taxid']]
        tbfile2.rename(columns={'中文名_x':'物种','taxid':'NCBI物种号'},inplace=True)
        tbfile2.sort_values('序列数量',ascending=False,inplace=True)
        tbfile2 = tbfile2.iloc[:5,:]
        tbfile2['相对丰度'] = tbfile2['相对丰度']*100
        tbfile2['相对丰度'] =  tbfile2['相对丰度'].round(2)
        tbfile2['相对丰度'] = tbfile2['相对丰度'].apply(lambda x:'<1%' if x < 1 else x)
        tbfile2.to_csv(f'{wkdir}/4.kraken_taxonomic/{Pre}.tb.report2.txt',index=False,sep='\t')
        subprocess.run(f'python {sc1} -r {wkdir}/4.kraken_taxonomic/{Pre}.bracken2.txt -o {wkdir}/4.kraken_taxonomic/{Pre}.krona_input.txt',shell=True)
        subprocess.run(f'perl /data/deploy/TB_soft/other_soft/KronaTools-2.8/bin/ktImportText -o {wkdir}/4.kraken_taxonomic/{Pre}.krona.html {wkdir}/4.kraken_taxonomic/{Pre}.krona_input.txt',shell=True)
        kran_dic = {}
        def kran_summ(ofn1,ofn2,kfile):
            with open(ofn1,'w') as f1:
                with open(ofn2,'w') as f2:
                    with open(kfile) as f3:
                        f1.write(f'物种\t序列数量\t比例\tNCBI物种号\t中文名\t致病源分类\t可>能引起的疾病\t分类\t属\t科\t基因组长度\n')
                        f2.write(f'亚种/菌株\t序列数量\t比例\tNCBI物种号\t中文名\t致病源分类\t可能引起的疾病\t分类\t种\t属\t基因组长度\n')
                        for line in f3.readlines():
                            line = line.strip().split('\t')
                            line[5] = line[5].strip()
                            if line[5] == 'unclassified':
                                #-----[物种英文，匹配到的reads数量,百分比]
                                f1.write(f'{line[5]}\t{line[1]}\t{line[0]}\t-\t-\t-\t-\t-\t-\t-\t-\n')
                                f2.write(f'{line[5]}\t{line[1]}\t{line[0]}\t-\t-\t-\t-\t-\t-\t-\t-\n')
                            else:
                                if line[3] == 'D': kran_dic['taD']=line[5]  
                                elif line[3] =='K': kran_dic['taK']=line[5]
                                elif line[3] =='P': kran_dic['taP']=line[5]
                                elif line[3] =='C': kran_dic['taC']=line[5]
                                elif line[3] =='O': kran_dic['taO']=line[5]
                                elif line[3] =='F': kran_dic['taF']=line[5]
                                elif line[3] =='G': kran_dic['taG']=line[5]
                                elif line[3] =='S' and int(line[1]) > 2:
                                    Chin = chin.get(line[4],'-')
                                    Info = info.get(line[4],f'-\t-')
                                    Lent = lth.get(line[4],'-')
                                    kran_dic['taS']=line[5]
                                    taS=kran_dic.get('taS','-')
                                    taD=kran_dic.get('taD','-')
                                    taF=kran_dic.get('taF','-')
                                    taG=kran_dic.get('taG','-')
                                    Line_new = f'{taS}\t{line[1]}\t{line[0]}\t{line[4]}\t{Chin}\t{Info}\t{taD}\t{taF}\t{taG}\t{Lent}\n'
                                    f1.write(Line_new)
                                elif line[3] in ['S1','S2','S3'] and  int(line[1]) > 2:
                                    Line_new = f'{line[-1]}\t{line[1]}\t{line[0]}\t{line[4]}\t{Chin}\t{Info}\t{taD}\t{taF}\t{taG}\t{Lent}\n'
                                    f2.write(Line_new)
        kran_summ(f'{wkdir}/4.kraken_taxonomic/{Pre}.list.txt',f'{wkdir}/4.kraken_taxonomic/{Pre}.list2.txt',f'{wkdir}/4.kraken_taxonomic/{Pre}.kraken2.txt')
        tafile = pd.read_table(f'{wkdir}/4.kraken_taxonomic/{Pre}.list2.txt')
        tbfile = pd.read_table('/data/deploy/TB_soft/ref/TB/TB_TAXA_info.tsv')
        tcfile = pd.merge(tafile,tbfile,left_on='亚种/菌株',right_on = 'taxname',how='left')[['中文名_y','NCBI物种号','序列数量','比例']]
        tcfile.fillna('-',inplace=True)
        ta1file = pd.read_table(f'{wkdir}/4.kraken_taxonomic/{Pre}.list.txt')
        tdfile = pd.merge(ta1file,tbfile,left_on='物种',right_on = 'taxname',how='left')[['中文名_y','NCBI物种号','序列数量','比例']]
        tdfile.fillna('-',inplace=True)
        ffile =  pd.concat([tdfile,tcfile]).sort_values('序列数量',ascending=False).rename(columns={'中文名_y':'物种'})[['物种','序列数量','比例','NCBI物种号']]
        ffile = ffile[ffile['比例']>0.01]
        mains = ffile.iloc[0,0]
        open(f'{wkdir}/MainSpe.txt','a').write(f'{Pre}\t{mains}\n')
        if ffile.shape[0]>5:
            ffile = ffile.iloc[:5,:]
        ffile.to_csv(f'{wkdir}/4.kraken_taxonomic/{Pre}_TB_prop.list.tsv',sep='\t',index=0)
        subprocess.run(f'mv {wkdir}/4.kraken_taxonomic/{Pre}.list2.txt {wkdir}/4.kraken_taxonomic/t_{Pre}.list2.txt',shell=True)
        subprocess.run(f'mv {wkdir}/4.kraken_taxonomic/{Pre}_TB_prop.list.tsv {wkdir}/4.kraken_taxonomic/{Pre}.list2.txt',shell=True)
        open(f'{wkdir}/4.kraken_taxonomic/{Pre}_ok','w').write('')

#-------------summary----------------
#---fq----
def summ_fq(indir,outdir):
    fplsit = [i for i in os.listdir(indir) if 'json' in i ]
    mappat = re.compile('(\S+)\((\S+)\)')
    fp_dic = {}
    fp_dic2 = {}
    open(f'{outdir}/result_base.txt','w').write(f'样本名称\t总碱基数\t总序列数\tgc含量\t平均长度\tQ20\tQ30\t过滤后数量\n')
    for i in fplsit: 
        samname=i.split('.')[0]
        tmp_dict = json.load(open(f'{indir}/{i}'))
        tmpbase = tmp_dict['summary']['before_filtering']['total_bases']
        tmpread = tmp_dict['summary']['before_filtering']['total_reads']
        tmpgc = tmp_dict['summary']['before_filtering']['gc_content']
        tmplen = tmp_dict['summary']['before_filtering']['read1_mean_length']
        tmpq20 = tmp_dict['summary']['before_filtering']['q20_rate']
        tmpq30 = tmp_dict['summary']['before_filtering']['q30_rate']
        tmpfp = tmp_dict['filtering_result']['passed_filter_reads']
        open(f'{outdir}/{samname}_result_base.txt','w').write(f'样本名称\t总碱基数\t总序列数\tgc含量\t平均长度\tQ20\tQ30\t过滤后数量\n')
        open(f'{outdir}/{samname}_result_base.txt','a').write(f'{samname}\t{tmpbase}\t{tmpread}\t{tmpgc}\t{tmplen}\t{tmpq20}\t{tmpq30}\t{tmpfp}\n')
        open(f'{outdir}/result_base.txt','a').write(f'{samname}\t{tmpbase}\t{tmpread}\t{tmpgc}\t{tmplen}\t{tmpq20}\t{tmpq30}\t{tmpfp}\n')
        #with open(f'{wkdir}/multiqc.log','w') as mlf:
        #    subprocess.run(f'multiqc {indir}/{samname}.fastp.json -i {samname} -o {outdir} -f --quiet',shell=True,stdout=mlf,stderr=mlf)
    #with open(f'{wkdir}/multiqc.log','w') as mlf:
    #    subprocess.run(f'multiqc {indir}/*.fastp.json -o {outdir} -f --quiet',shell=True,stdout=mlf,stderr=mlf)
#------mapping-----
def summ_map(indir,outdir):
    mapping_dic ={}
    qua_list = [i for i in os.listdir(indir) if 'qua' in i]
    for i in qua_list:
        name=i.replace('_quaimap','')
        #with open(f'{wkdir}/multiqc.log','w') as mlf:
        #    subprocess.run(f'multiqc {indir}/{name}_quaimap -o {outdir} -i {name} -f --quiet',shell=True,stdout=mlf,stderr=mlf)
        mappat = re.compile('(\S+) \((\S+)\)')
        qfile = os.path.abspath(f'{indir}/{i}/genome_results.txt')
        with open(qfile) as f:
            mapping_dic[name]={}
            for line in f:
                line = line.strip()
                if 'number of reads' in  line:
                    mapping_dic[name]['clean_fq_reads_nums'] = line.split('=')[-1]
                elif 'number of mapped reads' in line:
                    mapping_dic[name]['mapping_reads_nums'] = re.match(mappat,line.split('=')[-1].strip()).group(1)
                    mapping_dic[name]['mapping_rate'] = re.match(mappat,line.split('=')[-1].strip()).group(2)
                elif 'GC percentage' in line:
                    mapping_dic[name]['GC_rate'] = line.split('=')[-1]
                elif 'mean coverageData' in line:
                    mapping_dic[name]['mean_depth'] = line.split('=')[-1]
                    
    df1 = pd.DataFrame(mapping_dic)
    if df1.shape[0] != 0:
        df2 = pd.DataFrame(df1.values.T, index=df1.columns, columns=df1.index) 
        df2['samplename'] = df2.index
        df2 = df2[['samplename','mapping_reads_nums','mapping_rate','GC_rate','mean_depth']]
        df2.rename(columns={'samplename':'样本名称','mapping_reads_nums':'比对序列数','mapping_rate':'比对率','GC_rate':'GC含量','mean_depth':'平均深度'},inplace=True)
        for Sam in df2['样本名称'].to_list():
            tmpdf2 = df2[df2['样本名称']==Sam]
            tmpdf2.to_csv(f'{outdir}/{Sam}_mapping_summary.tsv',sep='\t',index=False)
#-----tbp------
def summ_tbp(indir,outdir):
    os.chdir(indir)
    twkdir = indir.replace(' 3.Tb_profiler','')
    with open('tbc.log','w') as tbcf:
        subprocess.run(f'/home/dell/miniconda3/bin/conda run -no-catpure-output -n TB-profiler tb-profiler collate',shell=True,stdout=tbcf,stderr=tbcf)
    if os.path.isfile('tbprofiler.txt') and os.path.getsize('tbprofiler.txt') != 0:
        tbpdb = pd.read_table('tbprofiler.txt',dtype={'sample':str})
        tbpdb.fillna('_',inplace=True)
        tbpdb.to_csv('tbprofiler.txt',index=False,sep='\t')
        tbpdb = tbpdb[['sample','drtype','main_lineage']]
        tbpdb.rename(columns={'sample':'样本名称','drtype':'抗药类型','main_lineage':'主家系'},inplace=True)
        tbpdb['样本名称'] = tbpdb['样本名称'].astype('str')
        tbpdb['样本名称老'] = tbpdb['样本名称']
        for tname in tbpdb['样本名称'].tolist():
            tbpdb.loc[tbpdb['样本名称']==tname,'抗药类型'] = os.popen(f'''awk -F '\t'  '{{if($1=="{tname}")print $2}}' ../1.snp_calling/all_sample_dr.tsv''').read().strip()
        tbpdb = tbpdb[['样本名称老','样本名称','抗药类型','主家系']]
        tbpdb.to_csv('tbpro.tsv',sep='\t',index=False)
    else:
        open('tbpro.tsv','w').write(f'样本名称老\t样本名称\t抗药类型\t主家系\n')
        with open(f'{wkdir}/Samplelist.txt') as f:
            for line in f:
                Sam = line.strip()
                open(f'tbpro.tsv','a').write(f'{Sam}\t{Sam}\t-\t-\n')

#-----vcf-----
def summ_vcf(indir,outdir):
    nu = pd.DataFrame(columns=['sample_name','del','snp','complex','ins','mnp'])
    for i in os.listdir(indir):
        if os.path.isdir(f'{indir}/{i}'):
            samn = i
            filedir = os.path.abspath(f'{indir}/{i}')
            tabfile = f'{filedir}/snps.tab'
            tbf = pd.read_table(tabfile)
            tmp1 = pd.DataFrame(tbf.iloc[:,2].value_counts()).T
            tmp1['sample_name'] = samn
            nu = pd.concat([tmp1,nu],axis=0,sort=True)
            nu = nu.sort_values(by = "sample_name")
    nu.to_csv(f'{outdir}/vcf_stat.csv',sep='\t',index=False,columns=['sample_name','del','snp','complex','ins','mnp'],na_rep='0')

 
def fast5_func(indir,outdir,mode,Pre):
    if mode == 'hac':
        gconf='dna_r9.4.1_450bps_hac.cfg'
    else:
        gconf='dna_r9.4.1_450bps_fast.cfg'
    subprocess.run(f'guppy_basecaller -i {indir} -s {outdir} -c {gconf} -x auto',shell=True)
    subprocess.run(f'cat {outdir}/pass/*.fastq > {outdir}/{Pre}.fq',shell=True)


def All_lin_Spe(ind):
    os.chdir(ind)
    if os.path.isfile('Samplelist.txt'):
        with open('Samplelist.txt') as f:
            if os.path.isfile(f'3.Tb_profiler/tbprofiler.txt') and os.path.getsize(f'3.Tb_profiler/tbprofiler.txt') != 0:
                tbprodb = pd.read_table('3.Tb_profiler/tbprofiler.txt', dtype={'sample': str})
                tbprodb['sample'] = tbprodb['sample'].astype('str')
            else:
                tbprodb = 0
            #allsmdb = pd.read_table('1.snp_calling/all_sample_dr.tsv',header=None,names=['sample','Drtype'])
            #allsmdb['sample'] = allsmdb['sample'].astype('str')
            allvfdb = pd.read_table('1.snp_calling/vcf_stat.csv', dtype={'sample_name': str})
            allvfdb['sample_name'] = allvfdb['sample_name'].astype('str')
            for Sam in f:
                #try:
                Sam = Sam.strip()
                MainSpe,readsnum,Abun = os.popen(f'head -n 2 4.kraken_taxonomic/{Sam}.tb.report1.txt|tail -n1|cut -f1,3,4').read().strip().split('\t')
                if isinstance(tbprodb, pd.DataFrame):
                    if Sam in tbprodb['sample'].tolist():
                        MainLin = tbprodb[tbprodb['sample']==Sam]['main_lineage'].tolist()[0]
                        Sublin = tbprodb[tbprodb['sample']==Sam]['sub_lineage'].tolist()[0]
                    else:
                        MainLin = '-'
                        Sublin = '-'
                else:
                    MainLin = '-'
                    Sublin = '-'
                Drtype = os.popen(f'''awk -F '\t'  '{{if($1=="{Sam}")print $2}}' 1.snp_calling/all_sample_dr.tsv''').read().strip()
                DelNum = allvfdb[allvfdb['sample_name']==Sam]['del'].tolist()[0]
                SnpNum = allvfdb[allvfdb['sample_name']==Sam]['snp'].tolist()[0]
                InsNum = allvfdb[allvfdb['sample_name']==Sam]['ins'].tolist()[0]
                if os.path.isfile(f'0.QC/{Sam}_mapping_summary.tsv'):
                    maprate,covDep = os.popen(f'grep {Sam} 0.QC/{Sam}_mapping_summary.tsv|cut -f3,5').read().strip().split('\t')
                else:
                    maprate,covDep = ['-','-']
                open(f'{Sam}_sum_result.txt','w').write(f'''物种鉴定（序列数量/丰度）：\t{MainSpe}({readsnum}/{Abun})

家系（子家系）鉴定：\t{MainLin}({Sublin})
耐药类型：\t{Drtype}
比对率：\t{maprate}
平均覆盖深度：\t{covDep}
变异检测统计：\t插入{InsNum}、缺失{DelNum}、单碱基突变{SnpNum}
                ''')


def QC_sum(ind):
    os.chdir(ind)
    if os.path.isfile('Samplelist.txt'):
        with open('Samplelist.txt') as f:
            for Sam in f:
                Sam = Sam.strip()
                if os.path.isfile(f'fq_file/{Sam}_result_base.txt'):
                    fqdb = pd.read_table(f'fq_file/{Sam}_result_base.txt')
                    mapdb = pd.read_table(f'0.QC/{Sam}_mapping_summary.tsv')
                    numfqampdb = fqdb.merge(mapdb,on='样本名称')
                    numfqampdb = numfqampdb.T
                    numfqampdb.drop('GC含量',inplace=True)
                    numfqampdb['质控标准'] = ['-','-','-','0.63~0.68','-','>0.95','>0.9','-','-','>0.95','>20x']
                    numfqampdb['质控结果'] = '通过'
                    #- 
                    if numfqampdb.loc[numfqampdb.index=='gc含量',0].tolist()[0] < 0.63 or numfqampdb.loc[numfqampdb.index=='gc含量',0].tolist()[0] > 0.68:
                        numfqampdb.loc[numfqampdb.index=='gc含量','质控结果'] = '失败'
                    if numfqampdb.loc[numfqampdb.index=='Q20',0].tolist()[0] < 0.95:
                        numfqampdb.loc[numfqampdb.index=='Q20','质控结果']  = '失败'
                    if numfqampdb.loc[numfqampdb.index=='Q30',0].tolist()[0] < 0.90:
                        numfqampdb.loc[numfqampdb.index=='Q30','质控结果']  = '失败'
                    if float(numfqampdb.loc[numfqampdb.index=='比对率',0].tolist()[0].replace('%','')) < 95:
                        numfqampdb.loc[numfqampdb.index=='比对率','质控结果']  = '失败'
                    if float(numfqampdb.loc[numfqampdb.index=='平均深度',0].tolist()[0].replace('X','').replace(',','')) < 0.95:
                        numfqampdb.loc[numfqampdb.index=='平均深度','质控结果']  = '失败'
                    numfqampdb['质控项'] =numfqampdb.index
                    numfqampdb.rename(columns={0:'测量值'},inplace=True)
                    numfqampdb = numfqampdb[['质控项','质控标准','测量值','质控结果']]
                    numfqampdb.to_csv(f'{Sam}_QCout.tsv',sep='\t',index=False)





#-------主流程-------
#----质控-----
#-----整理kraken2注释文件-----
Stime = time.time()
chin={}
lth={}
info={}
faslist=['.fasta','.fna','.fa','.fsa','.mpfa','.fas']
with open(listt) as f:
    for line in f.readlines():
        line = line.strip().split('\t')
        info[line[0]] = f'{line[3]}\t{line[4]}'
        lth[line[0]] = line[-1]
        chin[line[0]] = line[2]
#------------建立所有文件夹---------------
newstart = 1
if os.path.isfile('0.QC'):
    newstart = 0
for i in ['0.QC','1.snp_calling','2.Tree','3.Tb_profiler','4.kraken_taxonomic','5.virulence_genes','6.Snpit','7.snpcluster','fq_file','ref']:
    if not os.path.exists(i):
        os.mkdir(i)

#------snp call 质控-------
Pro_fun('ref_process(ref,gff,f\'{wkdir}/ref\')','1.准备参考基因组')
Allnum = os.popen(f'cat {Listn}|wc -l').read().strip()
Samnum = 1

stime = time.time()
s1time = 0
if not os.path.isfile('vcflist.txt'):
    open('vcflist.txt','w').write(f'样本\t数据路径')
if not os.path.isfile('bamlist.txt'):
    open('bamlist.txt','w').write(f'样本\t数据路径')

if not os.path.isfile(f'{wkdir}/Samplelist.txt'):
    open(f'{wkdir}/Samplelist.txt','w').write('')
if not os.path.isfile(f'{wkdir}/trim_fqlist'):
    open(f'{wkdir}/trim_fqlist','w').write(f'样本\t数据路径左\t数据路径右\t数据类型\n')

with open(Listn) as listf:
    for line in listf:
        line = line.strip().split('\t')
        if len(line) == 3:
            gs=2
            Pre,fq1,fq2 = line
            print(f'task_step：0/5\t样本进度：{Samnum}/{Allnum}\t样本：{Pre}\t开始分析')
            sys.stdout.flush()
            open(f'{wkdir}/{Pre}_gs_info.txt','w').write('illumina二代测序平台，依据WHO2023年发布的结核耐药点位标准文件')
            print(f'{Pre}为双端测序结果')
            #qc_funt(f'{wkdir}/0.QC',Pre,nT,fq1,fq2)
            if not os.path.isfile(f'{wkdir}/fq_file/{Pre}_OK'):
                print(f'{Pre} not ok')
                try:
                    Pro_fun('trim_fun(f\'{wkdir}/fq_file\',nT,Pre,2,fq1,fq2)',f'{Pre}_数据质控')
                except:
                    print(f'{Pre} 样本trim失败，请检测原始文件输入是否正确')
                    sys.stdout.flush()
            #qc_funt(f'{wkdir}/fq_file',f'{Pre}_clean',nT,f'{wkdir}/fq_file/{Pre}_clean_1.gz',f'{wkdir}/fq_file/{Pre}_clean_2.gz')
        else:
            Pre,fq1 = line
            print(f'task_step：0/5\t样本进度：{Samnum}/{Allnum}\t样本：{Pre}\t开始分析')
            sys.stdout.flush()
            open(f'{wkdir}/{Pre}_gs_info.txt','w').write('三代测序平台，依据WHO2021年发布的结核耐药点位标准文件')
            if not any([i for i in faslist if fq1.endswith(i)]):
                if os.path.isfile(fq1):
                    if fq1.endswith('fq') or fq1.endswith('fq.gz') or fq1.endswith('fastq') or fq1.endswith('fastq.gz'):
                        gs=3
                        print(f'{Pre}为单端测序结果')
                        #qc_funt(f'{wkdir}/0.QC',Pre,nT,fq1)
                        #ofn = f'{wkdir}/ont_fq'
                        #if not os.path.isdir(ofn):
                        #    os.makedirs(ofn)
                        #subprocess.run(f'cat {fq1} |seqkit sliding -s 1000 -W 1000 -g > {ofn}/{Pre}.fq',shell=True)
                        #if not os.path.isfile(f'{wkdir}/0.QC/{Pre}_OK'):
                        #    try:
                        trim_fun(f'{wkdir}/fq_file',nT,Pre,3,fq1)
                        #    except:
                        #        print(f'{Pre} 样本trim失败，请检测原始文件输入是否正确')
                        #        sys.stdout.flush()
                    else:
                        if fq1.endswith('vcf') or fq1.endswith('vcf.gz'):
                            print(f'{Pre} 为vcf格式文件')
                            sys.stdout.flush()
                            vcflistdb = pd.read_table(f'{wkdir}/trim_fqlist')
                            tmpvcfdict =  {'样本':f'{Pre}_2','数据路径左':fq1,'数据路径右':'-','数据类型':'vcf'}
                            if not f'{Pre}_2'  in vcflistdb['样本'].tolist():
                                vcflistdb.loc[len(vcflistdb)] = tmpvcfdict
                                vcflistdb.to_csv(f'{wkdir}/trim_fqlist',sep='\t',index=False)
                        elif fq1.endswith('bam'):
                            print(f'{Pre} 为bam格式文件')
                            sys.stdout.flush()
                            bamlistdb = pd.read_table(f'{wkdir}/trim_fqlist')
                            tmpbamdict = {'样本':f'{Pre}_2','数据路径左':fq1,'数据路径右':'-','数据类型':'bam'}
                            if not f'{Pre}_2'  in bamlistdb['样本'].tolist():
                                bamlistdb.loc[len(bamlistdb)] = tmpbamdict
                                bamlistdb.to_csv(f'{wkdir}/trim_fqlist',sep='\t',index=False)






                    #qc_funt(f'{wkdir}/fq_file',f'{Pre}_clean',nT,f'{wkdir}/fq_file/{Pre}_clean.gz')
                elif os.path.isdir(fq1):
                    gs=3
                    print(f'{Pre}为三代测序结果')
                    Ont_dir=os.path.abspath(fq1)
                    filetmp = os.listdir(Ont_dir)[0]
                    ofn=f'{wkdir}/fast5_output/'
                    fq_new=f'{ofn}/{Pre}.fq'
                    if not os.path.isdir(ofn):
                        os.mkdir(ofn)
                    if 'fast5' in filetmp:
                        fast5_func(Ont_dir,ofn,Mode,Pre)
                        if not os.path.isfile(f'{wkdir}/0.QC/{Pre}_OK'):
                            try:
                                trim_fun(f'{wkdir}/fq_file',nT,Pre,3,fq_new)
                            except:
                                print(f'{Pre} 样本trim失败，请检测原始文件输入是否正确')
                                sys.stdout.flush()
                    else:
                        if not os.path.isfile(f'{wkdir}/0.QC/{Pre}_OK'):
                            try:
                                subprocess.run(f'cat {Ont_dir}/* |seqkit sliding -s 1000 -W 1000 -g > {ofn}/{Pre}.fq',shell=True)
                                trim_fun(f'{wkdir}/fq_file',nT,Pre,3,fq_new)
                            except:
                                print(f'{Pre} 样本trim失败，请检测原始文件输入是否正确')
                                sys.stdout.flush()
            else:
                gs=2
                print(f'{Pre}是fasta文件')
                ofn=f'{wkdir}/fake_fastq_output/'
                if not os.path.isdir(ofn):
                    os.mkdir(ofn)
                if not os.path.isfile(f'{wkdir}/0.QC/{Pre}_OK'):
                    try:
                        fa_process(fq1,ofn,Pre)
                        Pro_fun('trim_fun(f\'{wkdir}/fq_file\',nT,Pre,3,f\'{ofn}/{Pre}.fq\')',f'{Pre}_数据质控')
                    except:
                        print(f'{Pre} 样本trim失败，请检测原始文件输入是否正确')
                        sys.stdout.flush()
        if s1time == 0:
            s1run = time.time() - stime
            s1time = time.time()
        else:
            s1run = time.time() - s1time
            s1time = time.time()
        s1runtime = format_seconds(s1run)
        print(f'task_step：1/5\t样本进度：{Samnum}/{Allnum}\t样本：{Pre}\t质控分析已结束\t运行时间: {s1runtime}')
        sys.stdout.flush()
        Samnum += 1 
#-----mapping callsnp quailmap----------
all_vcf = {}
all_vcf['pass'] = []
all_vcf['nopass'] = []
Samnum = 1
S2time = 0
S3time = 0 
with open(f'{wkdir}/trim_fqlist') as trimf:
    for line in trimf:
        if not line.startswith('样本'):
            line = line.strip().split('\t')
            if line[2] != '-':
                # 双端数据
                Pre_t,fq1,fq2,inputype = line
                Prelist = Pre_t.split('_')
                Prelist.pop()
                Pre = '_'.join(Prelist)
                gs = int(Pre_t.split('_')[-1])
                if not os.path.isfile(f'{wkdir}/1.snp_calling/{Pre}/map_OK'):
                    print('not OK？')
                    try:
                        fq_process(gs,Pre,inputype,f'{wkdir}/ref/ref.fa',nT,vcfdep,f'{wkdir}/1.snp_calling',fq1,fq2)
                    except:
                        print(f'{Pre} 样本生成VCF失败，请核实原因')
                        sys.stdout.flush()
                if S2time == 0:
                    s2run = time.time()-s1time
                    S2time = time.time()
                else:
                    s2run = time.time()-S2time
                    S2time = time.time()
                s2runtime = format_seconds(s2run)
                print(f'task_step：2/5\t样本进度：{Samnum}/{Allnum}\t样本：{Pre}\t突变分析已结束\t运行时间: {s2runtime}')
                sys.stdout.flush()
                with open('tb.log','w') as f1:
                    try:
                        if not os.path.isfile(f'{wkdir}/3.Tb_profiler/results/{Pre}.results.json'):
                            subprocess.run(f'/home/dell/miniconda3/bin/conda run -n TB-profiler tb-profiler profile --vcf {wkdir}/1.snp_calling/{Pre}/snps.vcf.gz -d {wkdir}/3.Tb_profiler --prefix {Pre}',shell=True,stdout=f1,stderr=f1)
                    except:
                        print(f'{Pre} 家系鉴定失败')
                        sys.stdout.flush()
                if not os.path.isfile(f'{wkdir}/4.kraken_taxonomic/{Pre}_ok'):
                    kk2(nT,Pre,fq1,fq2)
                S3time = time.time()
                if S3time == 0:
                    s3runtime = format_seconds(S3time-S2time)
                    S3time = time.time()
                else:
                    s3runtime = format_seconds(time.time()-S3time)
                    S3time = time.time()
                print(f'task_step：3/5\t样本进度：{Samnum}/{Allnum}\t样本：{Pre}\t物种鉴定已结束\t运行时间: {s3runtime}')
                sys.stdout.flush()
            else:
                # 单端或者其他类型数据
                Pre_t,fq1,fq2,inputype = line
                Prelist = Pre_t.split('_')
                Prelist.pop()
                Pre = '_'.join(Prelist)
                gs = int(Pre_t.split('_')[-1])
                try:
                    fq_process(gs,Pre,inputype,f'{wkdir}/ref/ref.fa',nT,vcfdep,f'{wkdir}/1.snp_calling',fq1)
                except:
                    print(f'{Pre} 样本生成VCF失败，请核实原因')
                    sys.stdout.flush()
                if inputype != 'fastq':
                    if int(os.popen(f'cat {wkdir}/Samplelist.txt|wc -l').read().strip()) != 0:
                        Slistdb = pd.read_table(f'{wkdir}/Samplelist.txt',header=None)
                        if Pre not in Slistdb[0].tolist():
                            stmpdict = {0:Pre}
                            Slistdb.loc[len(Slistdb)] = stmpdict
                            Slistdb.to_csv(f'{wkdir}/Samplelist.txt',sep='\t',index=False,header=None)
                    else:
                        open(f'{wkdir}/Samplelist.txt','a').write(f'{Pre}\n')
                if S2time == 0:
                    s2run = time.time()-s1time
                    S2time = time.time()
                else:
                    s2run = time.time()-S2time
                    S2time = time.time()
                s2runtime = format_seconds(s2run)
                print(f'task_step：2/5\t样本进度：{Samnum}/{Allnum}\t样本：{Pre}\t突变分析已结束\t运行时间: {s2runtime}')
                sys.stdout.flush()
                #subprocess.run(f'tb-profiler vcf_profile {wkdir}/1.snp_calling/{Pre}/snps.vcf.gz -d {wkdir}/3.Tb_profiler --reporting_af 0',shell=True)
                with open('tb.log','w') as f1:
                    try:
                        subprocess.run(f'/home/dell/miniconda3/bin/conda run -n TB-profiler tb-profiler profile --vcf {wkdir}/1.snp_calling/{Pre}/snps.vcf.gz -d {wkdir}/3.Tb_profiler --prefix {Pre}',shell=True,stdout=f1,stderr=f1)
                    except:
                        print(f'{Pre} 家系鉴定失败')
                if not os.path.isfile(f'{wkdir}/4.kraken_taxonomic/{Pre}_ok'):
                    if inputype == 'fastq':
                        kk2(nT,Pre,fq1)
                    else:
                        open(f'{wkdir}/4.kraken_taxonomic/{Pre}.tb.report1.txt','w').write(f'物种\t序列数量\t相对丰度\ttaxid\n结核分枝杆菌\t-\t-\tMycobacterium tuberculosis')
                        open(f'{wkdir}/4.kraken_taxonomic/{Pre}.list2.txt','w').write(f'物种\t序列数量\t比例\tNCBI物种号')
                        open(f'{wkdir}/4.kraken_taxonomic/{Pre}.tb.report2.txt','w').write(f'物种\t序列数量\t比例\tNCBI物种号')
                S3time = time.time()
                s3runtime = format_seconds(S3time-S2time)
                print(f'task_step：3/5\t样本进度：{Samnum}/{Allnum}\t样本：{Pre}\t物种鉴定已结束\t运行时间: {s3runtime}')
                sys.stdout.flush()
            if int(os.popen(f'wc -l {wkdir}/1.snp_calling/{Pre}/snps.drsubs.vcf').read().strip().split(' ')[0]) >= 0:
                all_vcf['pass'].append(Pre)
                subprocess.run(f'cp {wkdir}/1.snp_calling/{Pre}/snps.subs.drsubs.consensus.fa {wkdir}/1.snp_calling/{Pre}/snps.aligned.fa',shell=True)
            else:
                all_vcf['nopass'].append(Pre)
            with open(f'{wkdir}/0.QC/bamqc.log','a') as f:
                if not os.path.isfile(f'{wkdir}/0.QC/{Pre}_ok'):
                    if inputype != 'vcf':
                        try:
                            subprocess.run(f'qualimap bamqc -bam  {wkdir}/1.snp_calling/{Pre}/snps.bam --java-mem-size=8G --outdir {wkdir}/0.QC/{Pre}_quaimap',shell=True,stderr=f,stdout=f)
                            open(f'{wkdir}/0.QC/{Pre}_ok','w').write('')
                        except:
                            print(f'{Pre} 比对质控失败')
                            sys.stdout.flush()
            S4time = time.time()
            s4runtime = format_seconds(S4time-S3time)
            print(f'task_step：4/5\t样本进度：{Samnum}/{Allnum}\t样本：{Pre}\t比对质控已结束\t运行时间: {s3runtime}')
            sys.stdout.flush()
            Samnum += 1 
os.chdir(f'{wkdir}/1.snp_calling')
pass_vcflist = " ".join(all_vcf['pass'])
nopass_vcflist = " ".join(all_vcf['nopass'])
if len(all_vcf['nopass'])>0:
    print(f'{nopass_vcflist} 这些样本没有检测到变异位点,不做后续分析')
with open(f'{wkdir}/snippy.log','w') as snpf:
    subprocess.run(f'{snippydir}/bin/snippy-core --mask {black_list} --prefix all_sample --ref {wkdir}/ref/ref.fa {pass_vcflist}',shell=True,stdout=snpf,stderr=snpf)
    subprocess.run(f'{snippydir}/bin/snippy-core --prefix all_sample_raw --ref {wkdir}/ref/ref.fa {pass_vcflist}',shell=True,stdout=snpf,stderr=snpf)
subprocess.run(f'''cat {wkdir}/1.snp_calling/*/sam_dr.tsv > {wkdir}/1.snp_calling/all_sample_dr.tsv ''',shell=True)
#------毒力分析--------
#------abricate1.0.1 database.updata.2021-Mar-27------
#subprocess.run(f'ls  {wkdir}/1.snp_calling/*/snps.consensus.fa > {wkdir}/5.virulence_genes/samplelist',shell=True)
#open(f'{wkdir}/5.virulence_genes/ere_virulence_genes.tsv','w').write(f'文件名称\t染色体\t开始位置\t终止位置\t正负链\t基因名称\t覆盖度\t覆盖度展示\t空缺\t覆盖百分比\t%IDENTITY\t数据库\tACCESSION\t产物\n')
#with open(f'abricate.log','w') as abf:
#    subprocess.run(f'abricate --db vfdb --fofn {wkdir}/5.virulence_genes/samplelist --threads {nT} |sed \'1d\' >> {wkdir}/5.virulence_genes/tmp_Pre_virulence_genes.tsv',shell=True,stdout=abf,stderr=abf)
#with open(f'{wkdir}/5.virulence_genes/tmp_Pre_virulence_genes.tsv') as f:
#    for line in f:
#        line = line.strip().split('\t')
#        if line[0]!="文件名称":
#            sam = line[0].split('/')[-2]
#            line[0] = sam
#            line = '\t'.join(line)
#            line = f'{line}\n'
#            open(f'{wkdir}/5.virulence_genes/Pre_virulence_genes.tsv','a').write(line)
#subprocess.run(f'rm {wkdir}/5.virulence_genes/tmp_Pre_virulence_genes.tsv',shell=True)
#-----snpit-------
#subprocess.run(f'cp {wkdir}/1.snp_calling/all_sample.full.aln {wkdir}/6.Snpit/all_sample.full.fa',shell=True)
#subprocess.run(f'/data/wusihao/soft/snpit/.venv/bin/snpit -i {wkdir}/1.snp_calling/all_sample.vcf > {wkdir}/6.Snpit/Snpit_vcf.result',shell=True)
#subprocess.run(f'/data/wusihao/soft/snpit/.venv/bin/snpit -i {wkdir}/6.Snpit/all_sample.full.fa > {wkdir}/6.Snpit/Snpit_fa.result',shell=True)
#subprocess.run(f'rm {wkdir}/6.Snpit/all_sample.full.fa',shell=True)
#----def All
def SumAll(ind):
    sumalldict = {}
    os.chdir(ind)
    drug_dict = {'阿米卡星:amikacin':'amikacin(AMI)','贝达喹啉:bedaquiline':'bedaquiline(BDQ)','卷曲霉素:capreomycin':'Capreomycin(CAP)','德拉马尼:delamanid':'delamanid(DLM)','乙胺丁醇:ethambutol':'ethambutol(EMB)','乙硫异烟胺:ethionamide':'ethionamide(ETO)','异烟肼:isoniazid':'isoniazid(INH)','卡那霉素:kanamycin':'kanamycin(KAN)','利奈唑胺:linezolid':'inezolid(LSD)','莫西沙星:moxifloxacin':'moxifloxacin(MFX)','吡嗪酰胺:pyrazinamide':'pyrazinamide(PZA)','利福平:rifampicin':'rifampicin(RIF)','链霉素:streptomycin':'streptomycin(STM)','左氧氟沙星:levofloxacin':'levofloxacin(LEV)','环丝氨酸:cycloserine':'cycloserine(CS)','氯法齐明:Clofazimine':'clofazimine(CFZ)'}
    trimdb = pd.read_table('trim_fqlist')
    for Pre in trimdb['样本'].tolist():
        Pre = '_'.join(Pre.split('_')[:-1])
        print(Pre)
        sumalldict[Pre] = {'测序平台':'-','测序方案':'-','测序数据量(G)':'-','平均深度':'-','质控结果':'-','主家系':'-','子家系':'-','耐药类型':'-','rifampicin(RIF)':'','isoniazid(INH)':'','moxifloxacin(MFX)':'','ethambutol(EMB)':'','linezolid(LZD)':'','clofazimine(CFZ)':'','bedaquiline(BDQ)':'','cycloserine(CS)':'','kanamycin(KAN)':'','capreomycin(CAP)':'','para-aminosalicylic_acid':'','streptomycin(STM)':'','ethionamide(ETO)':'','pyrazinamide(PZA)':'','delamanid(DLM)':'','备注':''}
        #---1
        if os.path.isfile(f'{Pre}_QCout.tsv'):
            QCdb = pd.read_table(f'{Pre}_QCout.tsv')
            sumalldict[Pre]['测序数据量(G)']= round(int(QCdb.loc[QCdb['质控项']=='总碱基数','测量值'].tolist()[0])/1000000000,2)
            sumalldict[Pre]['平均深度']= QCdb.loc[QCdb['质控项']=='平均深度','测量值'].tolist()[0]
            failnum = QCdb[QCdb['质控结果']=='失败'].shape[0]
            if failnum >2:
                sumalldict[Pre]['质控结果']='fail'
            elif failnum >= 1 and failnum<=2:
                sumalldict[Pre]['质控结果']='warning'
            else:
                sumalldict[Pre]['质控结果']='pass'
            sumalldict[Pre]['GC质控']= QCdb.loc[QCdb['质控项']=='gc含量','测量值'].tolist()[0]
            sumalldict[Pre]['Q20质控']= QCdb.loc[QCdb['质控项']=='Q20','测量值'].tolist()[0]
            sumalldict[Pre]['Q30质控']= QCdb.loc[QCdb['质控项']=='Q30','测量值'].tolist()[0]

        #---2
        if os.path.isfile(f'{Pre}_sum_result.txt'):
            sumrdb =  pd.read_table(f'{Pre}_sum_result.txt',header=None)
            sumalldict[Pre]['主家系'] = sumrdb.iloc[1,1].split('(')[0]
            sumalldict[Pre]['子家系'] = sumrdb.iloc[1,1].split('(')[1].replace(')','')
            sumalldict[Pre]['耐药类型'] = sumrdb.iloc[2,1]
            sumalldict[Pre]['物种'] = sumrdb.iloc[0,1]
        #---3
        if os.path.isfile(f'1.snp_calling/{Pre}/Chin_snpdr.tsv'):
            samdrdb = pd.read_table(f'1.snp_calling/{Pre}/Chin_snpdr.tsv')
            samdrdb = samdrdb[samdrdb['相关性'].isin(['1) Assoc w R','2) Assoc w R - Interim'])]
            if samdrdb.shape[0]>0:
                for Drug in samdrdb['药物名'].unique().tolist():
                    newDrug = drug_dict.get(Drug,'-')
                    if newDrug!='-':
                        sumalldict[Pre][newDrug]=';'.join(samdrdb.loc[samdrdb['药物名']==Drug,'突变结果'].unique().tolist())
    sumalldb = pd.DataFrame(sumalldict).T.reset_index()
    sumalldb = sumalldb[['index','测序平台','测序方案','测序数据量(G)','平均深度','质控结果','GC质控','Q20质控','Q30质控','物种','主家系','子家系','耐药类型','rifampicin(RIF)','isoniazid(INH)','moxifloxacin(MFX)','ethambutol(EMB)','linezolid(LZD)','clofazimine(CFZ)','bedaquiline(BDQ)','cycloserine(CS)','kanamycin(KAN)','capreomycin(CAP)','para-aminosalicylic_acid','streptomycin(STM)','ethionamide(ETO)','pyrazinamide(PZA)','delamanid(DLM)','备注']].fillna('-')
    sumalldb.rename(columns={'index':'样本名称'},inplace=True)
    sumalldb.to_csv(f'CombineAll.tsv',sep='\t',index=False)

#-----------------总结结果信息--------------
summ_map(f'{wkdir}/0.QC',f'{wkdir}/0.QC')
summ_fq(f'{wkdir}/fq_file',f'{wkdir}/fq_file')
summ_tbp(f'{wkdir}/3.Tb_profiler',f'{wkdir}/3.Tb_profiler')
summ_vcf(f'{wkdir}/1.snp_calling',f'{wkdir}/1.snp_calling')
etime = time.time()
numSam = os.popen(f'cat {wkdir}/trim_fqlist|wc -l').read().strip()
numdb = pd.read_table(f'{wkdir}/fq_file/result_base.txt')
open(f'{wkdir}/sample_results.txt','w').write(f'''运行时间\t样本数量\t总数据量\t总序列数\tQ20数据占比
{etime-stime}\t{numSam}\t{numdb['总碱基数'].sum()}\t{numdb['总序列数'].sum()}\t{numdb['Q20'].mean()}\t{numdb['Q30'].mean()}
''')
subprocess.run(f'/home/dell/miniconda3/bin/conda run -n report_env Rscript /data/deploy/TB_soft/report.R {wkdir}',shell=True)
os.chdir(wkdir)
All_lin_Spe(wkdir)
QC_sum(wkdir)
SumAll(wkdir)
if not os.path.isdir('main_results'):
    os.makedirs('main_results')
subprocess.run(f'cp *.html main_results main_results',shell=True)
subprocess.run(f'cp 0.QC/*mapping_summary.tsv main_results',shell=True)
subprocess.run(f'cp fq_file/*result_base.txt main_results',shell=True)
if os.path.isfile(f'3.Tb_profiler/tbprofiler.txt') and os.path.getsize('3.Tb_profiler/tbprofiler.txt') !=0:
    tbprodb = pd.read_table('3.Tb_profiler/tbprofiler.txt',dtype={'sample': str})
    tbprodb = tbprodb[['sample','main_lineage','sub_lineage']]
    tbprodb.rename(columns={'sample':'样本名称','main_lineage':'主家系','sub_lineage':'子家系'},inplace=True)
    tbprodb.to_csv('main_results/lineage.tsv',sep='\t',index=False)
else:
    open(f'main_results/lineage.tsv','w').write(f'样本名称\t主家系\t子家系\n')
    with open('Samplelist.txt') as f:
        for line in f:
            Sam = line.strip()
            open(f'main_results/lineage.tsv','a').write(f'{Sam}\t-\t-\n')

for i in os.listdir('1.snp_calling'):
    if os.path.isdir(i):
        tsam = i.strip()
        subprocess.run(f'cp  1.snp_calling/{tsam}/Chin_snpdr.tsv main_results/{tsam}_Chin_snpdr.tsv',shell=True)
        subprocess.run(f'cp  1.snp_calling/{tsam}/snps.raw.vcf main_results/{tsam}_snp.vcf',shell=True)
        subprocess.run(f'cp 4.kraken_taxonomic/{tsam}.list2.txt  main_results/{tsam}_Spe.tsv',shell=True)
        subprocess.run(f'cp 1.snp_calling/{tsam}/snps.consensus.fa main_results/{tsam}.consensus.fasta',shell=True)

with open('combine.log','w') as cbl:
    subprocess.run(f'zip -r main_results.zip main_results',shell=True,stdout=cbl,stderr=cbl)
print(f'总运行时间:{etime-stime}')
sys.stdout.flush()
print(f'task_step：5/5\t样本进度：{Allnum}/{Allnum}\t样本：{Pre}\t分析已结束')
sys.stdout.flush()
