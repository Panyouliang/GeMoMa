import sys,glob
import os,argparse
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

GEMOMA='~/annotation/GeMoMa_pipeline/GeMoMa-1.9.jar'
java='~/miniconda3/bin/java'
mmseqs='~/miniconda3/bin/mmseqs'


currentdir = os.getcwd()
def run_shell(command):
    result = os.system(command)
    return f"{command} done, quit:{result}"


def Extractor(homolis):
    if os.path.exists('Extractor'):
        os.system('rm -rf Extractor')
    else:
        pass

    commands,works = [],[]

    with open(homolis,'r') as f:
        for line in f:
            Extract = [java, '-Xms20G', '-Xmx400G', '-jar', GEMOMA, 'CLI', 'Extractor', 'p=true', 'c=true', 'r=true', 'd=false', 'f=false']
            line = line.strip().split()
            name, gff, fna = line[0], line[1], line[2]
            works.append(name)
            Extract.extend(['a='+gff, 'g='+fna, 'outdir=Extractor/'+name])
            extrt = open(name+'.extract.sh', 'w')
            extrt.write(' '.join(x for x in Extract)+'\n')
            extrt.write('echo "$?: Extractor has done!"\n')
            extrt.close()
            os.chmod(name+'.extract.sh', 0o755)
            commands.append(name+'.extract.sh')

    with ProcessPoolExecutor(max_workers=len(commands)) as executor:
        futures = [executor.submit(run_shell, 'sh '+ command) for command in commands]


    for work in works:
        while True:
            if os.path.exists(currentdir+'/Extractor/'+work+'/assignment.tabular'):
                print(work+' extractor run done!')
                break
            else:
                time.sleep(30)


def searchhit(ref):

    if os.path.exists('Searchhit'):
        pass
    else:
        os.system('mkdir Searchhit')
    os.chdir('Searchhit')

    commands,works = [],[]
    for name in os.listdir(currentdir+'/Extractor/'):
        os.system('mkdir '+name)
        os.chdir(name)
        works.append(name)
        query = glob.glob(currentdir+'/Extractor/'+name+'/cds-parts.fasta')[0]

        refdb = [mmseqs, 'createdb', currentdir+'/'+ref, 'ref.DB']
        qrydb = [mmseqs, 'createdb', query,'query.DB']
        createindex = [mmseqs, 'convertalis', 'ref.DB', 'tmp']
        search = [mmseqs, 'search', 'query.DB', 'ref.DB', 'resultDB', 'align', '-a', '--threads', '8']
        convertalis = [mmseqs, 'convertalis', 'query.DB', 'ref.DB', 'resultDB', 'search.tabular', '--format-output', '"query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,empty,raw,nident,empty,empty,empty,qframe,tframe,qaln,taln,qlen,tlen"']

        mapped = open('mmseqs.sh','w')
        mapped.write('cd '+currentdir+'/Searchhit/'+name+'\n')
        mapped.write(' '.join(x for x in refdb)+'\n')
        mapped.write('echo "$?: make refDB"\n')
        mapped.write(' '.join(x for x in qrydb)+'\n')
        mapped.write('echo "$?: make queryDB"\n')
        mapped.write(' '.join(x for x in search)+'\n')
        mapped.write('echo "$?: search"\n')
        mapped.write(' '.join(x for x in convertalis)+'\n')
        mapped.write('echo "$?: convertalis"\n')
        mapped.close()
        os.chmod('mmseqs.sh',0o755)
        commands.append(currentdir+'/Searchhit/'+name+'/mmseqs.sh')
        os.chdir(currentdir+'/Searchhit')
    os.chdir(currentdir)

    with ProcessPoolExecutor(max_workers=3) as executor:
        futures = {executor.submit(run_shell, 'sh '+ command): command for command in commands}
        for future in as_completed(futures):
            command = futures[future]
            try:
                future.result()
                print(f"Completed: {command}")
            except Exception as exc:
                print(f"Error with {commands}: {exc}")

    for work in works:
        while True:
            if os.path.exists(currentdir+'/Searchhit/'+work+'/search.tabular'):
                print(work+' of mapping was ran done!')
                break
            else:
                time.sleep(30)


def GeneModelMapper(ref,maxintrons):

    mapper = [java, '-Xms20G', '-Xmx400G', '-jar', GEMOMA, 'CLI', 'GeMoMa', 'm='+maxintrons, 'sort=TRUE', 'outdir=output','s=search.tabular', 't=genome.fa', 'c=cds-parts.fasta','a=assignment.tabular']
    os.chdir(currentdir)
    if os.path.exists('GeneModelMapper'):
        os.system('rm -rf GeneModelMapper')

    os.system('mkdir GeneModelMapper')
    os.chdir('GeneModelMapper')

    commands,works = [],[]

    for name in os.listdir(currentdir+'/Extractor/'):
        works.append(name)
        os.system('mkdir '+name)
        os.chdir(name)
        os.symlink(currentdir+'/'+ref, currentdir+'/GeneModelMapper/'+name+'/genome.fa')
        os.symlink(currentdir+'/Searchhit/'+name+'/search.tabular', currentdir+'/GeneModelMapper/'+name+'/search.tabular')
        os.symlink(currentdir+'/Extractor/'+name+'/cds-parts.fasta', currentdir+'/GeneModelMapper/'+name+'/cds-parts.fasta')
        os.symlink(currentdir+'/Extractor/'+name+'/assignment.tabular', currentdir+'/GeneModelMapper/'+name+'/assignment.tabular')

        Gmapper = open('GeneModelMapper.sh', 'w')
        Gmapper.write('cd '+currentdir+'/GeneModelMapper/'+name+'\n')
        Gmapper.write(' '.join(x for x in mapper)+'\n')
        Gmapper.write('echo "$?: GeMoMe pipeline"\n')
        Gmapper.close()
        os.chmod('GeneModelMapper.sh', 0o755)
        commands.append(currentdir+'/GeneModelMapper/'+name+'/GeneModelMapper.sh')
        os.chdir(currentdir+'/GeneModelMapper')

    os.chdir(currentdir)

    with ProcessPoolExecutor(max_workers=3) as executor:
        futures = {executor.submit(run_shell, 'sh '+ command): command for command in commands}
        for future in as_completed(futures):
            command = futures[future]
            try:
                future.result()
                print(f"Completed: {command}")
            except Exception as exc:
                print(f"Error with {commands}: {exc}")


    for work in works:
        while True:
            if os.path.exists(currentdir+'/GeneModelMapper/'+work+'/output/predicted_annotation.gff'):
                print(work+' was ran done!')
                break
            else:
                time.sleep(30)


def GMMfilter():
    Filter = [java, '-Xms20G', '-Xmx400G', '-jar', GEMOMA, 'CLI', 'GAF', 'outdir=output', 'tf=true']

    os.chdir(currentdir)
    if os.path.exists('AnnotationFilter'):
        os.system('rm -rf AnnotationFilter')

    os.system('mkdir AnnotationFilter')

    os.chdir('AnnotationFilter')

    for name, gff in zip(os.listdir(currentdir+'/GeneModelMapper/'), glob.glob(currentdir+'/GeneModelMapper/*/output/predicted_annotation.gff')):
        GAF = ['p='+name, 'g='+gff, 'w=1']
        Filter.extend(GAF)



    clean = open('Combine_Filter.sh','w')
    clean.write(' '.join(x for x in Filter)+' f="score/aa>=1.0"'+'\n')
    clean.write(' '.join(x for x in Filter)+' f="score/aa>=2.0"'+'\n')
    clean.write(' '.join(x for x in Filter)+' f="score/aa>=3.0"'+'\n')
    clean.write(' '.join(x for x in Filter)+' f="score/aa>=4.0"'+'\n')
    clean.write(' '.join(x for x in Filter)+' f="score/aa>=5.0"'+'\n')
    clean.write('echo "$?: done!"\n')
    clean.close()
    os.chmod('Combine_Filter.sh',0o755)

    command = "sh Combine_Filter.sh >Combine_Filter.log"
    run_shell(command)
    os.chdir(currentdir)



def AnnotationFinalizer(ref,name):

    os.chdir(currentdir)
    if os.path.exists('AnnotationFinalizer'):
        os.system('rm -rf AnnotationFinalizer')

    os.system('mkdir AnnotationFinalizer')
    os.chdir('AnnotationFinalizer')


    os.symlink(currentdir+'/'+ref, currentdir+'/AnnotationFinalizer/genome.fa') 
    os.symlink(currentdir+'/AnnotationFilter/output/filtered_predictions.gff', currentdir+'/AnnotationFinalizer/filtered_predictions.gff')
    os.symlink(currentdir+'/AnnotationFilter/output/filtered_predictions_1.gff', currentdir+'/AnnotationFinalizer/filtered_predictions_1.gff')
    os.symlink(currentdir+'/AnnotationFilter/output/filtered_predictions_2.gff', currentdir+'/AnnotationFinalizer/filtered_predictions_2.gff')
    os.symlink(currentdir+'/AnnotationFilter/output/filtered_predictions_3.gff', currentdir+'/AnnotationFinalizer/filtered_predictions_3.gff')
    os.symlink(currentdir+'/AnnotationFilter/output/filtered_predictions_4.gff', currentdir+'/AnnotationFinalizer/filtered_predictions_4.gff')

    finalizes = [java, '-Xms20G', '-Xmx400G', '-jar', GEMOMA, 'CLI', 'AnnotationFinalizer', 'g=genome.fa', 'n=false', 'rename=SIMPLE', 'p='+name+'G']


    Final = open('Finalizer.sh','w')
    Final.write(' '.join(x for x in finalizes)+' a=filtered_predictions.gff outdir=Result_score2aa_1.0/'+'\n')
    Final.write(' '.join(x for x in finalizes)+' a=filtered_predictions_1.gff outdir=Result_score2aa_2.0/'+'\n')
    Final.write(' '.join(x for x in finalizes)+' a=filtered_predictions_2.gff outdir=Result_score2aa_3.0/'+'\n')
    Final.write(' '.join(x for x in finalizes)+' a=filtered_predictions_3.gff outdir=Result_score2aa_4.0/'+'\n')
    Final.write(' '.join(x for x in finalizes)+' a=filtered_predictions_4.gff outdir=Result_score2aa_5.0/'+'\n')
    Final.write('echo "$?: Finalizer done!"\n')
    Final.close()
    os.chmod('Finalizer.sh',0o755)
    command = "sh Finalizer.sh > Finalizer.log"
    run_shell(command)

    os.chdir(currentdir)


if __name__ == '__main__':

    version = "v1.0"
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''
======================================================================
This pipeline is used for call GeMoMa annotating pipeline.

Version: v1.0
Author: Panyouliang, panyouliang@genomics.cn
Date: 2024-01-06, yyyy-mm-dd
======================================================================''')
    parser.add_argument('-v', '--version', action='version', version=version)
    parser.add_argument('-ref', metavar='fasta', type=str, required=True, help='Please input the reference genome [.fa, .fasta].')
    parser.add_argument('-homolog_list', metavar='file', type=str, required=True, help='Please input the homolog_list.txt.')
    parser.add_argument('-threads', metavar='int', type=str, default='4', required=False, help='Please input the threads number, default=4, [OPTIONAL].')
    parser.add_argument('-maxintrons', metavar='int', type=str, default='100000',required=False, help='Please input the max intron [int], default=100000, [OPTIONAL].')
    parser.add_argument('-name', metavar='str', type=str,default='GeMoMa',required=False, help='Please input the specie name (OPTIONAL).')
    parser.add_argument('-MEM', metavar='int', type=str,default='20',required=False, help='Please assign the Memory size, default=20G, (OPTIONAL).')
    args = parser.parse_args()

    Extractor(args.homolog_list)
    searchhit(args.ref)
    GeneModelMapper(args.ref, args.maxintrons)
    GMMfilter()
    AnnotationFinalizer(args.ref,args.name)
