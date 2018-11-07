"""Read Alignment"""
import os
import sys
import logging
import datetime

import lib.star

import util.ddictfunc

SELF_FILE_PATH = os.path.realpath(__file__)
SELF_DIR_PATH = os.path.dirname(SELF_FILE_PATH)
SELF_FILE = os.path.basename(SELF_FILE_PATH)


def align_read(args, sw_cfg, task_cfg):
    """
    Read Alignment
    
    :parm args: an argparse.Namespace object of main argumens {.log_file}
    :parm sw_cfg: a dictionary of corresponding software configureation {star|bwa}
    :parm task_cfg: a dictionary of corresponding task configureation {main, star|bwa, references}
    :returns: returns corresponding code snippets, in which will be written to shell_script_path if provided
    :raises keyError: NA
    """
    # logging
    if args.log_file is None:
        log_file_path = '{}.{}.log'.format(
            SELF_FILE_PATH,
            datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S"))
    else:
        log_file_path = args.log_file

    formatter = "%(asctime)-15s %(levelname)-8s %(message)s"
    logging.basicConfig(
        level=[logging.NOTSET, logging.DEBUG, logging.INFO,
               logging.WARNING, logging.ERROR, logging.CRITICAL][1],
        format=formatter,
        filename=log_file_path)

    for aligner in task_cfg.keys():
        if aligner in ['main', 'references']: continue
        
        logging.info("[ {} ] Aligner: {}\n".format(SELF_FILE, aligner))

        for library in task_cfg[aligner].keys():
            if library in ['main']: continue
            
            for read_group in task_cfg[aligner][library].keys():
                task_cfg_rg = {
                    'main': task_cfg['main'][library][read_group],
                    'references': task_cfg['references'],
                    aligner: task_cfg[aligner][library][read_group]
                }
    
                out_dir_t_path = task_cfg_rg[aligner]['out_dir_path']
                if not os.path.exists(out_dir_t_path):
                    logging.debug("MKDIR: [ $out_dir_t_path ]")
                    os.makedirs(out_dir_t_path, exist_ok=True)                    
                
                
                # util.ddictfunc.pprint_ddicts(task_cfg_rg); sys.exit()
                '''
                {'main': {'in_file_path_list': ['/group/bioinformatics/CRI_RNAseq_2018/example/data/KO01.test.fastq.gz']},
                 'references': {'grch38': {'anno_bed': None,
                                           'anno_bed12': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/gencode.v24.primary_assembly.annotation.chr11.bed12',
                                           'anno_gtf': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/gencode.v24.primary_assembly.annotation.chr11.gtf',
                                           'anno_refflat': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/gencode.v24.primary_assembly.annotation.chr11.refFlat.txt',
                                           'anno_ribosome_rna_bed': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/gencode.v24.primary_assembly.annotation.chr11.rRNA.bed',
                                           'anno_ribosome_rna_interval': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/gencode.v24.primary_assembly.annotation.chr11.rRNA.interval_list',
                                           'chrom_size': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/GRCh38.primary_assembly.genome.chr11.chrom.size',
                                           'chrs': 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM',
                                           'genome': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/GRCh38.primary_assembly.genome.chr11.fa',
                                           'genomedict': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/GRCh38.primary_assembly.genome.chr11.dict',
                                           'star_index': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11'}},
                 'star': {'log_file_path': '/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/Aln/star/KO01/SRR1205282/run.alignRead.star.SRR1205282.log',
                          'out_dir_path': '/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/Aln/star/KO01/SRR1205282',
                          'out_file_path_list': ['/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/Aln/star/KO01/SRR1205282/SRR1205282.star.bam'],
                          'shell_script_path': '/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/Shell/run.alignRead.star.SRR1205282.sh'}}
                '''

                sw_cfg_tool = util.ddictfunc.subset(sw_cfg, [aligner, 'sambamba'])

                if(aligner == "star"):
                    lib.star.star(args, sw_cfg_tool, task_cfg_rg)
                else:
                    logging.info("[ {} ] cannot recognize {}\n".format(SELF_FILE, aligner))                    

        logging.debug("[ {} ] Aligner: {} - DONE\n".format(SELF_FILE, aligner))

    return
