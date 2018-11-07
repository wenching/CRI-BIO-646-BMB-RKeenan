"""Loci Calling"""
import os
import sys
import logging
import datetime

import lib.deseq2
import lib.edger
import lib.limma

import util.ddictfunc

SELF_FILE_PATH = os.path.realpath(__file__)
SELF_DIR_PATH = os.path.dirname(SELF_FILE_PATH)
SELF_FILE = os.path.basename(SELF_FILE_PATH)


def call(args, sw_cfg, task_cfg):
    """
    Loci Calling
    
    :parm args: an argparse.Namespace object of main argumens {.log_file}
    :parm sw_cfg: a dictionary of corresponding software configureation {deseq2|edger|limma}
    :parm task_cfg: a dictionary of corresponding task configureation {main, deseq2|edger|limma, references}
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
    
    for caller in task_cfg.keys():
        if caller in ["main", "references"]: continue
        logging.info("[ {} ] Caller: {}".format(SELF_FILE, caller))
        
        for quantifier in task_cfg[caller].keys():
            if quantifier in ["main"]: continue
            
            logging.info("[ {} ] Qauntifier: {}\n".format(SELF_FILE, quantifier))
    
            for aligner in task_cfg[caller][quantifier].keys():
                if aligner in ['main']: continue
                
                logging.info("[ {} ] Aligner: {}".format(SELF_FILE, aligner))

                out_dir_t_path = task_cfg[caller][quantifier][aligner]['out_dir_path']

                if not os.path.exists(out_dir_t_path):
                    logging.debug("MKDIR: [ $out_dir_t_path ]")
                    os.makedirs(out_dir_t_path, exist_ok=True)


                task_cfg_proj = {
                    'references': task_cfg['references'],
                    caller: task_cfg['main'][quantifier][aligner]
                }    
                task_cfg_proj[caller].update(
                    util.ddictfunc.subset(task_cfg['main'], [quantifier], invert = True)
                )
                task_cfg_proj[caller].update(task_cfg[caller][quantifier][aligner])
                
                #util.ddictfunc.pprint_ddicts(task_cfg_proj)
                '''
                {'edger': {'application': 'RNAseq',
                           'comparison': {'KO_vs_WT': {'criteria': None,
                                                       'filter_4_low_expr': None,
                                                       'fold_change': None,
                                                       'pval': None}},
                           'criteria': 'fc_q',
                           'filter_4_low_expr': True,
                           'fold_change': 1.5,
                           'in_dir_path': '/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/Quantification/featurecounts/star',
                           'in_file_path_list': ['/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/Quantification/featurecounts/star/KO01/KO01.star.featurecounts.count',
                                                 '/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/Quantification/featurecounts/star/KO02/KO02.star.featurecounts.count',
                                                 '/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/Quantification/featurecounts/star/KO03/KO03.star.featurecounts.count',
                                                 '/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/Quantification/featurecounts/star/WT01/WT01.star.featurecounts.count',
                                                 '/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/Quantification/featurecounts/star/WT02/WT02.star.featurecounts.count',
                                                 '/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/Quantification/featurecounts/star/WT03/WT03.star.featurecounts.count'],
                           'log_file_path': '/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/DEG/edger/featurecounts/star/run.call.edger.featurecounts.star.DLBC.log',
                           'metadata_path': '/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC.metadata.txt',
                           'out_dir_path': '/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/DEG/edger/featurecounts/star',
                           'out_file_path_list': ['/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/DEG/edger/featurecounts/star/DLBC.star.featurecounts.edger.count.txt'],
                           'pval': 0.1,
                           'shell_script_path': '/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/shell_scripts/run.call.edger.featurecounts.star.DLBC.sh'},
                 'references': {'grch38': {'anno_bed': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/gencode.v24.primary_assembly.annotation.chr11.bed12',
                                           'anno_bed12': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/gencode.v24.primary_assembly.annotation.chr11.bed12',
                                           'anno_gtf': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/gencode.v24.primary_assembly.annotation.chr11.gtf',
                                           'anno_refflat': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/gencode.v24.primary_assembly.annotation.chr11.refFlat.txt',
                                           'anno_ribosome_rna_bed': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/gencode.v24.primary_assembly.annotation.chr11.rRNA.bed',
                                           'anno_ribosome_rna_interval': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/gencode.v24.primary_assembly.annotation.chr11.rRNA.interval_list',
                                           'chrom_size': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/GRCh38.primary_assembly.genome.chr11.chrom.size',
                                           'chrs': 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM',
                                           'genome': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/GRCh38.primary_assembly.genome.chr11.fa',
                                           'genomedict': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/GRCh38.primary_assembly.genome.chr11.dict',
                                           'star_index': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11'}}}
                '''

                sw_cfg_tool = util.ddictfunc.subset(sw_cfg, [caller])

                if(caller == "deseq2"):
                    lib.deseq2.deseq2(args, sw_cfg_tool, task_cfg_proj)
                elif(caller == "edger"):
                    lib.edger.edger(args, sw_cfg_tool, task_cfg_proj)
                elif(caller == "limma"):
                    lib.limma.limma(args, sw_cfg_tool, task_cfg_proj)

                logging.debug("[ {} ] Aligner: {} - DONE\n".format(SELF_FILE, aligner))
            logging.debug("[ {} ] Quantifier: {} - DONE\n".format(SELF_FILE, quantifier))
        logging.debug("[ {} ] Caller: {} - DONE\n".format(SELF_FILE, caller))

    return
