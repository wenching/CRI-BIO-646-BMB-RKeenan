"""Alignment Merging"""
import os
import sys
import logging
import datetime

import lib.sambamba

import util.ddictfunc

SELF_FILE_PATH = os.path.realpath(__file__)
SELF_DIR_PATH = os.path.dirname(SELF_FILE_PATH)
SELF_FILE = os.path.basename(SELF_FILE_PATH)


def merge_aln(args, sw_cfg, task_cfg):
    """
    Alignment Merging
    
    :parm args: an argparse.Namespace object of main argumens {.log_file}
    :parm sw_cfg: a dictionary of corresponding software configureation {sambamba}
    :parm task_cfg: a dictionary of corresponding task configureation {sambamba::star}
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
        logging.info("[ {} ] Aligner: {}\n".format(SELF_FILE, aligner))

        for library in task_cfg[aligner].keys():
            if library in ['main']: continue

            out_dir_t_path = task_cfg[aligner][library]['out_dir_path']
            if not os.path.exists(out_dir_t_path):
                logging.debug("MKDIR: [ $out_dir_t_path ]")
                os.makedirs(out_dir_t_path, exist_ok=True)            

            task_cfg_lb = {'sambamba': task_cfg[aligner][library]}
            
            # util.ddictfunc.pprint_ddicts(task_cfg_lb); sys.exit()
            '''
            {'sambamba': {'in_file_path_list': ['/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/Aln/star/KO01/SRR1205282/SRR1205282.star.bam'],
                          'log_file_path': '/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/Aln/star/KO01/run.mergeAln.star.KO01.log',
                          'out_dir_path': '/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/Aln/star/KO01',
                          'out_file_path_list': ['/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/Aln/star/KO01/KO01.star.bam'],
                          'shell_script_path': '/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/Shell/run.mergeAln.star.KO01.sh'}}
            
            '''

            lib.sambamba.merge(args, sw_cfg, task_cfg_lb)            

        logging.debug("[ {} ] Aligner: {} - DONE\n".format(SELF_FILE, aligner))
    return
