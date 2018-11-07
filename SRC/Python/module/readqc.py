"""Read QC"""
import os
import sys
import logging
import datetime

import lib.fastqc

import util.ddictfunc

SELF_FILE_PATH = os.path.realpath(__file__)
SELF_DIR_PATH = os.path.dirname(SELF_FILE_PATH)
SELF_FILE = os.path.basename(SELF_FILE_PATH)


def read_qc(args, sw_cfg, task_cfg):
    """
    Read QC
    
    :parm args: an argparse.Namespace object of main argumens {.log_file}
    :parm sw_cfg: a dictionary of corresponding software configureation {fastqc}
    :parm task_cfg: a dictionary of corresponding task configureation {fastqc::library}
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

    for library in task_cfg.keys():
        if library in ['main']: continue
        logging.info("[ {} ] Library: {}\n".format(SELF_FILE, library))

        for read_group in task_cfg[library].keys():
            logging.info("[ {} ] ReadGroup: {}".format(SELF_FILE, read_group))

            out_dir_t_path = task_cfg[library][read_group]['out_dir_path']
            if not os.path.exists(out_dir_t_path):
                logging.debug("MKDIR: [ $out_dir_t_path ]")
                os.makedirs(out_dir_t_path, exist_ok=True)
            
            task_cfg_rg =  {'fastqc': task_cfg[library][read_group]}

            # util.ddictfunc.pprint_ddicts(task_cfg_rg); sys.exit()
            '''
            {'fastqc': {'in_file_path_list': ['/group/bioinformatics/CRI_RNAseq_2018/example/data/KO01.test.fastq.gz'],
                        'log_file_path': '/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/RawReadQC/KO01/SRR1205282/run.RawReadQC.FastQC.SRR1205282.log',
                        'out_dir_path': '/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/RawReadQC/KO01/SRR1205282',
                        'out_file_path_list': ['/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/RawReadQC/KO01/SRR1205282/KO01.test_fastqc.zip'],
                        'shell_script_path': '/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/Shell/run.RawReadQC.FastQC.SRR1205282.sh'}}
            
            '''

            lib.fastqc.fastqc(args, sw_cfg, task_cfg_rg)
            
            logging.debug("[ {} ] ReadGroup: {} - DONE\n".format(SELF_FILE, read_group))             
    
        logging.debug("[ {} ] Library: {} - DONE\n".format(SELF_FILE, library))
    return
