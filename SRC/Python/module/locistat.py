"""Loci Statistics"""
import os
import sys
import logging
import datetime

import lib.venn

import util.ddictfunc

SELF_FILE_PATH = os.path.realpath(__file__)
SELF_DIR_PATH = os.path.dirname(SELF_FILE_PATH)
SELF_FILE = os.path.basename(SELF_FILE_PATH)


def loci_stat(args, sw_cfg, task_cfg):
    """
    Loci Statistics
    
    :parm args: an argparse.Namespace object of main argumens {.log_file}
    :parm sw_cfg: a dictionary of corresponding software configureation {locistattool}
    :parm task_cfg: a dictionary of corresponding task configureation {main, featurecounts}
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

    for quantifier in task_cfg.keys():
        if quantifier in ["main", "references"]: continue
        logging.info("[ {} ] Quantifier: {}".format(SELF_FILE, quantifier))
        
        for aligner in task_cfg[quantifier].keys():
            if aligner in ["main"]: continue
            logging.info("[ {} ] Aligner: {}\n".format(SELF_FILE, aligner))

            for proj_comp in task_cfg[quantifier][aligner].keys():
                if proj_comp in ["main"]: continue
                logging.info("[ {} ] Comparison: {}".format(SELF_FILE, proj_comp))                                

                out_dir_t_path = task_cfg[quantifier][aligner][proj_comp]['out_dir_path']
                if not os.path.exists(out_dir_t_path):
                    logging.debug("MKDIR: [ $out_dir_t_path ]")
                    os.makedirs(out_dir_t_path, exist_ok=True)
    
                task_cfg_comp = {
                    proj_comp: task_cfg[quantifier][aligner][proj_comp]
                }                    
                task_cfg_comp[proj_comp].update(
                    util.ddictfunc.subset(task_cfg[quantifier][aligner]['main'], ['anchor_file_path'])
                )
                
                
                #util.ddictfunc.pprint_ddicts(task_cfg_comp)

                if(task_cfg['main']['application'] == "RNAseq"):
                    sw_cfg_tool = util.ddictfunc.subset(sw_cfg, ['venn'])
                    lib.venn.venn(args, sw_cfg_tool, task_cfg_comp)              

            
                logging.debug("[ {} ] Comparison: {} - DONE\n".format(SELF_FILE, proj_comp))                
            logging.debug("[ {} ] Aligner: {} - DONE\n".format(SELF_FILE, aligner))
        logging.debug("[ {} ] Quantifier: {} - DONE\n".format(SELF_FILE, quantifier))

    return
