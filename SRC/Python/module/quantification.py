"""Quantification"""
import os
import sys
import logging
import datetime

import lib.subread

import util.ddictfunc

SELF_FILE_PATH = os.path.realpath(__file__)
SELF_DIR_PATH = os.path.dirname(SELF_FILE_PATH)
SELF_FILE = os.path.basename(SELF_FILE_PATH)


def quant(args, sw_cfg, task_cfg):
    """
    Quantification
    
    :parm args: an argparse.Namespace object of main argumens {.log_file}
    :parm sw_cfg: a dictionary of corresponding software configureation {featurecounts}
    :parm task_cfg: a dictionary of corresponding task configureation {main, featurecounts, references}
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
    
    for quantifer in task_cfg.keys():
        if quantifer in ["main", "references"]: continue
        logging.info("[ {} ] Quantifier: {}".format(SELF_FILE, quantifer))
        
        for aligner in task_cfg[quantifer].keys():
            if aligner in ["main"]: continue
            
            logging.info("[ {} ] Aligner: {}\n".format(SELF_FILE, aligner))
    
            for library in task_cfg[quantifer][aligner].keys():
                if library in ['main']: continue
                
                logging.info("[ {} ] Library: {}".format(SELF_FILE, library))

                out_dir_t_path = task_cfg[quantifer][aligner][library]['out_dir_path']

                if not os.path.exists(out_dir_t_path):
                    logging.debug("MKDIR: [ $out_dir_t_path ]")
                    os.makedirs(out_dir_t_path, exist_ok=True)


                task_cfg_lb = {
                    'references': task_cfg['references'],
                    quantifer: task_cfg['main'][aligner][library]
                }                    
                task_cfg_lb[quantifer].update(task_cfg[quantifer][aligner][library])
                
                #util.ddictfunc.pprint_ddicts(task_cfg_lb)

                sw_cfg_tool = util.ddictfunc.subset(sw_cfg, [quantifer])

                if(quantifer == "featurecounts"):
                    lib.subread.featurecounts(args, sw_cfg_tool, task_cfg_lb)              

                logging.debug("[ {} ] Library: {} - DONE\n".format(SELF_FILE, library))
            logging.debug("[ {} ] Aligner: {} - DONE\n".format(SELF_FILE, aligner))
        logging.debug("[ {} ] Quantifier: {} - DONE\n".format(SELF_FILE, quantifer))

    return
