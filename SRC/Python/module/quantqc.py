"""Quantification QC"""
import os
import sys
import logging
import datetime

import lib.pca

import util.ddictfunc

SELF_FILE_PATH = os.path.realpath(__file__)
SELF_DIR_PATH = os.path.dirname(SELF_FILE_PATH)
SELF_FILE = os.path.basename(SELF_FILE_PATH)


def quant_qc(args, sw_cfg, task_cfg):
    """
    Quantification QC
    
    :parm args: an argparse.Namespace object of main argumens {.log_file}
    :parm sw_cfg: a dictionary of corresponding software configureation {tool}
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
    
    for tool in task_cfg.keys():
        if tool in ["main", "references"]: continue
        logging.info("[ {} ] Quantification QC tool: {}".format(SELF_FILE, tool))    
        
        for quantifier in task_cfg[tool].keys():
            if quantifier in ["main"]: continue
            logging.info("[ {} ] Quantifier: {}".format(SELF_FILE, quantifier))
            
            for aligner in task_cfg[tool][quantifier].keys():
                if aligner in ["main"]: continue
                
                logging.info("[ {} ] Aligner: {}\n".format(SELF_FILE, aligner))

                out_dir_t_path = task_cfg['main'][quantifier][aligner]['out_dir_path']
                if not os.path.exists(out_dir_t_path):
                    logging.debug("MKDIR: [ $out_dir_t_path ]")
                    os.makedirs(out_dir_t_path, exist_ok=True)
    
                task_cfg_proj = {
                    tool: task_cfg[tool][quantifier][aligner]
                }                    
                task_cfg_proj[tool].update(
                    task_cfg['main'][quantifier][aligner]
                )
                task_cfg_proj[tool].update(
                    {
                        'library_info': task_cfg['main']['library_info']
                    }
                )
                
                
                #util.ddictfunc.pprint_ddicts(task_cfg_proj)
    
                sw_cfg_tool = util.ddictfunc.subset(sw_cfg, [tool])
    
                if(task_cfg['main']['application'] == "RNAseq"):
                    if(tool == "pca"):
                        lib.pca.pca(args, sw_cfg_tool, task_cfg_proj)              
    
                logging.debug("[ {} ] Aligner: {} - DONE\n".format(SELF_FILE, aligner))
            logging.debug("[ {} ] Quantifier: {} - DONE\n".format(SELF_FILE, quantifier))
        logging.debug("[ {} ] Quantification QC tool: {} - DONE\n".format(SELF_FILE, tool))

    return
