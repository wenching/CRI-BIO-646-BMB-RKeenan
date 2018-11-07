"""CRI-RNAseq-2018: Build RNAseq Project."""
import os
import sys
import re
import logging
import datetime
import argparse
from collections import defaultdict

import pandas
import yaml
#import IPython

import module.makemasterbds
import module.readqc
import module.alignread
import module.mergealn
import module.alnqc
import module.quantification
import module.callloci
import module.quantqc
import module.locistat
import module.postana

import util.ddictfunc
import util.hashlib


print("PYTHON_EXCUTABLE_VERSION = [{}]".format(os.path.basename(sys.executable)))


VERSION = "0.0.1"


SELF_FILE_PATH = os.path.realpath(__file__)
SELF_DIR_PATH = os.path.dirname(SELF_FILE_PATH)
SELF_FILE = os.path.basename(SELF_FILE_PATH)
print("SELF_DIR: {}".format(SELF_DIR_PATH))
print("SELF_PATH: {}".format(SELF_FILE_PATH))


def main(args = None):
    """CRI-RNAseq-2018: Build RNAseq Project."""

    parser = argparse.ArgumentParser(
      description='''
Program     : {}
Description : Build project directory and config files for running RNAseq pipeline
Version     : 2018-05-31 (release {})
License     : LGPLv3
Contact     : Wen-Ching Chan <wchan10@bsd.uchicago.edu>
  '''.format(SELF_FILE, VERSION),
     epilog='''
%(prog)s Version {}
--projdir /Volumes/bioinformatics/cri_rnaseq_2018/example
--metadata /Volumes/bioinformatics/cri_rnaseq_2018/example/DLBC.metadata.txt
--config /Volumes/bioinformatics/cri_rnaseq_2018/example/DLBC.pipeline.yaml
--systype cluster
--log_file /gpfs/data/bioinformatics/cri_rnaseq_2018/Build_RNAseq.DLBC.YYYY-MM-DD_HH_MM_SS.log
    '''.format(VERSION)
    )

    parser.add_argument(
        '-v', '--version', action='version', version='%(prog)s Version {}'.format(VERSION))
    parser.add_argument(
        '-l', '--log', default=sys.stdout, type=argparse.FileType('w'),
        help='Log all tasks (do not delete tmp files).')

    essential_argument = parser.add_argument_group('required named arguments')
    essential_argument.add_argument(
        '-j', '--projdir', dest='proj_dir',
        required=True, help='Project Directory.')
    essential_argument.add_argument(
        '-m', '--metadata', dest='metadata_path',
        required=True, help='Sample Metadata Table in tab-delimited format.')
    essential_argument.add_argument(
        '-c', '--config', dest='config_file',
        required=True, help='Pipeline Configuration File in YAML format.')
    essential_argument.add_argument(
        '-s', '--systype', dest='system_type',
        required=True, help='System Type.')
    essential_argument.add_argument(
        '-t', '--threads', dest='threads', type=int, default=1,
        required=False, help='Number of Threads.')
    essential_argument.add_argument(
        '-g', '--log_file', dest='log_file',
        required=False, help='Number of Threads.')

    # parser.print_help()
    print("PARSE ARGUMENTS...")
    args = parser.parse_args()
    # print(type(args)); sys.exit()

    # logging
    if args.log_file is None:
        log_file_path = '{}.{}.log'.format(
            SELF_FILE_PATH,
            datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S"))
    else:
        log_file_path = args.log_file
        if os.path.isfile(args.log_file):
            os.remove(args.log_file)


    formatter = "%(asctime)-15s %(levelname)-8s %(message)s"
    logging.basicConfig(
        level=[logging.NOTSET, logging.DEBUG, logging.INFO,
               logging.WARNING, logging.ERROR, logging.CRITICAL][1],
        format=formatter,
        filename=log_file_path)

    console = logging.StreamHandler()
    console.setLevel([logging.NOTSET, logging.DEBUG, logging.INFO,
                      logging.WARNING, logging.ERROR, logging.CRITICAL][2])
    console.setFormatter(logging.Formatter(formatter))
    logging.getLogger('').addHandler(console)


    # Build_RNAseq.pl > GetOpt.pm
    logging.info("METADATA:        = [ {} ]".format(args.metadata_path))
    logging.info("CONFIG:          = [ {} ]".format(args.config_file))
    logging.info("PROJECT_RESULT:  = [ {} ]".format(args.proj_dir))
    args.proj = os.path.basename(args.proj_dir)
    logging.info("PROJECT_ID:      = [ {} ]".format(args.proj))
    logging.info("SYSTEM_TYPE:     = [ {} ]".format(args.system_type))
    logging.info("REQUIRED_NODE:   = [ {} ]".format(args.threads))
    args.job_script = "{}/Submit_{}.sh".format(args.proj_dir, args.proj)
    args.submitter = "{}/Submit_{}.bds".format(args.proj_dir, args.proj)
    logging.info("SUBMISSION_SCRIPT: = [ {} ]".format(args.job_script))
    args.submitter_cfg = "{}/{}.cfg".format(args.proj_dir, re.sub(".bds$", "", os.path.basename(args.submitter)))
    logging.info("SUBMISSION_CFG: = [ {} ]".format(args.submitter_cfg))
    args.tree_output = "{}/{}.tree.txt".format(args.proj_dir, args.proj)
    logging.info("PROJECT_DIR_TREE: = [ {} ]".format(args.tree_output))

    out_dir_t_path = args.proj_dir
    if not os.path.exists(out_dir_t_path):
        logging.debug("MKDIR: [ {} ]".format(out_dir_t_path))
        os.makedirs(out_dir_t_path, exist_ok=True)


    # Build_RNAseq.pl > BuildProject.pm > ParseMetadata.pm
    logging.info("LOAD Sample Information\n")

    if not os.path.isfile(args.metadata_path):
        sys.exit("cannot find {}".format(args.metadata_path))
    else:
        curr_metadata_md5 = util.hashlib.md5(args.metadata_path)

    # Build_RNAseq.pl > BuildProject.pm > ParseMetadata.pm > SampleList.pm
    df_metadata = pandas.read_csv(args.metadata_path, sep="\t", comment='#')
    logging.debug(type(df_metadata.index)) # pandas.core.frame.DataFrame
    logging.debug(list(map(type, list(df_metadata))))
    logging.debug(df_metadata['Sample'])

    logging.debug("ROW_NAMES: [ {} ]".format(' | '.join(str(v) for v in df_metadata.index)))
    logging.debug("COL_NAMES: [ {} ]".format(' | '.join(df_metadata.columns)))
    logging.debug(df_metadata[['Sample', 'ReadGroup']])
    logging.debug(df_metadata.loc[:, ['Sample', 'ReadGroup']]) # ~ df_metadata.loc[df_metadata.index, ['Sample', 'ReadGroup']]
    logging.debug(df_metadata.loc[0:2, ['Sample', 'ReadGroup']])


    logging.info("VERIFIY of Sample Information\n")

    genome_list = df_metadata['Genome'].drop_duplicates().values.tolist()

    if(len(genome_list) > 1):
        logging.error(
            "Now only allow one genome per sample list, but multiple were provided [ {} ]".format(
                ', '.join(genome_list)
            )
        )
    else:
        target_genome = genome_list[0]
        logging.info(
            "TARGET GENOME: [ {} ]".format(
                target_genome
            )
        )

    df_lib_type_count_per_lib = \
        df_metadata.loc[:,['Library', 'LibType']].drop_duplicates().groupby(by = "Library").agg({'LibType': 'count'})

    if((df_lib_type_count_per_lib != 1).any(axis = {'column-wise': 0, 'row-wise': 1}['column-wise']).bool()):
        logging.error(
            "Now only allow one strand specificity per library, but multiple were provided [ {} ]".format(
                ', '.join(df_lib_type_count_per_lib.loc[:, (df_lib_type_count_per_lib != 1).any(axis = 0)].index)
            )
        )

    df_Flavor_count_per_lib = \
        df_metadata.loc[:,['Library', 'Flavor']].drop_duplicates().groupby(by = "Library").agg({'Flavor': 'count'})

    if((df_Flavor_count_per_lib != 1).any(axis = {'column-wise': 0, 'row-wise': 1}['column-wise']).bool()):
        logging.error(
            "Now only allow one Flavor (i.e., SEQUENCING_LAYOUTxSEQUENCING_LENGTH, in which SEQUENCING_LAYOUT = 1/2; e.g., 1x49) per library, but multiple were provided [ {} ]".format(
                ', '.join(df_Flavor_count_per_lib.loc[:, (df_Flavor_count_per_lib != 1).any(axis = 0)].index)
            )
        )

    if not os.path.isabs(df_metadata.iloc[0]['Location']):
        logging.info("COMPLETE relative path to absolute path")
        df_metadata['Location'] = \
            args.proj_dir + "/" + df_metadata['Location'].astype(str)
    

    logging.debug("VERIFIY of Sample Information - DONE\n")


    logging.debug("LOAD Sample Information - DONE\n")


    # Build_RNAseq.pl > BuildProject.pm > YAML::Tiny::LoadFile
    logging.info("LOAD Project Configuration\n")

    if not os.path.isfile(args.config_file):
        sys.exit("cannot find {}".format(args.config_file))

    with open(args.config_file) as stream:
        try:
            project_cfg = yaml.load(stream)
            logging.debug("CFG_IN_YAML:")
            logging.debug(project_cfg)
            logging.debug(type(project_cfg))
        except yaml.YAMLError as exc:
            print(exc)


    project_cfg['project']['name'] = args.proj

    logging.info("DETERMINE whether running as practice\n")

    if project_cfg['project']['ex_metadata_md5'] == curr_metadata_md5:
        logging.info("RUN AS PRACTICE\n")
        project_cfg['project']['run_as_practice'] = True
    else:
        project_cfg['project']['run_as_practice'] = False

    logging.debug("DETERMINE whether running as practice - DONE\n")

    logging.debug("LOAD Project Configuration - DONE\n")


    #util.ddictfunc.pprint_ddicts(project_cfg)

    logging.info("PROCESS Project Configuration\n")

    # Build_RNAseq.pl > BuildProject.pm > BuildSampleCfg.pm > PrintSample.pm > Util.pm::Add_Key()
    project_cfg['aligners'] = []
    for k, v in project_cfg['pipeline']['flags']['aligners'].items():
        if v:
            project_cfg['aligners'].append(re.sub(r"""^run_""", "", k))
    # util.ddictfunc.pprint_ddicts(project_cfg['aligners'])

    # Build_RNAseq.pl > BuildProject.pm > BuildSampleCfg.pm > PrintSample.pm > Util.pm::Add_Key()
    project_cfg['quantifiers'] = []
    for k, v in project_cfg['pipeline']['flags']['quantifiers'].items():
        if v:
            project_cfg['quantifiers'].append(re.sub(r"""^run_""", "", k))
    # util.ddictfunc.pprint_ddicts(project_cfg['quantifiers'])

    project_cfg['callers'] = []
    for k, v in project_cfg['pipeline']['flags']['callers'].items():
        if v:
            project_cfg['callers'].append(re.sub(r"""^run_""", "", k))
    # util.ddictfunc.pprint_ddicts(project_cfg['callers'])

    if target_genome not in project_cfg['pipeline']['references']:
        logging.error(
            "cannot find {} in configuration file with [ {} ]".format(
                target_genome,
                ', '.join(project_cfg['pipeline']['references'].keys())
            )
        )
    else:
        project_cfg['pipeline']['references']['target'] = target_genome

    # util.ddictfunc.pprint_ddicts(project_cfg['pipeline']['references'])

    project_cfg['project']['proj_dir_path'] = args.proj_dir

    project_cfg['project']['appl_dir_path'] = \
        '''{}/{}'''.format(
            project_cfg['project']['proj_dir_path'],
            project_cfg['project']['application']
        )
    out_dir_t_path = project_cfg['project']['appl_dir_path']
    if not os.path.exists(out_dir_t_path):
        logging.debug("MKDIR: [ {} ]".format(out_dir_t_path))
        os.makedirs(out_dir_t_path, exist_ok=True)

    project_cfg['project']['shell_dir_path'] = \
        '{}/shell_scripts'.format(
            project_cfg['project']['appl_dir_path']
        )
    out_dir_t_path = project_cfg['project']['shell_dir_path']
    logging.info("SHELL_DIR: [ {} ]\n".format(out_dir_t_path))
    if not os.path.exists(out_dir_t_path):
        logging.debug("MKDIR: [ {} ]".format(out_dir_t_path))
        os.makedirs(out_dir_t_path, exist_ok=True)

    project_cfg['project']['metadata_path'] = args.metadata_path
    project_cfg['pipeline']['software']['main']['r_dir_path'] = \
        '''{}/{}'''.format(
            os.path.dirname(os.path.dirname(SELF_DIR_PATH)),
            project_cfg['pipeline']['software']['main']['r_dir_path']
        )

    for script_tool in project_cfg['pipeline']['software']:
        if 'script' in project_cfg['pipeline']['software'][script_tool] :
            if re.search("\.R$", project_cfg['pipeline']['software'][script_tool]['script']) is not None:
                project_cfg['pipeline']['software'][script_tool]['script'] = \
                    '''{}/{}'''.format(
                        project_cfg['pipeline']['software']['main']['r_dir_path'],
                        project_cfg['pipeline']['software'][script_tool]['script']
                    )

    logging.debug("PROCESS Project Configuration - DONE\n")


    logging.info("Make The Master BigDataScript\n")

    task_cfg = util.ddictfunc.ddicts()
    module.makemasterbds.make_master_bds(
        args=args,
        metadata = df_metadata,
        project_cfg = project_cfg,
        task_cfg = task_cfg
    )

    task_cfg = util.ddictfunc.ddicts_2_dict(task_cfg)

    logging.debug("Make The Master BigDataScript - DONE\n")


    logging.info("Raw Read QC\n")

    # util.ddictfunc.pprint_ddicts(task_cfg, ['read_qc'])

    module.readqc.read_qc(
        args=args,
        sw_cfg=util.ddictfunc.subset(project_cfg['pipeline']['software'], ["fastqc"]),
        task_cfg=task_cfg['read_qc']
    )

    logging.debug("Raw Read QC - DONE\n")


    logging.info("Read Alignment\n")

    # util.ddictfunc.pprint_ddicts(task_cfg, ['align_read'])

    module.alignread.align_read(
        args=args,
        sw_cfg=util.ddictfunc.subset(project_cfg['pipeline']['software'], project_cfg['aligners'] + ["sambamba"]),
        task_cfg=task_cfg['align_read']
    )

    logging.debug("Read Alignment - DONE\n")


    logging.info("Alignment Merging\n")

    module.mergealn.merge_aln(
        args=args,
        sw_cfg=util.ddictfunc.subset(project_cfg['pipeline']['software'], ['sambamba']),
        task_cfg=task_cfg['merge_aln']
    )

    logging.debug("Alignment Merging - DONE\n")


    logging.info("Alignment QC\n")

    if("RNAseq" in project_cfg['project']['application']):
        dict_aln_qc_tool = \
            {
                "picard":
                [
                    "CollectRnaSeqMetrics"
                ],
                "rseqc":
                [
                    "clipping_profile.py",
                    "geneBody_coverage.py", # DISABLE at makemasterbds.py due to long running time
                    "RPKM_saturation.py", # DISABLE at makemasterbds.py due to long running time
                    "infer_experiment.py"
                ]
            }
    else:
        logging.error(
            "PLEASE specify tools for alignment QC in [{}]\n".format(
                project_cfg['project']['application']
            )
        )
        sys.exit()

    #util.ddictfunc.pprint_ddicts(project_cfg['pipeline']['software'], dict_aln_qc_tool.keys())

    module.alnqc.aln_qc(
        args=args,
        sw_cfg=util.ddictfunc.subset(project_cfg['pipeline']['software'], dict_aln_qc_tool.keys()),
        task_cfg=task_cfg['aln_qc']
    )

    logging.debug("Alignment QC - DONE\n")


    logging.info("Quantification\n")

    #util.ddictfunc.pprint_ddicts(project_cfg['pipeline']['software'], project_cfg['quantifiers'])

    module.quantification.quant(
        args=args,
        sw_cfg=util.ddictfunc.subset(project_cfg['pipeline']['software'], project_cfg['quantifiers']),
        task_cfg=task_cfg['quant']
    )

    logging.debug("Quantification - DONE\n")


    logging.info("Loci Calling\n")

    #util.ddictfunc.pprint_ddicts(project_cfg['pipeline']['software'], project_cfg['callers'])

    module.callloci.call(
        args=args,
        sw_cfg=util.ddictfunc.subset(project_cfg['pipeline']['software'], project_cfg['callers']),
        task_cfg=task_cfg['call']
    )

    logging.debug("Loci Calling - DONE\n")


    logging.info("Quantification QC\n")

    if("RNAseq" in project_cfg['project']['application']):
        quantqctool_list = \
            [
                "pca"
            ]
    else:
        logging.error(
            "PLEASE specify tools for quantification QC in [{}]\n".format(
                project_cfg['project']['application']
            )
        )
        sys.exit()

    #util.ddictfunc.pprint_ddicts(project_cfg['pipeline']['software'], quantqctool_list)

    module.quantqc.quant_qc(
        args=args,
        sw_cfg=util.ddictfunc.subset(project_cfg['pipeline']['software'], quantqctool_list),
        task_cfg=task_cfg['quant_qc']
    )

    logging.debug("Quantification QC - DONE\n")


    logging.info("Loci Statistics\n")

    if("RNAseq" in project_cfg['project']['application']):
        locistattool_list = \
            [
                "venn"
            ]
    else:
        logging.error(
            "PLEASE specify tools for loci statistics in [{}]\n".format(
                project_cfg['project']['application']
            )
        )
        sys.exit()

    #util.ddictfunc.pprint_ddicts(project_cfg['pipeline']['software'], locistattool_list)

    module.locistat.loci_stat(
        args=args,
        sw_cfg=util.ddictfunc.subset(project_cfg['pipeline']['software'], locistattool_list),
        task_cfg=task_cfg['loci_stat']
    )

    logging.debug("Loci Statistics - DONE\n")


    logging.info("Post Analysis - Heat Map & GSEA\n")

    if("RNAseq" in project_cfg['project']['application']):
        dict_post_ana_tool = \
            {
                "pheatmap":
                [
                    "heatmap.pdf"
                ],
                "clusterprofiler":
                [
                    "enrichGO.ALL.txt"
                ]
            }
    else:
        logging.error(
            "PLEASE specify tools for post analysis in [{}]\n".format(
                project_cfg['project']['application']
            )
        )
        sys.exit()

    #util.ddictfunc.pprint_ddicts(project_cfg['pipeline']['software'], dict_post_ana_tool.keys())
    #util.ddictfunc.pprint_ddicts(task_cfg['post_ana'])

    module.postana.post_ana(
        args=args,
        sw_cfg=util.ddictfunc.subset(project_cfg['pipeline']['software'], dict_post_ana_tool.keys()),
        task_cfg=task_cfg['post_ana']
    )

    logging.debug("Post Analysis - Heat Map & GSEA - DONE\n")


    logging.info("Make Submitter Shell Script of Master BigDataScript Script\n")

    # Build_RNAseq.pl > BuildProject.pm > BuildPipelineJob.pm::Print_Job()
    try:
        with open(args.job_script, "w") as outfile:
            outfile.write('#!/bin/bash\n')

            outfile.write('\n\n\n## set up environment\n')
            # Build_RNAseq.pl > BuildProject.pm > BuildPipelineJob.pm::Print_Job() > SetBDSenv.pm
            outfile.write('\n'.join(list()))

            outfile.write('\n## submit job \n')
            outfile.write('''\n

module purge; module load gcc/6.2.0 java-jdk/1.8.0_92 bds; module update

bds \\
-config {} \\
-reportHtml \\
-retry 0 \\
-verbose \\
-log \\
-system {} \\
{} \\
&> {}/{}.log &

    \n'''.format(
             args.submitter_cfg,
             args.system_type,
             args.submitter,
             args.proj_dir,
             os.path.basename(args.submitter)
         ))
            outfile.close()
    except IOError as exc:
        print(exc)

    logging.debug("Make Submitter Shell Script of Master BigDataScript Script - DONE\n")


    logging.info("Make Configuration File of Submitter Shell Script\n")

    # Build_RNAseq.pl > BuildProject.pm > BuildPipelineJob.pm > PrintBDScfg.pm
    # Build_RNAseq.pl > BuildProject.pm > BuildPipelineJob.pm > PrintBDScfg.pm::BDS_Cfg()
    try:
        with open(args.submitter_cfg, "w") as outfile:
            outfile.write('''\n
####--------------------------------------------------####
####                                                  ####
####         BigDataScript configuration file         ####
####                                                  ####
####             BDS documentation website            ####
####      http://pcingola.github.io/BigDataScript     ####
####                                                  ####
####--------------------------------------------------####

#---
# Mesos parameters
#---

#mesos.master = 127.0.0.1:5050

##---
## Default parameters
##---

## Default number of retries
retry = 0

## Wait time in between job submission (milli seconds)
waitAfterTaskRun = 1000

## Set task shell and sys shell env
taskShell = /bin/bash -e
sysShell = /bin/bash -e -c

## Default memory (-1 = unrestricted)
# mem = -1

## Default execution node (none)
# node = ""

## Add default queue name (if any)
# queue = ""

## Task timeout in seconds (default is one day)
# timeout = 86400

##---
## SGE parameters
##---

## Parallel environment in SGE (e.g. 'qsub -pe mpi 4')
## Custom CRI-openmp sge parallel environment added (allocation_rule $pe_slots)
# sge.pe = orte
sge.pe = smp

## Parameter for requesting amount of memory in qsub (e.g. 'qsub -l mem 4G')
## Note on sge, mem_free is per slot!
sge.mem = mem_free

## Parameter for timeout in qsub (e.g. 'qsub -l h_rt 24:00:00')
sge.timeout = h_rt
    \n''')
            outfile.close()
    except IOError as exc:
        print(exc)

    logging.debug("Make Configuration File of Submitter Shell Script - DONE\n")


    logging.info("CURRENT WORKING DIRECTORY: " + os.getcwd())

    #print(IPython.sys_info())

    sys.exit(1)


if __name__ == "__main__":
    # execute only if run as a script
    main()
