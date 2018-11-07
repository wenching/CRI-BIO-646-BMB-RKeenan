"""Make Master BigDataScript."""
import os
import sys
import re
import logging
import datetime
import itertools

import util.ddictfunc


SELF_FILE_PATH = os.path.realpath(__file__)
SELF_DIR_PATH = os.path.dirname(SELF_FILE_PATH)
SELF_FILE = os.path.basename(SELF_FILE_PATH)


lib_type_2_strand_specificity_dict = \
    {
        "NS": {
            "picard": "NONE",
            "rseqc__RPKM_saturation.py": None,
            "featurecounts": 0,
            "kallisto": ""
        },
        "RF": {
            "picard": "SECOND_READ_TRANSCRIPTION_STRAND",
            "rseqc__RPKM_saturation.py": "1+-,1-+,2++,2--",
            "featurecounts": 2,
            "kallisto": "--rf-stranded"
        },
        "FR": {
            "picard": "FIRST_READ_TRANSCRIPTION_STRAND",
            "rseqc__RPKM_saturation.py": "1++,1--,2+-,2-+",
            "featurecounts": 1,
            "kallisto": "--fr-stranded"
        }
    }

flavor_2_sequencing_layout_dict = \
    {
        "1": {
            "rseqc__clipping_profile.py": "SE",
            "featurecounts": False # featurecounts::isPairedEnd
        },
        "2": {
            "rseqc__clipping_profile.py": "PE",
            "featurecounts": True # featurecounts::isPairedEnd
        }
    }


def make_master_bds(args, metadata, project_cfg, task_cfg):
    """
    Make Master BigDataScript.

    :parm args: an argparse.Namespace object of main argumens {.log_file}
    :parm metadata: an pandas.core.frame.DataFrame object of sample information {DF::metadata}
    :parm project_cfg: a dictionary of project configureation {project_configuration}
    :parm task_cfg: a dictionary of task configureation {task_configuration}
    :returns: None
    :raises keyError: None
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

    out_file_path = args.submitter
    try:
        with open(out_file_path, "w") as outfile:
            outfile.write('''#!/usr/bin/env bds\n''')
            outfile.close()
    except IOError as exc:
        print(exc)


    logging.info("TRANSFORM Sample Information\n")

    ddicts_metadata = util.ddictfunc.ddicts()
    dict_library_info = {}

    for index, row in metadata.iterrows():
        ddicts_metadata[row['Group']][row['Sample']][row['Library']][row['ReadGroup']] = \
            util.ddictfunc.subset(row.to_dict(), ["Group", "Sample", "Library", "ReadGroup"], invert = True)
        dict_library_info[row['Library']] = row['Group']

    ddicts_metadata = util.ddictfunc.ddicts_2_dict(ddicts_metadata)
    #util.ddictfunc.pprint_ddicts(ddicts_metadata)
    #util.ddictfunc.pprint_ddicts(dict_library_info)

    logging.debug("TRANSFORM Sample Information - DONE\n")


    target_genome = project_cfg['pipeline']['references']['target']
    logging.debug(
        "TARGET GENOME: [ {} ]".format(
            target_genome
        )
    )

    proj_appl = project_cfg['project']['application']
    appl_dir_path = project_cfg['project']['appl_dir_path']
    shell_dir_path = project_cfg['project']['shell_dir_path']



    logging.info("[ {} ] Raw Read QC\n".format(SELF_FILE))

    out_file_path = args.submitter
    try:
        with open(out_file_path, "a") as outfile:
            outfile.write("\n// Raw Read QC >>>\n\n\n")
            outfile.close()
    except IOError as exc:
        print(exc)

    task_cfg['read_qc']['main']['out_step_dir'] = \
        '{}/RawReadQC'.format(
            appl_dir_path
        )
    out_dir_t_path = task_cfg['read_qc']['main']['out_step_dir']
    if not os.path.exists(out_dir_t_path):
        logging.debug("MKDIR: [ {} ]".format(out_dir_t_path))
        os.makedirs(out_dir_t_path, exist_ok=True)

    for group in ddicts_metadata.keys():
        for sample in ddicts_metadata[group].keys():
            for library in ddicts_metadata[group][sample].keys():
                for read_group in ddicts_metadata[group][sample][library].keys():
                    task_cfg['read_qc'][library][read_group]['out_dir_path'] = \
                        '{}/{}/{}'.format(
                            task_cfg['read_qc']['main']['out_step_dir'],
                            library,
                            read_group)

                    task_cfg['read_qc'][library][read_group]['in_file_path_list'] = \
                        [
                            '{}/{}'.format(
                                ddicts_metadata[group][sample][library][read_group]['Location'],
                                ddicts_metadata[group][sample][library][read_group]['Seqfile1']
                            )
                        ]

                    task_cfg['read_qc'][library][read_group]['out_file_path_list']= \
                        [
                            '{}/{}_fastqc.zip'.format(
                                task_cfg['read_qc'][library][read_group]['out_dir_path'],
                                re.sub(
                                    r'''.fastq.gz$''',
                                    "",
                                    ddicts_metadata[group][sample][library][read_group]['Seqfile1']
                                )
                            )
                        ]

                    if 'Seqfile2' in row:
                        task_cfg['read_qc'][library][read_group]['in_file_path_list'].append(
                            '{}/{}'.format(
                                ddicts_metadata[group][sample][library][read_group]['Location'],
                                ddicts_metadata[group][sample][library][read_group]['Seqfile2']
                            )
                        )
                        task_cfg['read_qc'][library][read_group]['out_file_path_list'].append(
                            '{}/{}_fastqc.zip'.format(
                                task_cfg['read_qc'][library][read_group]['out_dir_path'],
                                re.sub(
                                    r'''.fastq.gz$''',
                                    "",
                                    ddicts_metadata[group][sample][library][read_group]['Seqfile2']
                                )
                            )
                        )

                    task_cfg['read_qc'][library][read_group]['shell_script_path'] = \
                        '{}/run.RawReadQC.FastQC.{}.sh'.format(shell_dir_path, read_group)

                    task_cfg['read_qc'][library][read_group]['log_file_path'] = \
                        '{}/run.RawReadQC.FastQC.{}.log'.format(
                            task_cfg['read_qc'][library][read_group]['out_dir_path'],
                            read_group)

                    #util.ddictfunc.pprint_ddicts(task_cfg, ['read_qc'])
                    '''
                    {'read_qc': {'KO01': {'SRR1205282': {'in_file_path_list': ['/gpfs/data/bioinformatics/cri_rnaseq_2018/example/data/WT01.test.fastq.gz'],
                                                         'log_file_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/RawReadQC/KO01/SRR1205282/run.RawReadQC.FastQC.SRR1205282.log',
                                                         'out_dir_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/RawReadQC/KO01/SRR1205282',
                                                         'out_file_path_list': ['/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/RawReadQC/KO01/SRR1205282/SRR1205282_fastqc.zip'],
                                                         'shell_script_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/Shell/run.RawReadQC.FastQC.SRR1205282.sh'}},
                                 'main': {'out_step_dir': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/RawReadQC'}}}
                    '''

                    local_resource = ''
                    if args.system_type == 'cluster':
                        pbs_ppn = min([len(task_cfg['read_qc'][library][read_group]['in_file_path_list']), args.threads, 6])
                        local_resource = ', cpus := {}, mem := {}*G, timeout := {}*hour, taskName := "FastQC.{}"'.format(
                            pbs_ppn,
                            pbs_ppn * 16,
                            72,
                            read_group
                        )

                    out_file_path = args.submitter
                    try:
                        with open(out_file_path, "a") as outfile:
                            outfile.write('''
// Raw Read QC: [ {} ]

dep( [ '{}' ] <- [ '{}' ]{}) sys bash {}; sleep {}
goal( [ '{}' ] )
\n'''.format(
         read_group,
         "', '".join(task_cfg['read_qc'][library][read_group]['out_file_path_list']),
         "', '".join(task_cfg['read_qc'][library][read_group]['in_file_path_list']),
         local_resource,
         task_cfg['read_qc'][library][read_group]['shell_script_path'],
         project_cfg['pipeline']['software']['bigdatascript']['safeSleep'],
         "', '".join(task_cfg['read_qc'][library][read_group]['out_file_path_list'])
     ))
                            outfile.close()
                    except IOError as exc:
                        print(exc)

    out_file_path = args.submitter
    try:
        with open(out_file_path, "a") as outfile:
            outfile.write("\n// <<< Raw Read QC\n\n\n")
            outfile.close()
    except IOError as exc:
        print(exc)

    #util.ddictfunc.pprint_ddicts(task_cfg, ['read_qc'])

    logging.debug("[ {} ] Raw Read QC - DONE\n".format(SELF_FILE))



    logging.info("[ {} ] Read Alignment\n".format(SELF_FILE))

    out_file_path = args.submitter
    try:
        with open(out_file_path, "a") as outfile:
            outfile.write("\n// Read Alignment >>>\n\n\n")
            outfile.close()
    except IOError as exc:
        print(exc)

    for group in ddicts_metadata.keys():
        for sample in ddicts_metadata[group].keys():
            for library in ddicts_metadata[group][sample].keys():
                for read_group in ddicts_metadata[group][sample][library].keys():
                    task_cfg['align_read']['main'][library][read_group]['in_file_path_list'] = \
                        task_cfg['read_qc'][library][read_group]['in_file_path_list'].copy()

    for aligner in project_cfg['aligners']:
        logging.info("[ {} ] Read Alignment:: Aligner: {}\n".format(SELF_FILE, aligner))

        task_cfg['align_read'][aligner]['main']['out_step_dir'] = \
            '{}/Aln/{}'.format(
                appl_dir_path,
                aligner
            )
        out_dir_t_path = task_cfg['align_read'][aligner]['main']['out_step_dir']
        if not os.path.exists(out_dir_t_path):
            logging.debug("MKDIR: [ {} ]".format(out_dir_t_path))
            os.makedirs(out_dir_t_path, exist_ok=True)

        task_cfg['align_read']['references'][target_genome] = project_cfg['pipeline']['references'][target_genome]

        for group in ddicts_metadata.keys():
            for sample in ddicts_metadata[group].keys():
                for library in ddicts_metadata[group][sample].keys():
                    for read_group in ddicts_metadata[group][sample][library].keys():
                        task_cfg['align_read'][aligner][library][read_group]['out_dir_path'] = \
                            '{}/{}/{}'.format(
                                task_cfg['align_read'][aligner]['main']['out_step_dir'],
                                library,
                                read_group)
                        out_dir_t_path = task_cfg['align_read'][aligner][library][read_group]['out_dir_path']
                        if not os.path.exists(out_dir_t_path):
                            logging.debug("MKDIR: [ {} ]".format(out_dir_t_path))
                            os.makedirs(out_dir_t_path, exist_ok=True)

                        task_cfg['align_read'][aligner][library][read_group]['out_file_path_list'] = \
                            [
                                '{}/{}.{}.bam'.format(
                                    task_cfg['align_read'][aligner][library][read_group]['out_dir_path'],
                                    read_group,
                                    aligner
                                )
                            ]

                        task_cfg['align_read'][aligner][library][read_group]['shell_script_path'] = \
                            '{}/run.alignRead.{}.{}.sh'.format(shell_dir_path, aligner, read_group)

                        task_cfg['align_read'][aligner][library][read_group]['log_file_path'] = \
                            '{}/run.alignRead.{}.{}.log'.format(
                                task_cfg['align_read'][aligner][library][read_group]['out_dir_path'],
                                aligner,
                                read_group)

                        #util.ddictfunc.pprint_ddicts(task_cfg, ['align_read'])
                        '''
                        {'align_read': {'main': {'KO01': {'SRR1205282': {'in_file_path_list': ['/group/bioinformatics/CRI_RNAseq_2018/example/data/KO01.test.fastq.gz']}}},
                                        'references': {'grch38': {'anno_bed': None,
                                                                  'anno_bed12': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.gtf.bed12',
                                                                  'anno_gtf': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.gtf',
                                                                  'anno_refflat': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.refFlat.txt',
                                                                  'anno_ribosome_rna_bed': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/hg38_rRNA.bed',
                                                                  'anno_ribosome_rna_interval': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/hg38_rRNA.bed.locations.txt',
                                                                  'chrom_size': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/GRCh38.primary_assembly.genome.chr11.chrom.sizes',
                                                                  'chrs': 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM',
                                                                  'genome': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genome.fa',
                                                                  'genomedict': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genome.dict',
                                                                  'star_index': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/STAR'}},
                                        'star': {'KO01': {'SRR1205282': {'log_file_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/Aln/star/KO01/SRR1205282/run.alignRead.star.SRR1205282.log',
                                                                         'out_dir_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/Aln/star/KO01/SRR1205282',
                                                                         'out_file_path_list': ['/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/Aln/star/KO01/SRR1205282/SRR1205282.star.bam'],
                                                                         'shell_script_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/Shell/run.alignRead.star.SRR1205282.sh'}},
                                                 'main': {'out_step_dir': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/Aln/star'}}}}
                        '''

                        local_resource = ''
                        if args.system_type == 'cluster':
                            pbs_ppn = min([project_cfg['pipeline']['software'][aligner]['threads'], args.threads])
                            local_resource = ', cpus := {}, mem := {}*G, timeout := {}*hour, taskName := "{}.{}"'.format(
                                pbs_ppn,
                                pbs_ppn * 16,
                                72,
                                aligner,
                                read_group
                            )

                        out_file_path = args.submitter
                        try:
                            with open(out_file_path, "a") as outfile:
                                outfile.write('''
// Read Alignment: [ {} ] using [ {} ]

dep( [ '{}' ] <- [ '{}' ]{}) sys bash {}; sleep {}
goal( [ '{}' ] )
\n'''.format(
                 read_group,
                 aligner,
                 "', '".join(task_cfg['align_read'][aligner][library][read_group]['out_file_path_list']),
                 "', '".join(task_cfg['align_read']['main'][library][read_group]['in_file_path_list']),
                 local_resource,
                 task_cfg['align_read'][aligner][library][read_group]['shell_script_path'],
                 project_cfg['pipeline']['software']['bigdatascript']['safeSleep'],
                 "', '".join(task_cfg['align_read'][aligner][library][read_group]['out_file_path_list'])
             ))
                                outfile.close()
                        except IOError as exc:
                            print(exc)


        if 'max_multiple_hit' not in project_cfg['pipeline']['software'][aligner]:
            project_cfg['pipeline']['software'][aligner]['max_multiple_hit'] =  \
                project_cfg['pipeline']['software']['main']['max_multiple_hit']

        logging.debug("[ {} ] Read Alignment:: Aligner: {} - DONE\n".format(SELF_FILE, aligner))

    out_file_path = args.submitter
    try:
        with open(out_file_path, "a") as outfile:
            outfile.write("\n// <<< Read Alignment\n\n\n")
            outfile.close()
    except IOError as exc:
        print(exc)

    #util.ddictfunc.pprint_ddicts(task_cfg, ['align_read'])

    logging.debug("[ {} ] Read Alignment - DONE\n".format(SELF_FILE))



    logging.info("[ {} ] Alignment Merging\n".format(SELF_FILE))

    out_file_path = args.submitter
    try:
        with open(out_file_path, "a") as outfile:
            outfile.write("\n// Alignment Merging >>>\n\n")
            outfile.close()
    except IOError as exc:
        print(exc)

    for aligner in project_cfg['aligners']:
        logging.info("[ {} ] Alignment Merging:: Aligner: {}\n".format(SELF_FILE, aligner))

        task_cfg['merge_aln'][aligner]['main']['out_step_dir'] = \
            '{}/Aln/{}'.format(
                appl_dir_path,
                aligner
            )
        out_dir_t_path = task_cfg['merge_aln'][aligner]['main']['out_step_dir']
        if not os.path.exists(out_dir_t_path):
            logging.debug("MKDIR: [ {} ]".format(out_dir_t_path))
            os.makedirs(out_dir_t_path, exist_ok=True)

        for group in ddicts_metadata.keys():
            for sample in ddicts_metadata[group].keys():
                for library in ddicts_metadata[group][sample].keys():
                    task_cfg['merge_aln'][aligner][library]['in_file_path_list'] = []

                    for read_group in ddicts_metadata[group][sample][library].keys():
                        task_cfg['merge_aln'][aligner][library]['in_file_path_list'].extend(
                            task_cfg['align_read'][aligner][library][read_group]['out_file_path_list'].copy()
                        )

                    task_cfg['merge_aln'][aligner][library]['out_dir_path'] = \
                        '{}/{}'.format(
                            task_cfg['merge_aln'][aligner]['main']['out_step_dir'],
                            library
                        )
                    out_dir_t_path = task_cfg['merge_aln'][aligner][library]['out_dir_path']
                    if not os.path.exists(out_dir_t_path):
                        logging.debug("MKDIR: [ {} ]".format(out_dir_t_path))
                        os.makedirs(out_dir_t_path, exist_ok=True)

                    task_cfg['merge_aln'][aligner][library]['out_file_path_list'] = \
                        [
                            '{}/{}.{}.bam'.format(
                                task_cfg['merge_aln'][aligner][library]['out_dir_path'],
                                library,
                                aligner
                            )
                        ]

                    task_cfg['merge_aln'][aligner][library]['shell_script_path'] = \
                        '{}/run.mergeAln.{}.{}.sh'.format(shell_dir_path, aligner, library)

                    task_cfg['merge_aln'][aligner][library]['log_file_path'] = \
                        '{}/run.mergeAln.{}.{}.log'.format(
                            task_cfg['merge_aln'][aligner][library]['out_dir_path'],
                            aligner,
                            library
                        )

                    #util.ddictfunc.pprint_ddicts(task_cfg, ['merge_aln'])
                    '''
                    {'merge_aln': {'star': {'KO01': {'in_file_path_list': ['/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/Aln/star/KO01/SRR1205282/SRR1205282.star.bam'],
                                                     'log_file_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/Aln/star/KO01/run.mergeAln.star.KO01.log',
                                                     'out_dir_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/Aln/star/KO01',
                                                     'out_file_path_list': ['/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/Aln/star/KO01/KO01.star.bam'],
                                                     'shell_script_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/Shell/run.mergeAln.star.KO01.sh'},
                                            'main': {'out_step_dir': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/Aln/star'}}}}
                    '''

                    local_resource = ''
                    if args.system_type == 'cluster':
                        pbs_ppn = min([project_cfg['pipeline']['software'][aligner]['threads'], args.threads])
                        local_resource = ', cpus := {}, mem := {}*G, timeout := {}*hour, taskName := "{}.{}"'.format(
                            pbs_ppn,
                            pbs_ppn * 8,
                            72,
                            aligner,
                            library
                        )

                    out_file_path = args.submitter
                    try:
                        with open(out_file_path, "a") as outfile:
                            outfile.write('''\n
// Alignment Merging: [ {} ] using [ {} ]

dep( [ '{}' ] <- [ '{}' ]{}) sys bash {}; sleep {}
goal( [ '{}' ] )
\n'''.format(
                 library,
                 aligner,
                 "', '".join([re.sub(r'\.bam$', '.bai', s) for s in task_cfg['merge_aln'][aligner][library]['out_file_path_list']]),
                 "', '".join(task_cfg['merge_aln'][aligner][library]['in_file_path_list']),
                 local_resource,
                 task_cfg['merge_aln'][aligner][library]['shell_script_path'],
                 project_cfg['pipeline']['software']['bigdatascript']['safeSleep'],
                 "', '".join([re.sub(r'\.bam$', '.bai', s) for s in task_cfg['merge_aln'][aligner][library]['out_file_path_list']])
             )
                            )
                            outfile.close()
                    except IOError as exc:
                        print(exc)

        logging.debug("[ {} ] Alignment Merging:: Aligner: {} - DONE\n".format(SELF_FILE, aligner))

    out_file_path = args.submitter
    try:
        with open(out_file_path, "a") as outfile:
            outfile.write("\n// <<< Alignment Merging\n\n\n")
            outfile.close()
    except IOError as exc:
        print(exc)

    #util.ddictfunc.pprint_ddicts(task_cfg, ['merge_aln'])

    logging.debug("[ {} ] Alignment Merging - DONE\n".format(SELF_FILE))



    logging.info("[ {} ] Alignment QC\n".format(SELF_FILE))

    out_file_path = args.submitter
    try:
        with open(out_file_path, "a") as outfile:
            outfile.write("\n// Alignment QC>>>\n\n\n")
            outfile.close()
    except IOError as exc:
        print(exc)

    if("RNAseq" == proj_appl):
        dict_aln_qc_tool = \
            {
                "picard":
                [
                    "CollectRnaSeqMetrics"
                ],
                "rseqc":
                [
                    "clipping_profile.py",
                    #"geneBody_coverage.py", # DISABLE due to long running time
                    #"RPKM_saturation.py", # DISABLE due to long running time
                    "infer_experiment.py"
                ]
            }
    else:
        logging.error(
            "[ {} ] cannot recognize [{}]\n".format(
                SELF_FILE,
                proj_appl
            )
        )
        sys.exit()

    task_cfg['aln_qc']['references'][target_genome] = project_cfg['pipeline']['references'][target_genome]

    for aligner in project_cfg['aligners']:
        for group in ddicts_metadata.keys():
            for sample in ddicts_metadata[group].keys():
                for library in ddicts_metadata[group][sample].keys():
                    task_cfg['aln_qc']['main'][aligner][library]['in_file_path_list'] = \
                        task_cfg['merge_aln'][aligner][library]['out_file_path_list'].copy()

    for alnqctool,module_list in dict_aln_qc_tool.items():
        task_cfg['aln_qc'][alnqctool]['main']['module'] = module_list

        logging.info(
            "[ {} ] Alignment QC:: Alignment QC Tool: {} with modules: {}\n".format(
                SELF_FILE,
                alnqctool,
                ', '.join(module_list)
            )
        )

        for aligner in project_cfg['aligners']:
            logging.info("[ {} ] Alignment QC:: Aligner: {}\n".format(SELF_FILE, aligner))

            task_cfg['aln_qc'][alnqctool][aligner]['main']['out_step_dir'] = \
                '{}/AlnQC/{}/{}'.format(
                    appl_dir_path,
                    alnqctool,
                    aligner
                )
            out_dir_t_path = task_cfg['aln_qc'][alnqctool][aligner]['main']['out_step_dir']
            if not os.path.exists(out_dir_t_path):
                logging.debug("MKDIR: [ {} ]".format(out_dir_t_path))
                os.makedirs(out_dir_t_path, exist_ok=True)

            for group in ddicts_metadata.keys():
                for sample in ddicts_metadata[group].keys():
                    for library in ddicts_metadata[group][sample].keys():
                        logging.info("[ {} ] Alignment QC:: Library: {}\n".format(SELF_FILE, library))

                        task_cfg['aln_qc'][alnqctool][aligner][library]['main']['out_dir_path'] = \
                            '{}/{}'.format(
                                task_cfg['aln_qc'][alnqctool][aligner]['main']['out_step_dir'],
                                library
                            )
                        out_dir_t_path = task_cfg['aln_qc'][alnqctool][aligner][library]['main']['out_dir_path']
                        if not os.path.exists(out_dir_t_path):
                            logging.debug("MKDIR: [ {} ]".format(out_dir_t_path))
                            os.makedirs(out_dir_t_path, exist_ok=True)

                        task_cfg['aln_qc'][alnqctool][aligner][library]['main']['tmp_dir_path'] = \
                            '{}/{}'.format(
                                task_cfg['aln_qc'][alnqctool][aligner][library]['main']['out_dir_path'],
                                "tmp"
                            )
                        out_dir_t_path = task_cfg['aln_qc'][alnqctool][aligner][library]['main']['tmp_dir_path']
                        if not os.path.exists(out_dir_t_path):
                            logging.debug("MKDIR: [ {} ]".format(out_dir_t_path))
                            os.makedirs(out_dir_t_path, exist_ok=True)

                        for module in module_list:
                            task_cfg['aln_qc'][alnqctool][aligner][library][module]['out_file_base'] = \
                                '{}/{}.{}.{}'.format(
                                    task_cfg['aln_qc'][alnqctool][aligner][library]['main']['out_dir_path'],
                                    library,
                                    aligner,
                                    alnqctool
                                )

                            if(alnqctool == "picard"):
                                if(module == "CollectRnaSeqMetrics"):
                                    task_cfg['aln_qc'][alnqctool][aligner][library]['main']['strand_specificity'] = \
                                        lib_type_2_strand_specificity_dict[metadata.loc[metadata['Library'] == library]['LibType'].drop_duplicates().values.tolist()[0]][alnqctool]

                                    suffix_list = [
                                        "RNA_Metrics",
                                        "RNA_Metrics.pdf"
                                    ]
                            elif(alnqctool == "rseqc"):
                                if(module == "clipping_profile.py"):
                                    task_cfg['aln_qc'][alnqctool][aligner][library]['main']['sequencing_layout'] = \
                                        flavor_2_sequencing_layout_dict[str(metadata.loc[metadata['Library'] == library]['Flavor'].drop_duplicates().values.tolist()[0][0])]['''{}__{}'''.format(alnqctool, module)]

                                    suffix_list = [
                                        "clipping_profile.xls",
                                        "clipping_profile.r"
                                    ]
                                    
                                    if(str(metadata.loc[metadata['Library'] == library]['Flavor'].drop_duplicates().values.tolist()[0][0]) == "1"):
                                        suffix_list = \
                                        suffix_list + ["clipping_profile.pdf"]
                                    else:
                                        suffix_list = \
                                        suffix_list + [
                                            "clipping_profile.R1.pdf",
                                            "clipping_profile.R2.pdf"
                                        ]
                                elif(module == "geneBody_coverage.py"):
                                    suffix_list = [
                                        "geneBodyCoverage.txt",
                                        "geneBodyCoverage.r",
                                        "geneBodyCoverage.curves.pdf"
                                    ]
                                elif(module == "infer_experiment.py"):
                                    suffix_list = [
                                        "infer_experiment.txt"
                                    ]
                                elif(module == "RPKM_saturation.py"):
                                    task_cfg['aln_qc'][alnqctool][aligner][library]['main']['strand_specificity'] = \
                                        lib_type_2_strand_specificity_dict[metadata.loc[metadata['Library'] == library]['LibType'].drop_duplicates().values.tolist()[0]]['''{}__{}'''.format(alnqctool, module)]

                                    suffix_list = [
                                        "eRPKM.xls",
                                        "rawCount.xls",
                                        "saturation.r"#,
                                        #"saturation.pdf" # due to the hard coding in RSeQC::RPKM_saturation.py, skipping the checking of this file
                                    ]

                            task_cfg['aln_qc'][alnqctool][aligner][library][module]['out_file_path_list'] = \
                                [
                                    "{}.{}".format(p, s) for p, s in list(zip(itertools.cycle([task_cfg['aln_qc'][alnqctool][aligner][library][module]['out_file_base']]), suffix_list))
                                ]

                            task_cfg['aln_qc'][alnqctool][aligner][library][module]['shell_script_path'] = \
                                '{}/run.alnQC.{}.{}.{}.{}.sh'.format(
                                    shell_dir_path,
                                    alnqctool,
                                    aligner,
                                    library,
                                    module
                                )

                            task_cfg['aln_qc'][alnqctool][aligner][library][module]['log_file_path'] = \
                                '{}/run.alnQC.{}.{}.{}.{}.log'.format(
                                    task_cfg['aln_qc'][alnqctool][aligner][library]['main']['out_dir_path'],
                                    alnqctool,
                                    aligner,
                                    library,
                                    module
                                )

                        #util.ddictfunc.pprint_ddicts(task_cfg, ['aln_qc'])
                        '''
                        {'aln_qc': {'main': {'star': {'KO01': {'in_file_path_list': ['/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/Aln/star/KO01/KO01.star.bam']},
                                                      'WT01': {'in_file_path_list': ['/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/Aln/star/WT01/WT01.star.bam']}}},
                                    'picard': {'main': {'module': ['CollectRnaSeqMetrics']},
                                               'star': {'KO01': {'CollectRnaSeqMetrics': {'log_file_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/AlnQC/picard/star/KO01/run.alnQC.picard.star.KO01.CollectRnaSeqMetrics.log',
                                                                                          'out_file_base': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/AlnQC/picard/star/KO01/KO01.star.picard',
                                                                                          'out_file_path_list': ['/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/AlnQC/picard/star/KO01/KO01.star.picard.RNA_Metrics',
                                                                                                                 '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/AlnQC/picard/star/KO01/KO01.star.picard.RNA_Metrics.pdf'],
                                                                                          'shell_script_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/Shell/run.alnQC.picard.star.KO01.CollectRnaSeqMetrics.sh'},
                                                                 'main': {'out_dir_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/AlnQC/picard/star/KO01',
                                                                          'strand_specificity': 'NONE',
                                                                          'tmp_dir_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/AlnQC/picard/star/KO01/tmp'}},
                                                        'main': {'out_step_dir': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/AlnQC/picard/star'}}},
                                    'references': {'grch38': {'anno_bed': None,
                                                              'anno_bed12': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.gtf.bed12',
                                                              'anno_gtf': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.gtf',
                                                              'anno_refflat': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.refFlat.txt',
                                                              'anno_ribosome_rna_bed': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/hg38_rRNA.bed',
                                                              'anno_ribosome_rna_interval': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/hg38_rRNA.bed.locations.txt',
                                                              'chrom_size': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/GRCh38.primary_assembly.genome.chr11.chrom.sizes',
                                                              'chrs': 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM',
                                                              'genome': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genome.fa',
                                                              'genomedict': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genome.dict',
                                                              'star_index': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/STAR'}}}}
                        '''
                        
                        local_resource = ''
                        if args.system_type == 'cluster':
                            pbs_ppn = min([project_cfg['pipeline']['software'][alnqctool]['threads'], args.threads])
                            local_resource = ', cpus := {}, mem := {}*G, timeout := {}*hour, taskName := "{}.{}.{}"'.format(
                                pbs_ppn,
                                pbs_ppn * 8,
                                72,
                                alnqctool,
                                aligner,
                                library
                            )

                        for module in module_list:
                            out_file_path = args.submitter
                            try:
                                with open(out_file_path, "a") as outfile:
                                    outfile.write('''
// Alignment QC: [ {} ] aligned by [ {} ] using [ {}::{} ]

dep( [ '{}' ] <- [ '{}' ]{}) sys bash {}; sleep {}
goal( [ '{}' ] )
\n'''.format(
                 library,
                 aligner,
                 alnqctool,
                 module,
                 "', '".join(task_cfg['aln_qc'][alnqctool][aligner][library][module]['out_file_path_list']),
                 "', '".join([re.sub(r'\.bam$', '.bai', s) for s in task_cfg['aln_qc']['main'][aligner][library]['in_file_path_list']]),
                 local_resource,
                 task_cfg['aln_qc'][alnqctool][aligner][library][module]['shell_script_path'],
                 project_cfg['pipeline']['software']['bigdatascript']['safeSleep'],
                 "', '".join(task_cfg['aln_qc'][alnqctool][aligner][library][module]['out_file_path_list'])
             ))
                                    outfile.close()
                            except IOError as exc:
                                print(exc)

                        logging.debug("[ {} ] Alignment QC:: Library: {} - DONE\n".format(SELF_FILE, library))
            logging.debug("[ {} ] Alignment QC:: Aligner: {} - DONE\n".format(SELF_FILE, aligner))

        logging.debug(
            "[ {} ] Alignment QC:: Alignment QC Tool: {} with modules: {} - DONE\n".format(
                SELF_FILE,
                alnqctool,
                ', '.join(module_list)
            )
        )

    out_file_path = args.submitter
    try:
        with open(out_file_path, "a") as outfile:
            outfile.write("\n// <<< Alignment QC\n\n\n")
            outfile.close()
    except IOError as exc:
        print(exc)

    #util.ddictfunc.pprint_ddicts(task_cfg, ['aln_qc'])

    logging.debug("[ {} ] Alignment QC - DONE\n".format(SELF_FILE))



    logging.info("[ {} ] Quantification\n".format(SELF_FILE))

    out_file_path = args.submitter
    try:
        with open(out_file_path, "a") as outfile:
            outfile.write("\n// Quantification >>>\n\n\n")
            outfile.close()
    except IOError as exc:
        print(exc)

    task_cfg['quant']['references'][target_genome] = project_cfg['pipeline']['references'][target_genome]

    for aligner in project_cfg['aligners']:
        for group in ddicts_metadata.keys():
            for sample in ddicts_metadata[group].keys():
                for library in ddicts_metadata[group][sample].keys():
                    task_cfg['quant']['main'][aligner][library]['in_file_path_list'] = \
                        task_cfg['merge_aln'][aligner][library]['out_file_path_list'].copy()

    for quantifier in project_cfg['quantifiers']:
        logging.info(
            "[ {} ] Quantification:: Quantifier: {}\n".format(
                SELF_FILE,
                quantifier
            )
        )

        for aligner in project_cfg['aligners']:
            logging.info("[ {} ] Quantification:: Aligner: {}\n".format(SELF_FILE, aligner))

            task_cfg['quant'][quantifier][aligner]['main']['out_step_dir'] = \
                '{}/Quantification/{}/{}'.format(
                    appl_dir_path,
                    quantifier,
                    aligner
                )
            out_dir_t_path = task_cfg['quant'][quantifier][aligner]['main']['out_step_dir']
            if not os.path.exists(out_dir_t_path):
                logging.debug("MKDIR: [ {} ]".format(out_dir_t_path))
                os.makedirs(out_dir_t_path, exist_ok=True)

            for group in ddicts_metadata.keys():
                for sample in ddicts_metadata[group].keys():
                    for library in ddicts_metadata[group][sample].keys():
                        logging.info("[ {} ] Quantification:: Library: {}\n".format(SELF_FILE, library))

                        task_cfg['quant'][quantifier][aligner][library]['out_dir_path'] = \
                            '{}/{}'.format(
                                task_cfg['quant'][quantifier][aligner]['main']['out_step_dir'],
                                library
                            )
                        out_dir_t_path = task_cfg['quant'][quantifier][aligner][library]['out_dir_path']
                        if not os.path.exists(out_dir_t_path):
                            logging.debug("MKDIR: [ {} ]".format(out_dir_t_path))
                            os.makedirs(out_dir_t_path, exist_ok=True)

                        #task_cfg['quant'][quantifier][aligner][library]['tmp_dir_path'] = \
                            #'{}/{}'.format(
                                #task_cfg['quant'][quantifier][aligner][library]['out_dir_path'],
                                #"tmp"
                            #)
                        #out_dir_t_path = task_cfg['quant'][quantifier][aligner][library]['tmp_dir_path']
                        #if not os.path.exists(out_dir_t_path):
                            #logging.debug("MKDIR: [ {} ]".format(out_dir_t_path))
                            #os.makedirs(out_dir_t_path, exist_ok=True)

                        if quantifier == "featurecounts":
                            task_cfg['quant'][quantifier][aligner][library]['strandSpecific'] = \
                                lib_type_2_strand_specificity_dict[metadata.loc[metadata['Library'] == library]['LibType'].drop_duplicates().values.tolist()[0]][quantifier]

                            task_cfg['quant'][quantifier][aligner][library]['isPairedEnd'] = \
                                flavor_2_sequencing_layout_dict[str(metadata.loc[metadata['Library'] == library]['Flavor'].drop_duplicates().values.tolist()[0][0])][quantifier]

                        task_cfg['quant'][quantifier][aligner][library]['out_file_base'] = \
                            '{}/{}.{}.{}'.format(
                                task_cfg['quant'][quantifier][aligner][library]['out_dir_path'],
                                library,
                                aligner,
                                quantifier
                            )

                        if quantifier == "featurecounts":
                            ext_str = "count"

                        task_cfg['quant'][quantifier][aligner][library]['out_file_path_list'] = \
                            [
                                "{}.{}".format(
                                    task_cfg['quant'][quantifier][aligner][library]['out_file_base'],
                                    ext_str
                                )
                            ]

                        task_cfg['quant'][quantifier][aligner][library]['shell_script_path'] = \
                            '{}/run.quant.{}.{}.{}.sh'.format(
                                shell_dir_path,
                                quantifier,
                                aligner,
                                library
                            )

                        task_cfg['quant'][quantifier][aligner][library]['log_file_path'] = \
                            '{}/run.quant.{}.{}.{}.log'.format(
                                task_cfg['quant'][quantifier][aligner][library]['out_dir_path'],
                                quantifier,
                                aligner,
                                library
                            )

                        #util.ddictfunc.pprint_ddicts(task_cfg, ['quant'])
                        '''
                        {'quant': {'featurecounts': {'star': {'KO01': {'isPairedEnd': False,
                                                                       'log_file_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/Quantification/featurecounts/star/KO01/run.quant.featurecounts.star.KO01.log',
                                                                       'out_dir_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/Quantification/featurecounts/star/KO01',
                                                                       'out_file_base': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/Quantification/featurecounts/star/KO01/KO01.star.featurecounts',
                                                                       'out_file_path_list': ['/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/Quantification/featurecounts/star/KO01/KO01.star.featurecounts.count'],
                                                                       'shell_script_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/Shell/run.quant.featurecounts.star.KO01.sh',
                                                                       'strandSpecific': 0},
                                                              'main': {'out_step_dir': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/Quantification/featurecounts/star'}}},
                                   'main': {'star': {'KO01': {'in_file_path_list': ['/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/Aln/star/KO01/KO01.star.bam']},
                                                     'WT01': {'in_file_path_list': ['/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/Aln/star/WT01/WT01.star.bam']}}},
                                   'references': {'grch38': {'anno_bed': None,
                                                             'anno_bed12': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.bed12',
                                                             'anno_gtf': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.gtf',
                                                             'anno_refflat': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.refFlat.txt',
                                                             'anno_ribosome_rna_bed': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.rRNA.bed',
                                                             'anno_ribosome_rna_interval': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.rRNA.interval_list',
                                                             'chrom_size': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/GRCh38.primary_assembly.genome.chr11.chrom.size',
                                                             'chrs': 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM',
                                                             'genome': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/GRCh38.primary_assembly.genome.chr11.fa',
                                                             'genomedict': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/GRCh38.primary_assembly.genome.chr11.dict',
                                                             'star_index': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12'}}}}

                        '''

                        local_resource = ''
                        if args.system_type == 'cluster':
                            pbs_ppn = min([project_cfg['pipeline']['software'][quantifier]['threads'], args.threads])
                            local_resource = ', cpus := {}, mem := {}*G, timeout := {}*hour, taskName := "{}.{}.{}"'.format(
                                pbs_ppn,
                                pbs_ppn * 8,
                                72,
                                quantifier,
                                aligner,
                                library
                            )

                        out_file_path = args.submitter
                        try:
                            with open(out_file_path, "a") as outfile:
                                outfile.write('''
// Quantification: [ {} ] aligned by [ {} ] using [ {} ]

dep( [ '{}' ] <- [ '{}' ]{}) sys bash {}; sleep {}
goal( [ '{}' ] )
\n'''.format(
                 library,
                 aligner,
                 quantifier,
                 "', '".join(task_cfg['quant'][quantifier][aligner][library]['out_file_path_list']),
                 "', '".join([re.sub(r'\.bam$', '.bai', s) for s in task_cfg['quant']['main'][aligner][library]['in_file_path_list']]),
                 local_resource,
                 task_cfg['quant'][quantifier][aligner][library]['shell_script_path'],
                 project_cfg['pipeline']['software']['bigdatascript']['safeSleep'],
                 "', '".join(task_cfg['quant'][quantifier][aligner][library]['out_file_path_list'])
             ))
                                outfile.close()
                        except IOError as exc:
                            print(exc)

                        logging.debug("[ {} ] Quantification:: Library: {} - DONE\n".format(SELF_FILE, library))
            logging.debug("[ {} ] Quantification:: Aligner: {} - DONE\n".format(SELF_FILE, aligner))

        logging.debug(
            "[ {} ] Quantification:: Quantifier: {} - DONE\n".format(
                SELF_FILE,
                quantifier
            )
        )

    out_file_path = args.submitter
    try:
        with open(out_file_path, "a") as outfile:
            outfile.write("\n// <<< Quantification\n\n\n")
            outfile.close()
    except IOError as exc:
        print(exc)

    #util.ddictfunc.pprint_ddicts(task_cfg, ['quant'])

    logging.debug("[ {} ] Quantification - DONE\n".format(SELF_FILE))



    logging.info("[ {} ] Calling\n".format(SELF_FILE))

    out_file_path = args.submitter
    try:
        with open(out_file_path, "a") as outfile:
            outfile.write("\n// Calling >>>\n\n\n")
            outfile.close()
    except IOError as exc:
        print(exc)

    task_cfg['call']['references'][target_genome] = project_cfg['pipeline']['references'][target_genome]
    task_cfg['call']['main']['metadata_path'] = args.metadata_path
    task_cfg['call']['main'].update({'application': proj_appl})
    task_cfg['call']['main']['comparison'] = \
        project_cfg['project']['comparison']
    task_cfg['call']['main'].update(project_cfg['pipeline']['parameters'][proj_appl])

    for quantifier in project_cfg['quantifiers']:
        for aligner in project_cfg['aligners']:
            task_cfg['call']['main'][quantifier][aligner]['in_dir_path'] = \
                task_cfg['quant'][quantifier][aligner]['main']['out_step_dir']

            task_cfg['call']['main'][quantifier][aligner]['in_file_path_list'] = []

            for group in ddicts_metadata.keys():
                for sample in ddicts_metadata[group].keys():
                    for library in ddicts_metadata[group][sample].keys():
                        logging.info("[ {} ] Quantification:: Library: {}\n".format(SELF_FILE, library))

                        task_cfg['call']['main'][quantifier][aligner]['in_file_path_list'].extend(
                            task_cfg['quant'][quantifier][aligner][library]['out_file_path_list'].copy()
                        )

    for caller in project_cfg['callers']:
        logging.info(
            "[ {} ] Calling:: Caller: {}\n".format(
                SELF_FILE,
                caller
            )
        )

        for quantifier in project_cfg['quantifiers']:
            logging.info(
                "[ {} ] Calling:: Quantifier: {}\n".format(
                    SELF_FILE,
                    quantifier
                )
            )

            for aligner in project_cfg['aligners']:
                logging.info("[ {} ] Calling:: Aligner: {}\n".format(SELF_FILE, aligner))

                task_cfg['call'][caller][quantifier][aligner]['out_dir_path'] = \
                    '{}/DEG/{}/{}/{}'.format(
                        appl_dir_path,
                        caller,
                        quantifier,
                        aligner
                    )
                out_dir_t_path = task_cfg['call'][caller][quantifier][aligner]['out_dir_path']
                if not os.path.exists(out_dir_t_path):
                    logging.debug("MKDIR: [ {} ]".format(out_dir_t_path))
                    os.makedirs(out_dir_t_path, exist_ok=True)

                task_cfg['call'][caller][quantifier][aligner]['out_file_path_list'] = \
                    [
                        '{}/{}.{}.{}.{}.{}'.format(
                            task_cfg['call'][caller][quantifier][aligner]['out_dir_path'],
                            project_cfg['project']['name'],
                            aligner,
                            quantifier,
                            caller,
                            "count.txt"
                        )
                    ]

                task_cfg['call'][caller][quantifier][aligner]['shell_script_path'] = \
                    '{}/run.call.{}.{}.{}.{}.sh'.format(
                        shell_dir_path,
                        caller,
                        quantifier,
                        aligner,
                        project_cfg['project']['name']
                    )

                task_cfg['call'][caller][quantifier][aligner]['log_file_path'] = \
                    '{}/run.call.{}.{}.{}.{}.log'.format(
                        task_cfg['call'][caller][quantifier][aligner]['out_dir_path'],
                        caller,
                        quantifier,
                        aligner,
                        project_cfg['project']['name']
                    )

                #for comparison in project_cfg['project']['comparison']:
                    #case_group, ctrl_group = comparison.split('_vs_', 1)
                    #out_dir_t_path = \
                        #'''{}/{}'''.format(
                            #task_cfg['call'][caller][quantifier][aligner]['out_dir_path'],
                            #comparison
                        #)
                    #if not os.path.exists(out_dir_t_path):
                        #logging.debug("MKDIR: [ {} ]".format(out_dir_t_path))
                        #os.makedirs(out_dir_t_path, exist_ok=True)

                #util.ddictfunc.pprint_ddicts(task_cfg, ['call'])
                '''
                {'call': {'edger': {'featurecounts': {'star': {'log_file_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/DEG/edger/featurecounts/star/run.call.edger.featurecounts.star.DLBC.log',
                                                               'out_dir_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/DEG/edger/featurecounts/star',
                                                               'out_file_path_list': ['/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/DEG/edger/featurecounts/star/DLBC.star.featurecounts.edger.count.txt'],
                                                               'shell_script_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/shell_scripts/run.call.edger.featurecounts.star.DLBC.sh'}}},
                          'main': {'application': 'RNAseq',
                                   'comparison': {'KO_vs_WT': {'criteria': None,
                                                               'filter_4_low_expr': None,
                                                               'fold_change': None,
                                                               'pval': None}},
                                   'criteria': 'fc_q',
                                   'featurecounts': {'star': {'in_dir_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/Quantification/featurecounts/star',
                                                              'in_file_path_list': ['/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/Quantification/featurecounts/star/KO01/KO01.star.featurecounts.count',
                                                                                    '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/Quantification/featurecounts/star/KO02/KO02.star.featurecounts.count',
                                                                                    '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/Quantification/featurecounts/star/KO03/KO03.star.featurecounts.count',
                                                                                    '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/Quantification/featurecounts/star/WT01/WT01.star.featurecounts.count',
                                                                                    '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/Quantification/featurecounts/star/WT02/WT02.star.featurecounts.count',
                                                                                    '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/Quantification/featurecounts/star/WT03/WT03.star.featurecounts.count']}},
                                   'filter_4_low_expr': True,
                                   'fold_change': 1.5,
                                   'metadata_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC.metadata.txt',
                                   'pval': 0.1},
                          'references': {'grch38': {'anno_bed': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.bed12',
                                                    'anno_bed12': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.bed12',
                                                    'anno_gtf': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.gtf',
                                                    'anno_refflat': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.refFlat.txt',
                                                    'anno_ribosome_rna_bed': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.rRNA.bed',
                                                    'anno_ribosome_rna_interval': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.rRNA.interval_list',
                                                    'chrom_size': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/GRCh38.primary_assembly.genome.chr11.chrom.size',
                                                    'chrs': 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM',
                                                    'genome': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/GRCh38.primary_assembly.genome.chr11.fa',
                                                    'genomedict': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/GRCh38.primary_assembly.genome.chr11.dict',
                                                    'star_index': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12'}}}}
                '''

                local_resource = ''
                if args.system_type == 'cluster':
                    pbs_ppn = min([project_cfg['pipeline']['software'][caller]['threads'], args.threads])
                    local_resource = ', cpus := {}, mem := {}*G, timeout := {}*hour, taskName := "{}.{}.{}.{}"'.format(
                        pbs_ppn,
                        pbs_ppn * 8,
                        72,
                        caller,
                        quantifier,
                        aligner,
                        project_cfg['project']['name']
                    )

                out_file_path = args.submitter
                try:
                    with open(out_file_path, "a") as outfile:
                        outfile.write('''
// Calling: [ {} ] aligned by [ {} ], quantified by  [ {} ], and called by [ {} ]

dep( [ '{}' ] <- [ '{}' ]{}) sys bash {}; sleep {}
goal( [ '{}' ] )
\n'''.format(
         project_cfg['project']['name'],
         aligner,
         quantifier,
         caller,
         "', '".join(task_cfg['call'][caller][quantifier][aligner]['out_file_path_list'] + [re.sub(r'.count.txt$', '.test.DEG.txt', file) for file in task_cfg['call'][caller][quantifier][aligner]['out_file_path_list']]),
         "', '".join(task_cfg['call']['main'][quantifier][aligner]['in_file_path_list']),
         local_resource,
         task_cfg['call'][caller][quantifier][aligner]['shell_script_path'],
         project_cfg['pipeline']['software']['bigdatascript']['safeSleep'],
         "', '".join(task_cfg['call'][caller][quantifier][aligner]['out_file_path_list'])
     ))
                        outfile.close()
                except IOError as exc:
                    print(exc)

                logging.debug(
                    "[ {} ] Calling:: Aligner: {} - DONE\n".format(
                        SELF_FILE,
                        aligner
                    )
                )

            logging.debug(
                "[ {} ] Calling:: Quantifier: {} - DONE\n".format(
                    SELF_FILE,
                    quantifier
                )
            )

        logging.debug(
            "[ {} ] Calling:: Caller: {} - DONE\n".format(
                SELF_FILE,
                caller
            )
        )

    out_file_path = args.submitter
    try:
        with open(out_file_path, "a") as outfile:
            outfile.write("\n// <<< Calling\n\n\n")
            outfile.close()
    except IOError as exc:
        print(exc)

    #util.ddictfunc.pprint_ddicts(task_cfg, ['call'])

    logging.debug("[ {} ] Calling - DONE\n".format(SELF_FILE))



    logging.info("[ {} ] Quantification QC\n".format(SELF_FILE))

    out_file_path = args.submitter
    try:
        with open(out_file_path, "a") as outfile:
            outfile.write("\n// Quantification QC >>>\n\n\n")
            outfile.close()
    except IOError as exc:
        print(exc)

    if("RNAseq" == proj_appl):
        quant_qc_tool = \
            [
                "pca"
            ]

        caller = 'deseq2' if 'deseq2' in project_cfg['callers'] else project_cfg['callers'][1]
    else:
        logging.error(
            "[ {} ] cannot recognize [{}]\n".format(
                SELF_FILE,
                proj_appl
            )
        )
        sys.exit()

    task_cfg['quant_qc']['references'][target_genome] = project_cfg['pipeline']['references'][target_genome]
    task_cfg['quant_qc']['main']['metadata_path'] = args.metadata_path
    task_cfg['quant_qc']['main'].update(
        {
            'application': proj_appl,
            'library_info': dict_library_info
        }
    )

    for quantifier in project_cfg['quantifiers']:
        for aligner in project_cfg['aligners']:
            task_cfg['quant_qc']['main'][quantifier][aligner]['in_file_path_list'] = \
                task_cfg['call'][caller][quantifier][aligner]['out_file_path_list'].copy()

            if project_cfg['project']['run_as_practice']:
                task_cfg['quant_qc']['main'][quantifier][aligner]['in_file_path_list'] = \
                    [
                        re.sub(
                            r'DLBC/RNAseq',
                            'DLBC_full/RNAseq',
                            i
                        ) for i in task_cfg['quant_qc']['main'][quantifier][aligner]['in_file_path_list']
                    ]
                logging.debug("[RUN_AS_PRACTICE] using full set result instead")

            task_cfg['quant_qc']['main'][quantifier][aligner]['out_dir_path'] = \
                '{}/QuantQC/{}/{}'.format(
                    appl_dir_path,
                    quantifier,
                    aligner
                )
            out_dir_t_path = task_cfg['quant_qc']['main'][quantifier][aligner]['out_dir_path']
            if not os.path.exists(out_dir_t_path):
                logging.debug("MKDIR: [ {} ]".format(out_dir_t_path))
                os.makedirs(out_dir_t_path, exist_ok=True)

    for tool in quant_qc_tool:
        logging.info(
            "[ {} ] Quantification QC:: Quantification QC tool: {}\n".format(
                SELF_FILE,
                tool
            )
        )
        for quantifier in project_cfg['quantifiers']:
            logging.info(
                "[ {} ] Quantification QC:: Quantifier: {}\n".format(
                    SELF_FILE,
                    quantifier
                )
            )

            for aligner in project_cfg['aligners']:
                logging.info(
                    "[ {} ] Quantification QC:: Aligner: {}\n".format(
                        SELF_FILE,
                        aligner
                    )
                )

                task_cfg['quant_qc'][tool][quantifier][aligner]['out_file_path_list'] = \
                    [
                        '{}/{}.{}.{}.{}.{}'.format(
                            task_cfg['quant_qc']['main'][quantifier][aligner]['out_dir_path'],
                            project_cfg['project']['name'],
                            aligner,
                            quantifier,
                            tool,
                            "pdf"
                        )
                    ]

                task_cfg['quant_qc'][tool][quantifier][aligner]['shell_script_path'] = \
                    '{}/run.quantQC.{}.{}.{}.{}.sh'.format(
                        shell_dir_path,
                        tool,
                        quantifier,
                        aligner,
                        project_cfg['project']['name']
                    )

                task_cfg['quant_qc'][tool][quantifier][aligner]['log_file_path'] = \
                    '{}/run.quantQC.{}.{}.{}.{}.log'.format(
                        task_cfg['quant_qc']['main'][quantifier][aligner]['out_dir_path'],
                        tool,
                        quantifier,
                        aligner,
                        project_cfg['project']['name']
                    )

                #util.ddictfunc.pprint_ddicts(task_cfg, ['quant_qc'])
                '''
                {'quant_qc': {'main': {'application': 'RNAseq',
                                       'featurecounts': {'star': {'in_file_path_list': ['/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/DEG/deseq2/featurecounts/star/DLBC.star.featurecounts.deseq2.count.txt'],
                                                                  'out_dir_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/QuantQC/featurecounts/star'}},
                                       'library_info': {'KO01': 'KO',
                                                        'KO02': 'KO',
                                                        'KO03': 'KO',
                                                        'WT01': 'WT',
                                                        'WT02': 'WT',
                                                        'WT03': 'WT'},
                                       'metadata_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC.metadata.txt'},
                              'pca': {'featurecounts': {'star': {'log_file_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/QuantQC/featurecounts/star/run.quantQC.pca.featurecounts.star.DLBC.log',
                                                                 'out_file_path_list': ['/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/QuantQC/featurecounts/star/DLBC.star.featurecounts.pca.pdf'],
                                                                 'shell_script_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/shell_scripts/run.quantQC.pca.featurecounts.star.DLBC.sh'}}},
                              'references': {'grch38': {'anno_bed': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.bed12',
                                                        'anno_bed12': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.bed12',
                                                        'anno_gtf': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.gtf',
                                                        'anno_refflat': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.refFlat.txt',
                                                        'anno_ribosome_rna_bed': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.rRNA.bed',
                                                        'anno_ribosome_rna_interval': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.rRNA.interval_list',
                                                        'chrom_size': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/GRCh38.primary_assembly.genome.chr11.chrom.size',
                                                        'chrs': 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM',
                                                        'genome': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/GRCh38.primary_assembly.genome.chr11.fa',
                                                        'genomedict': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/GRCh38.primary_assembly.genome.chr11.dict',
                                                        'star_index': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12'}}}}
                '''

                local_resource = ''
                if args.system_type == 'cluster':
                    pbs_ppn = args.threads
                    local_resource = ', cpus := {}, mem := {}*G, timeout := {}*hour, taskName := "{}.{}.{}.{}"'.format(
                        pbs_ppn,
                        pbs_ppn * 8,
                        72,
                        tool,
                        quantifier,
                        aligner,
                        project_cfg['project']['name']
                    )

                out_file_path = args.submitter
                try:
                    with open(out_file_path, "a") as outfile:
                        outfile.write('''
// Quantification QC: [ {} ] aligned by [ {} ], quantified by  [ {} ], and called by [ {} ] using quantificaiton QC tool [ {} ]

dep( [ '{}' ] <- [ '{}' ]{}) sys bash {}; sleep {}
goal( [ '{}' ] )
\n'''.format(
         project_cfg['project']['name'],
         aligner,
         quantifier,
         caller,
         tool,
         "', '".join(task_cfg['quant_qc'][tool][quantifier][aligner]['out_file_path_list']),
         "', '".join(task_cfg['quant_qc']['main'][quantifier][aligner]['in_file_path_list']),
         local_resource,
         task_cfg['quant_qc'][tool][quantifier][aligner]['shell_script_path'],
         project_cfg['pipeline']['software']['bigdatascript']['safeSleep'],
         "', '".join(task_cfg['quant_qc'][tool][quantifier][aligner]['out_file_path_list'])
     ))
                        outfile.close()
                except IOError as exc:
                    print(exc)

                logging.debug(
                    "[ {} ] Quantification QC:: Aligner: {} - DONE\n".format(
                        SELF_FILE,
                        aligner
                    )
                )

            logging.debug(
                "[ {} ] Quantification QC:: Quantifier: {} - DONE\n".format(
                    SELF_FILE,
                    quantifier
                )
            )

        logging.debug(
            "[ {} ] Quantification QC:: Quantification QC tool: {} - DONE\n".format(
                SELF_FILE,
                tool
            )
        )

    out_file_path = args.submitter
    try:
        with open(out_file_path, "a") as outfile:
            outfile.write("\n// <<< Quantification QC\n\n\n")
            outfile.close()
    except IOError as exc:
        print(exc)

    #util.ddictfunc.pprint_ddicts(task_cfg, ['quant_qc'])

    logging.debug("[ {} ] Quantification QC - DONE\n".format(SELF_FILE))


    logging.info("[ {} ] Loci Statistics\n".format(SELF_FILE))

    out_file_path = args.submitter
    try:
        with open(out_file_path, "a") as outfile:
            outfile.write("\n// Loci Statistics >>>\n\n\n")
            outfile.close()
    except IOError as exc:
        print(exc)

    if("RNAseq" == proj_appl):
        loci_stat_tool = \
            [
                "venn"
            ]
        caller = 'deseq2' if 'deseq2' in project_cfg['callers'] else project_cfg['callers'][1]
    else:
        logging.error(
            "[ {} ] cannot recognize [{}]\n".format(
                SELF_FILE,
                proj_appl
            )
        )
        sys.exit()

    task_cfg['loci_stat']['references'][target_genome] = project_cfg['pipeline']['references'][target_genome]
    task_cfg['loci_stat']['main'].update({'application': proj_appl})

    for quantifier in project_cfg['quantifiers']:
        logging.info(
            "[ {} ] Loci Statistics:: Quantifier: {}\n".format(
                SELF_FILE,
                quantifier
            )
        )

        for aligner in project_cfg['aligners']:
            logging.info(
                "[ {} ] Loci Statistics:: Aligner: {}\n".format(
                    SELF_FILE,
                    aligner
                )
            )

            task_cfg['loci_stat'][quantifier][aligner]['main']['anchor_file_path'] = \
                '{}/{}.{}.{}.{}.{}'.format(
                    task_cfg['call'][caller][quantifier][aligner]['out_dir_path'],
                    project_cfg['project']['name'],
                    aligner,
                    quantifier,
                    caller,
                    "test.txt"
                )

            if project_cfg['project']['run_as_practice']:
                task_cfg['loci_stat'][quantifier][aligner]['main']['anchor_file_path'] = \
                    re.sub(
                        r'DLBC/RNAseq',
                        'DLBC_full/RNAseq',
                        task_cfg['loci_stat'][quantifier][aligner]['main']['anchor_file_path']
                    )
                logging.debug("[RUN_AS_PRACTICE] using full set result instead")

            if True:
                proj_comp = project_cfg['project']['name']

                logging.info(
                    "[ {} ] Loci Statistics:: Project: {}\n".format(
                        SELF_FILE,
                        proj_comp
                    )
                )

                task_cfg['loci_stat'][quantifier][aligner][proj_comp]['in_file_path_list'] = []

                for caller in project_cfg['callers']:
                    task_cfg['loci_stat'][quantifier][aligner][proj_comp]['in_file_path_list'].append(
                        '{}/{}.{}.{}.{}.{}'.format(
                            task_cfg['call'][caller][quantifier][aligner]['out_dir_path'],
                            proj_comp,
                            aligner,
                            quantifier,
                            caller,
                            "test.DEG.txt"
                        )
                    )

                if project_cfg['project']['run_as_practice']:
                    task_cfg['loci_stat'][quantifier][aligner][proj_comp]['in_file_path_list'] = \
                        [
                            re.sub(
                                r'DLBC/RNAseq',
                                'DLBC_full/RNAseq',
                                i
                            ) for i in task_cfg['loci_stat'][quantifier][aligner][proj_comp]['in_file_path_list']
                        ]
                    logging.debug("[RUN_AS_PRACTICE] useing full set result instead")

                task_cfg['loci_stat'][quantifier][aligner][proj_comp]['out_dir_path'] = \
                    '{}/LociStat/{}/{}/{}'.format(
                        appl_dir_path,
                        quantifier,
                        aligner,
                        proj_comp
                    )
                out_dir_t_path = task_cfg['loci_stat'][quantifier][aligner][proj_comp]['out_dir_path']
                if not os.path.exists(out_dir_t_path):
                    logging.debug("MKDIR: [ {} ]".format(out_dir_t_path))
                    os.makedirs(out_dir_t_path, exist_ok=True)

                task_cfg['loci_stat'][quantifier][aligner][proj_comp]['out_file_path_list'] = \
                    [
                        '{}/{}.{}.{}.{}'.format(
                            task_cfg['loci_stat'][quantifier][aligner][proj_comp]['out_dir_path'],
                            proj_comp,
                            aligner,
                            quantifier,
                            "overlap.txt"
                        )
                    ]

                task_cfg['loci_stat'][quantifier][aligner][proj_comp]['shell_script_path'] = \
                    '{}/run.lociStat.{}.{}.{}.sh'.format(
                        shell_dir_path,
                        quantifier,
                        aligner,
                        proj_comp
                    )

                task_cfg['loci_stat'][quantifier][aligner][proj_comp]['log_file_path'] = \
                    '{}/run.lociStat.{}.{}.{}.log'.format(
                        task_cfg['loci_stat'][quantifier][aligner][proj_comp]['out_dir_path'],
                        quantifier,
                        aligner,
                        proj_comp
                    )

                #util.ddictfunc.pprint_ddicts(task_cfg, ['loci_stat'])
                '''
                {'loci_stat': {'featurecounts': {'star': {'DLBC': {'in_file_path_list': ['/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC_full/RNAseq/DEG/edger/featurecounts/star/DLBC.star.featurecounts.edger.test.DEG.txt',
                                                                                         '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC_full/RNAseq/DEG/deseq2/featurecounts/star/DLBC.star.featurecounts.deseq2.test.DEG.txt',
                                                                                         '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC_full/RNAseq/DEG/limma/featurecounts/star/DLBC.star.featurecounts.limma.test.DEG.txt'],
                                                                   'log_file_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/LociStat/featurecounts/star/DLBC/run.call.featurecounts.star.DLBC.log',
                                                                   'out_dir_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/LociStat/featurecounts/star/DLBC',
                                                                   'out_file_path_list': ['/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/LociStat/featurecounts/star/DLBC/DLBC.star.featurecounts.overlap.txt'],
                                                                   'shell_script_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/shell_scripts/run.lociStat.featurecounts.star.DLBC.sh'},
                                                          'main': {'anchor_file_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC_full/RNAseq/DEG/deseq2/featurecounts/star/DLBC.star.featurecounts.deseq2.test.txt'}}},
                               'main': {'application': 'RNAseq'},
                               'references': {'grch38': {'anno_bed': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.bed12',
                                                         'anno_bed12': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.bed12',
                                                         'anno_gtf': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.gtf',
                                                         'anno_refflat': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.refFlat.txt',
                                                         'anno_ribosome_rna_bed': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.rRNA.bed',
                                                         'anno_ribosome_rna_interval': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.rRNA.interval_list',
                                                         'chrom_size': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/GRCh38.primary_assembly.genome.chr11.chrom.size',
                                                         'chrs': 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM',
                                                         'genome': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/GRCh38.primary_assembly.genome.chr11.fa',
                                                         'genomedict': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/GRCh38.primary_assembly.genome.chr11.dict',
                                                         'star_index': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12'}}}}

                '''

                local_resource = ''
                if args.system_type == 'cluster':
                    pbs_ppn = min([project_cfg['pipeline']['software'][caller]['threads'], args.threads])
                    local_resource = ', cpus := {}, mem := {}*G, timeout := {}*hour, taskName := "lociStat.{}.{}.{}"'.format(
                        pbs_ppn,
                        pbs_ppn * 8,
                        72,
                        quantifier,
                        aligner,
                        project_cfg['project']['name']
                    )

                out_file_path = args.submitter
                try:
                    with open(out_file_path, "a") as outfile:
                        outfile.write('''
// Loci Statistics: [ {} ] aligned by [ {} ] and quantified by  [ {} ]

dep( [ '{}' ] <- [ '{}' ]{}) sys bash {}; sleep {}
goal( [ '{}' ] )
\n'''.format(
         project_cfg['project']['name'],
         aligner,
         quantifier,
         "', '".join(task_cfg['loci_stat'][quantifier][aligner][proj_comp]['out_file_path_list']),
         "', '".join(task_cfg['loci_stat'][quantifier][aligner][proj_comp]['in_file_path_list']),
         local_resource,
         task_cfg['loci_stat'][quantifier][aligner][proj_comp]['shell_script_path'],
         project_cfg['pipeline']['software']['bigdatascript']['safeSleep'],
         "', '".join(task_cfg['loci_stat'][quantifier][aligner][proj_comp]['out_file_path_list'])
     ))
                        outfile.close()
                except IOError as exc:
                    print(exc)

                logging.debug(
                    "[ {} ] Loci Statistics:: Project: {} - DONE\n".format(
                        SELF_FILE,
                        proj_comp
                    )
                )

            logging.debug(
                "[ {} ] Loci Statistics:: Aligner: {} - DONE\n".format(
                    SELF_FILE,
                    aligner
                )
            )

        logging.debug(
            "[ {} ] Loci Statistics:: Quantifier: {} - DONE\n".format(
                SELF_FILE,
                quantifier
            )
        )

    out_file_path = args.submitter
    try:
        with open(out_file_path, "a") as outfile:
            outfile.write("\n// <<< Loci Statistics\n\n\n")
            outfile.close()
    except IOError as exc:
        print(exc)

    #util.ddictfunc.pprint_ddicts(task_cfg, ['loci_stat'])

    logging.debug("[ {} ] Loci Statistics - DONE\n".format(SELF_FILE))


    logging.info("[ {} ] Post Analysis - Heat Map & GSEA\n".format(SELF_FILE))

    out_file_path = args.submitter
    try:
        with open(out_file_path, "a") as outfile:
            outfile.write("\n// Post Analysis >>>\n\n\n")
            outfile.close()
    except IOError as exc:
        print(exc)

    if("RNAseq" == proj_appl):
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
        caller = 'deseq2' if 'deseq2' in project_cfg['callers'] else project_cfg['callers'][1]
    else:
        logging.error(
            "[ {} ] cannot recognize [{}]\n".format(
                SELF_FILE,
                proj_appl
            )
        )
        sys.exit()

    task_cfg['post_ana']['references'][target_genome] = project_cfg['pipeline']['references'][target_genome]
    task_cfg['post_ana']['main'].update({'application': proj_appl})

    for quantifier in project_cfg['quantifiers']:
        for aligner in project_cfg['aligners']:
            if True:
                proj_comp = project_cfg['project']['name']

                task_cfg['post_ana']['main'][quantifier][aligner][proj_comp]['in_file_path_list'] = []

                task_cfg['post_ana']['main'][quantifier][aligner][proj_comp]['in_file_path_list'].append(
                    '{}/{}.{}.{}.{}'.format(
                        task_cfg['loci_stat'][quantifier][aligner][proj_comp]['out_dir_path'],
                        proj_comp,
                        aligner,
                        quantifier,
                        "overlap.txt"
                    )
                )

                task_cfg['post_ana']['main'][quantifier][aligner][proj_comp]['expr_tbl_path'] = \
                    '{}/{}.{}.{}.{}.{}'.format(
                        task_cfg['call'][caller][quantifier][aligner]['out_dir_path'],
                        project_cfg['project']['name'],
                        aligner,
                        quantifier,
                        caller,
                        "count.txt"
                    )

                if project_cfg['project']['run_as_practice']:
                    task_cfg['post_ana']['main'][quantifier][aligner][proj_comp]['expr_tbl_path'] = \
                        re.sub(
                            r'DLBC/RNAseq',
                            'DLBC_full/RNAseq',
                            task_cfg['post_ana']['main'][quantifier][aligner][proj_comp]['expr_tbl_path']
                        )
                    logging.debug("[RUN_AS_PRACTICE] useing full set result instead")

    for tool,ext_list in dict_post_ana_tool.items():
        logging.info(
            "[ {} ] Post Analysis:: tool: {}\n".format(
                SELF_FILE,
                tool
            )
        )

        for quantifier in project_cfg['quantifiers']:
            logging.info(
                "[ {} ] Post Analysis:: Quantifier: {}\n".format(
                    SELF_FILE,
                    quantifier
                )
            )

            for aligner in project_cfg['aligners']:
                logging.info(
                    "[ {} ] Post Analysis:: Aligner: {}\n".format(
                        SELF_FILE,
                        aligner
                    )
                )

                if True:
                    proj_comp = project_cfg['project']['name']

                    task_cfg['post_ana'][tool][quantifier][aligner][proj_comp]['out_dir_path'] = \
                        '{}/PostAna/{}/{}/{}/{}'.format(
                            appl_dir_path,
                            tool,
                            quantifier,
                            aligner,
                            proj_comp
                        )

                    out_dir_t_path = task_cfg['post_ana'][tool][quantifier][aligner][proj_comp]['out_dir_path']
                    if not os.path.exists(out_dir_t_path):
                        logging.debug("MKDIR: [ {} ]".format(out_dir_t_path))
                        os.makedirs(out_dir_t_path, exist_ok=True)

                    out_file_base = \
                        '{}/{}.{}.{}'.format(
                            task_cfg['post_ana'][tool][quantifier][aligner][proj_comp]['out_dir_path'],
                            proj_comp,
                            aligner,
                            quantifier
                        )

                    task_cfg['post_ana'][tool][quantifier][aligner][proj_comp]['out_file_path_list'] = \
                            [
                                '{}.{}'.format(
                                    out_file_base,
                                    ext
                                ) for ext in ext_list
                            ]

                    task_cfg['post_ana'][tool][quantifier][aligner][proj_comp]['shell_script_path'] = \
                        '{}/run.postAna.{}.{}.{}.{}.sh'.format(
                            shell_dir_path,
                            tool,
                            quantifier,
                            aligner,
                            proj_comp
                        )

                    task_cfg['post_ana'][tool][quantifier][aligner][proj_comp]['log_file_path'] = \
                        '{}/run.postAna.{}.{}.{}.{}.log'.format(
                            task_cfg['post_ana'][tool][quantifier][aligner][proj_comp]['out_dir_path'],
                            tool,
                            quantifier,
                            aligner,
                            proj_comp
                        )

                    #util.ddictfunc.pprint_ddicts(task_cfg, ['post_ana'])
                    '''
                    {'post_ana': {'main': {'application': 'RNAseq',
                                           'featurecounts': {'star': {'DLBC': {'expr_tbl_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC_full/RNAseq/DEG/deseq2/featurecounts/star/DLBC.star.featurecounts.deseq2.count.txt',
                                                                               'in_file_path_list': ['/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/LociStat/featurecounts/star/DLBC/DLBC.star.featurecounts.overlap.txt']}}}},
                                  'pheatmap': {'featurecounts': {'star': {'DLBC': {'log_file_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/PostAna/pheatmap/featurecounts/star/DLBC/run.postAna.pheatmap.featurecounts.star.DLBC.log',
                                                                                   'out_dir_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/PostAna/pheatmap/featurecounts/star/DLBC',
                                                                                   'out_file_path_list': ['/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/PostAna/pheatmap/featurecounts/star/DLBC/DLBC.star.featurecounts.heatmap.pdf'],
                                                                                   'shell_script_path': '/gpfs/data/bioinformatics/cri_rnaseq_2018/example/DLBC/RNAseq/shell_scripts/run.postAna.pheatmap.featurecounts.star.DLBC.sh'}}}},
                                  'references': {'grch38': {'anno_bed': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.bed12',
                                                            'anno_bed12': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.bed12',
                                                            'anno_gtf': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.gtf',
                                                            'anno_refflat': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.refFlat.txt',
                                                            'anno_ribosome_rna_bed': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.rRNA.bed',
                                                            'anno_ribosome_rna_interval': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/genes.rRNA.interval_list',
                                                            'chrom_size': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/GRCh38.primary_assembly.genome.chr11.chrom.size',
                                                            'chrs': 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM',
                                                            'genome': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/GRCh38.primary_assembly.genome.chr11.fa',
                                                            'genomedict': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12/GRCh38.primary_assembly.genome.chr11.dict',
                                                            'star_index': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/v28_92_GRCh38.p12'}}}}
                    '''


                    local_resource = ''
                    if args.system_type == 'cluster':
                        pbs_ppn = min([project_cfg['pipeline']['software'][caller]['threads'], args.threads])
                        local_resource = ', cpus := {}, mem := {}*G, timeout := {}*hour, taskName := "postAna.{}.{}.{}.{}"'.format(
                            pbs_ppn,
                            pbs_ppn * 8,
                            72,
                            tool,
                            quantifier,
                            aligner,
                            project_cfg['project']['name']
                        )

                    out_file_path = args.submitter
                    try:
                        with open(out_file_path, "a") as outfile:
                            outfile.write('''
// Post Analysis: [ {} ] aligned by [ {} ] and quantified by  [ {} ] using post analysis tool [ {} ]

dep( [ '{}' ] <- [ '{}' ]{}) sys bash {}; sleep {}
goal( [ '{}' ] )
\n'''.format(
         project_cfg['project']['name'],
         aligner,
         quantifier,
         tool,
         "', '".join(task_cfg['post_ana'][tool][quantifier][aligner][proj_comp]['out_file_path_list']),
         "', '".join(task_cfg['post_ana']['main'][quantifier][aligner][proj_comp]['in_file_path_list']),
         local_resource,
         task_cfg['post_ana'][tool][quantifier][aligner][proj_comp]['shell_script_path'],
         project_cfg['pipeline']['software']['bigdatascript']['safeSleep'],
         "', '".join(task_cfg['post_ana'][tool][quantifier][aligner][proj_comp]['out_file_path_list'])
     ))
                            outfile.close()
                    except IOError as exc:
                        print(exc)

                    logging.debug(
                        "[ {} ] Post Analysis:: Project: {} - DONE\n".format(
                            SELF_FILE,
                            proj_comp
                        )
                    )

                logging.debug(
                    "[ {} ] Post Analysis:: Aligner: {} - DONE\n".format(
                        SELF_FILE,
                        aligner
                    )
                )

            logging.debug(
                "[ {} ] Post Analysis:: Quantifier: {} - DONE\n".format(
                    SELF_FILE,
                    quantifier
                )
            )

        logging.debug(
            "[ {} ] Post Analysis:: tool: {} - DONE\n".format(
                SELF_FILE,
                tool
            )
        )

    out_file_path = args.submitter
    try:
        with open(out_file_path, "a") as outfile:
            outfile.write("\n// <<< Post Analysis\n\n\n")
            outfile.close()
    except IOError as exc:
        print(exc)

    #util.ddictfunc.pprint_ddicts(task_cfg, ['post_ana'])

    logging.debug("[ {} ] Post Analysis - Heat Map & GSEA - DONE\n".format(SELF_FILE))


    return
