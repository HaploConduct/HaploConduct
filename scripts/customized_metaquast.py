#!/usr/bin/env python

############################################################################
# Copyright (c) 2015-2016 Saint Petersburg State University
# Copyright (c) 2011-2015 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

from __future__ import with_statement
import sys
import os
import shutil

from libs import qconfig
from libs.metautils import remove_from_quast_py_args, Assembly, correct_meta_references, correct_assemblies, \
    get_downloaded_refs_with_alignments, partition_contigs
from libs.options_parser import parse_options

qconfig.check_python_version()
from libs import contigs_analyzer, reads_analyzer, search_references_meta
from libs import qutils
from libs.qutils import cleanup

from libs.log import get_logger
logger = get_logger(qconfig.LOGGER_META_NAME)
logger.set_up_console_handler()

from site import addsitedir
addsitedir(os.path.join(qconfig.LIBS_LOCATION, 'site_packages'))


def _start_quast_main(
        args, assemblies, reference_fpath=None,
        output_dirpath=None, num_notifications_tuple=None, is_first_run=None):
    args = args[:]

    args.extend([asm.fpath for asm in assemblies])

    if reference_fpath:
        args.append('-R')
        args.append(reference_fpath)

    if output_dirpath:
        args.append('-o')
        args.append(output_dirpath)

    args.append('--labels')

    def quote(line):
        if ' ' in line:
            line = '"%s"' % line
        return line

    args.append(quote(', '.join([asm.label for asm in assemblies])))

    import quast
    reload(quast)
    quast.logger.set_up_console_handler(indent_val=1, debug=qconfig.debug)
    quast.logger.set_up_metaquast()
    logger.info_to_file('(logging to ' +
                        os.path.join(output_dirpath,
                                     qconfig.LOGGER_DEFAULT_NAME + '.log)'))
    return_code = quast.main(args)
    if num_notifications_tuple:
        cur_num_notifications = quast.logger.get_numbers_of_notifications()
        num_notifications_tuple = map(sum, zip(num_notifications_tuple,cur_num_notifications))

    if is_first_run:
        labels = [qconfig.assembly_labels_by_fpath[fpath] for fpath in qconfig.assemblies_fpaths]
        assemblies = [Assembly(fpath, qconfig.assembly_labels_by_fpath[fpath]) for fpath in qconfig.assemblies_fpaths]
        return return_code, num_notifications_tuple, assemblies, labels
    else:
        return return_code, num_notifications_tuple


def main(args):
    if ' ' in qconfig.QUAST_HOME:
        logger.error('QUAST does not support spaces in paths. \n'
                     'You are trying to run it from ' + str(qconfig.QUAST_HOME) + '\n'
                     'Please, put QUAST in a different directory, then try again.\n',
                     to_stderr=True,
                     exit_with_code=3)

    if not args:
        qconfig.usage(meta=True)
        sys.exit(0)

    metaquast_path = [os.path.realpath(__file__)]
    quast_py_args, contigs_fpaths = parse_options(logger, metaquast_path + args, is_metaquast=True)
    output_dirpath, ref_fpaths, labels = qconfig.output_dirpath, qconfig.reference, qconfig.labels
    html_report = qconfig.html_report
    test_mode = qconfig.test

    # Directories
    output_dirpath, _, _ = qutils.set_up_output_dir(
        output_dirpath, None, not output_dirpath,
        save_json=False)

    corrected_dirpath = os.path.join(output_dirpath, qconfig.corrected_dirname)

    qconfig.set_max_threads(logger)
    qutils.logger = logger

    ########################################################################

    from libs import reporting
    reload(reporting)
    from libs import plotter

    if os.path.isdir(corrected_dirpath):
        shutil.rmtree(corrected_dirpath)
    os.mkdir(corrected_dirpath)

    # PROCESSING REFERENCES
    if ref_fpaths:
        logger.main_info()
        logger.main_info('Reference(s):')

        corrected_ref_fpaths, combined_ref_fpath, chromosomes_by_refs, ref_names =\
            correct_meta_references(ref_fpaths, corrected_dirpath)

    # PROCESSING CONTIGS
    logger.main_info()
    logger.main_info('Contigs:')
    assemblies, labels = correct_assemblies(contigs_fpaths, output_dirpath, labels)
    if not assemblies:
        logger.error("None of the assembly files contains correct contigs. "
                     "Please, provide different files or decrease --min-contig threshold.")
        return 4

    # Running QUAST(s)
    quast_py_args += ['--meta']
    downloaded_refs = False

    # SEARCHING REFERENCES
    if not ref_fpaths:
        logger.main_info()
        if qconfig.max_references == 0:
            logger.notice("Maximum number of references (--max-ref-number) is set to 0, search in SILVA 16S rRNA database is disabled")
        else:
            if qconfig.references_txt:
                logger.main_info("List of references was provided, starting to download reference genomes from NCBI...")
            else:
                logger.main_info("No references are provided, starting to search for reference genomes in SILVA 16S rRNA database "
                        "and to download them from NCBI...")
            downloaded_dirpath = os.path.join(output_dirpath, qconfig.downloaded_dirname)
            if not os.path.isdir(downloaded_dirpath):
                os.mkdir(downloaded_dirpath)
            ref_fpaths = search_references_meta.do(assemblies, labels, downloaded_dirpath, qconfig.references_txt)
            if ref_fpaths:
                search_references_meta.is_quast_first_run = True
                if not qconfig.references_txt:
                    downloaded_refs = True
                logger.main_info()
                logger.main_info('Downloaded reference(s):')
                corrected_ref_fpaths, combined_ref_fpath, chromosomes_by_refs, ref_names =\
                    correct_meta_references(ref_fpaths, corrected_dirpath)
            elif test_mode and not ref_fpaths:
                logger.error('Failed to download or setup SILVA 16S rRNA database for working without '
                             'references on metagenome datasets!', to_stderr=True, exit_with_code=4)

    if not ref_fpaths:
        # No references, running regular quast with MetaGenemark gene finder
        logger.main_info()
        logger.notice('No references are provided, starting regular QUAST with MetaGeneMark gene finder')
        _start_quast_main(quast_py_args, assemblies=assemblies, output_dirpath=output_dirpath)
        exit(0)

    # Running combined reference
    combined_output_dirpath = os.path.join(output_dirpath, qconfig.combined_output_name)

    reads_fpaths = []
    if qconfig.forward_reads:
        reads_fpaths.append(qconfig.forward_reads)
    if qconfig.reverse_reads:
        reads_fpaths.append(qconfig.reverse_reads)
    if (reads_fpaths or qconfig.sam or qconfig.bam) and ref_fpaths:
        bed_fpath, cov_fpath, _ = reads_analyzer.do(combined_ref_fpath, contigs_fpaths, reads_fpaths, corrected_ref_fpaths,
                                                    os.path.join(combined_output_dirpath, qconfig.variation_dirname),
                                                    external_logger=logger, sam_fpath=qconfig.sam, bam_fpath=qconfig.bam, bed_fpath=qconfig.bed)
        qconfig.bed = bed_fpath

    if qconfig.bed:
        quast_py_args += ['--sv-bed']
        quast_py_args += [qconfig.bed]
    if qconfig.sam:
        quast_py_args += ['--sam']
        quast_py_args += [qconfig.sam]
    if qconfig.bam:
        quast_py_args += ['--bam']
        quast_py_args += [qconfig.bam]

    for arg in args:
        if arg in ('-s', "--scaffolds"):
            quast_py_args.remove(arg)

    quast_py_args += ['--combined-ref']
    if qconfig.draw_plots or qconfig.html_report:
        if plotter.dict_color_and_ls:
            colors_and_ls = [plotter.dict_color_and_ls[asm.label] for asm in assemblies]
            quast_py_args += ['--colors']
            quast_py_args += [','.join([style[0] for style in colors_and_ls])]
            quast_py_args += ['--ls']
            quast_py_args += [','.join([style[1] for style in colors_and_ls])]
    run_name = 'for the combined reference'
    logger.main_info()
    logger.main_info('Starting quast.py ' + run_name + '...')
    total_num_notices = 0
    total_num_warnings = 0
    total_num_nf_errors = 0
    total_num_notifications = (total_num_notices, total_num_warnings, total_num_nf_errors)
    if qconfig.html_report:
        from libs.html_saver import json_saver
        json_texts = []
    else:
        json_texts = None
    return_code, total_num_notifications, assemblies, labels = \
        _start_quast_main(quast_py_args + ([] if qconfig.unique_mapping else ["--ambiguity-usage", 'one']),
        assemblies=assemblies,
        reference_fpath=combined_ref_fpath,
        output_dirpath=combined_output_dirpath,
        num_notifications_tuple=total_num_notifications, is_first_run=True)

    if json_texts is not None:
        json_texts.append(json_saver.json_text)
    search_references_meta.is_quast_first_run = False

    genome_info_dirpath = os.path.join(output_dirpath, qconfig.combined_output_name, 'genome_stats')
    genome_info_fpath = os.path.join(genome_info_dirpath, 'genome_info.txt')
    if not os.path.exists(genome_info_fpath):
        logger.main_info('')
        logger.main_info('Failed aligning the contigs for all the references. ' + ('Try to restart MetaQUAST with another references.'
                                                        if not downloaded_refs else 'Try to use option --max-ref-number to change maximum number of references '
                                                                                    '(per each assembly) to download.'))
        logger.main_info('')
        cleanup(corrected_dirpath)
        logger.main_info('MetaQUAST finished.')
        logger.finish_up(numbers=tuple(total_num_notifications), check_test=test_mode)
        return

    if downloaded_refs:
        logger.main_info()
        logger.main_info('Excluding downloaded references with low genome fraction from further analysis..')
        corr_ref_fpaths = get_downloaded_refs_with_alignments(genome_info_fpath, ref_fpaths, chromosomes_by_refs)
        if corr_ref_fpaths and corr_ref_fpaths != ref_fpaths:
            logger.main_info()
            logger.main_info('Filtered reference(s):')
            os.remove(combined_ref_fpath)
            contigs_analyzer.ref_labels_by_chromosomes = {}
            corrected_ref_fpaths, combined_ref_fpath, chromosomes_by_refs, ref_names = \
                correct_meta_references(corr_ref_fpaths, corrected_dirpath)
            run_name = 'for the corrected combined reference'
            logger.main_info()
            logger.main_info('Starting quast.py ' + run_name + '...')
            return_code, total_num_notifications, assemblies, labels = \
                _start_quast_main(quast_py_args + ([] if qconfig.unique_mapping else ["--ambiguity-usage", 'one']),
                assemblies=assemblies,
                reference_fpath=combined_ref_fpath,
                output_dirpath=combined_output_dirpath,
                num_notifications_tuple=total_num_notifications, is_first_run=True)
            if json_texts is not None:
                json_texts = json_texts[:-1]
                json_texts.append(json_saver.json_text)
        elif corr_ref_fpaths == ref_fpaths:
            logger.main_info('All downloaded references have genome fraction more than 10%. Nothing was excluded.')
        else:
            logger.main_info('All downloaded references have low genome fraction. Nothing was excluded for now.')

    quast_py_args += ['--no-check-meta']
    qconfig.contig_thresholds = ','.join([str(threshold) for threshold in qconfig.contig_thresholds if threshold > qconfig.min_contig])
    if not qconfig.contig_thresholds:
        qconfig.contig_thresholds = 'None'
    quast_py_args = remove_from_quast_py_args(quast_py_args, '--contig-thresholds', qconfig.contig_thresholds)
    quast_py_args += ['--contig-thresholds']
    quast_py_args += [qconfig.contig_thresholds]
    quast_py_args.remove('--combined-ref')

    logger.main_info()
    logger.main_info('Partitioning contigs into bins aligned to each reference..')

    assemblies_by_reference, not_aligned_assemblies = partition_contigs(
        assemblies, corrected_ref_fpaths, corrected_dirpath,
        os.path.join(combined_output_dirpath, 'contigs_reports', 'alignments_%s.tsv'), labels)

    ref_names = []
    output_dirpath_per_ref = os.path.join(output_dirpath, qconfig.per_ref_dirname)
    for ref_fpath, ref_assemblies in assemblies_by_reference:
        ref_name = qutils.name_from_fpath(ref_fpath)
        logger.main_info('')
        if not ref_assemblies:
            logger.main_info('No contigs were aligned to the reference ' + ref_name + ', skipping..')
        else:
            ref_names.append(ref_name)
            run_name = 'for the contigs aligned to ' + ref_name
            logger.main_info('Starting quast.py ' + run_name)

            return_code, total_num_notifications = _start_quast_main(quast_py_args,
                assemblies=ref_assemblies,
                reference_fpath=ref_fpath,
                output_dirpath=os.path.join(output_dirpath_per_ref, ref_name),
                num_notifications_tuple=total_num_notifications)
            if json_texts is not None:
                json_texts.append(json_saver.json_text)

    # Finally running for the contigs that has not been aligned to any reference
    no_unaligned_contigs = True
    for assembly in not_aligned_assemblies:
        if os.path.isfile(assembly.fpath) and os.stat(assembly.fpath).st_size != 0:
            no_unaligned_contigs = False
            break

    run_name = 'for the contigs not aligned anywhere'
    logger.main_info()
    if no_unaligned_contigs:
        logger.main_info('Skipping quast.py ' + run_name + ' (everything is aligned!)')
    else:
        logger.main_info('Starting quast.py ' + run_name + '...')

        return_code, total_num_notifications = _start_quast_main(quast_py_args,
            assemblies=not_aligned_assemblies,
            output_dirpath=os.path.join(output_dirpath, qconfig.not_aligned_name),
            num_notifications_tuple=total_num_notifications)

        if return_code not in [0, 4]:
            logger.error('Error running quast.py for the contigs not aligned anywhere')
        elif return_code == 4:  # no unaligned contigs, i.e. everything aligned
            no_unaligned_contigs = True
        if not no_unaligned_contigs:
            if json_texts is not None:
                json_texts.append(json_saver.json_text)

    if ref_names:
        logger.print_timestamp()
        logger.main_info("Summarizing results...")

        summary_output_dirpath = os.path.join(output_dirpath, qconfig.meta_summary_dir)
        if not os.path.isdir(summary_output_dirpath):
            os.makedirs(summary_output_dirpath)
        if html_report and json_texts:
            from libs.html_saver import html_saver
            html_summary_report_fpath = html_saver.init_meta_report(output_dirpath)
        else:
            html_summary_report_fpath = None
        from libs import create_meta_summary
        metrics_for_plots = reporting.Fields.main_metrics
        misassembl_metrics = [reporting.Fields.MIS_RELOCATION, reporting.Fields.MIS_TRANSLOCATION, reporting.Fields.MIS_INVERTION,
                           reporting.Fields.MIS_ISTRANSLOCATIONS]
        create_meta_summary.do(html_summary_report_fpath, summary_output_dirpath, combined_output_dirpath, output_dirpath_per_ref, metrics_for_plots, misassembl_metrics,
                               ref_names if no_unaligned_contigs else ref_names + [qconfig.not_aligned_name])
        if html_report and json_texts:
            html_saver.save_colors(output_dirpath, contigs_fpaths, plotter.dict_color_and_ls, meta=True)
            if qconfig.create_icarus_html:
                icarus_html_fpath = html_saver.create_meta_icarus(output_dirpath, ref_names)
                logger.main_info('  Icarus (contig browser) is saved to %s' % icarus_html_fpath)
            html_saver.create_meta_report(output_dirpath, json_texts)

    cleanup(corrected_dirpath)
    logger.main_info('')
    logger.main_info('MetaQUAST finished.')
    return logger.finish_up(numbers=tuple(total_num_notifications), check_test=test_mode)


if __name__ == '__main__':
    try:
        return_code = main(sys.argv[1:])
        exit(return_code)
    except Exception:
        _, exc_value, _ = sys.exc_info()
        logger.exception(exc_value)
        logger.error('exception caught!', exit_with_code=1, to_stderr=True)


















