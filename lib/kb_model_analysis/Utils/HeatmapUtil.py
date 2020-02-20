
import errno
import os
import uuid
import logging
import pandas as pd
from xlrd.biffh import XLRDError
import json
import shutil
from sklearn import preprocessing
import traceback
import sys

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.fba_toolsClient import fba_tools


class HeatmapUtil:

    def _mkdir_p(self, path):
        """
        _mkdir_p: make directory for given path
        """
        if not path:
            return
        try:
            os.makedirs(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise

    def _compute_cluster_label_order(self, values, labels):

        # values = [[1, 0, 21, 50, 1], [20, 0, 60, 80, 30], [30, 60, 1, -10, 20]]
        # labels = ['model_1', 'model_2', 'model_3']
        if len(labels) == 1:
            return labels
        dist_matrix = pdist(values)
        linkage_matrix = linkage(dist_matrix, 'ward')

        dn = dendrogram(linkage_matrix, labels=labels, distance_sort='ascending')

        ordered_label = dn['ivl']

        return ordered_label

    def _read_csv_file(self, file_path):
        try:
            df = pd.read_excel(file_path, dtype='str',  index_col=0)
        except XLRDError:
            reader = pd.read_csv(file_path, sep=None, iterator=True)
            inferred_sep = reader._engine.data.dialect.delimiter
            df = pd.read_csv(file_path, sep=inferred_sep, index_col=0)

        df.fillna('', inplace=True)

        return df

    def _get_model_overall_stats(self, model_obj, model_df, model_ref):
        """
        overall_data: [model_name, len(modelreactions), genecount, len(modelcompounds),
                       core_gapfilling, baseline_gapfilling, base_atp, auxo_biomass,
                       [auxotrophs with 1]]
        """
        overall_stats = list()
        reaction_stats = list()

        model_data = model_obj['data']
        model_info = model_obj['info']
        model_name = model_info[1]

        attributes = model_data.get('attributes', {})
        modelreactions = model_data.get('modelreactions', [])
        modelcompounds = model_data.get('modelcompounds', [])

        fbas = attributes.get('fbas', {})
        auxo_biomass = fbas.get('auxomedia', {}).get('biomass', 0)

        # build overall_stats
        overall_stats.append(model_name)  # model_name
        overall_stats.append(len(modelreactions))  # Total reactions
        overall_stats.append(attributes.get('gene_count', 0))  # Total genes
        overall_stats.append(len(modelcompounds))  # Total compounds
        overall_stats.append(attributes.get('core_gapfilling', 0))  # Core gapfilling
        overall_stats.append(attributes.get('base_gapfilling', 0))  # Rich media gapfilling
        overall_stats.append(attributes.get('base_atp', 0))  # ATP per mol glucose
        overall_stats.append(auxo_biomass)  # Rich biomass yield

        auxotrophy = attributes.get('auxotrophy', {})
        auxotrophy_count = 0
        for auxo_id, auxo_data in auxotrophy.items():
            if auxo_data.get('is_auxotrophic', 0):
                auxotrophy_count += 1

        overall_stats.append(auxotrophy_count)  # Predicted auxotrophies

        meta_data = model_df.loc[model_ref, :].tolist()
        overall_stats.extend(meta_data)

        # build reaction_stats
        reaction_stats.append(model_name)  # model_name

        fbas_auxo = fbas.get('auxomedia', {})
        auxo_class_Negative = fbas_auxo.get('Negative', 0)
        auxo_class_Positive = fbas_auxo.get('Positive', 0)

        reaction_stats.append(auxo_class_Negative + auxo_class_Positive)  # Defined media essential
        reaction_stats.append(fbas_auxo.get('NegativeVariable', 0) +
                              fbas_auxo.get('PositiveVariable', 0) +
                              fbas_auxo.get('Variable', 0))  # Defined media functional
        reaction_stats.append(fbas_auxo.get('Blocked', 0))  # Defined media blocked

        fbas_complete = fbas.get('complete', {})
        complete_class_Negative = fbas_complete.get('Negative', 0)
        complete_class_Positive = fbas_complete.get('Positive', 0)

        reaction_stats.append(complete_class_Negative + complete_class_Positive)  # Rich media essential
        reaction_stats.append(fbas_complete.get('NegativeVariable', 0) +
                              fbas_complete.get('PositiveVariable', 0) +
                              fbas_complete.get('Variable', 0))  # Rich media functional
        reaction_stats.append(fbas_complete.get('Blocked', 0))  # Rich media blocked

        return overall_stats, reaction_stats

    def _check_model_obj_version(self, model_df, workspace_name):
        model_refs = model_df.index.tolist()

        for model_ref in model_refs:
            model_obj = self.dfu.get_objects({'object_refs': [model_ref]})['data'][0]
            model_info = model_obj['info']
            model_data = model_obj['data']

            model_type = model_info[2]

            if 'KBaseFBA.FBAModel' not in model_type:
                raise ValueError("Please provide only KBaseFBA.FBAModel objects")

            obj_version = float(model_type.split('-')[-1])

            if 'ci.kbase' in self.endpoint:
                latest_version = 14
            elif 'appdev.kbase' in self.endpoint:
                latest_version = 12
            else:
                latest_version = 12

            if obj_version < latest_version:
                err_msg = "Please provde KBaseFBA.FBAModel with version greater than {}".format(
                                                                                    latest_version)
                raise ValueError(err_msg)

            if not model_data.get('attributes', {}).get('pathways', {}):
                try:
                    logging.warning('Found empty attributes and pathways in {}'.format(model_ref))
                    logging.warning('Trying to run model characterization')
                    ret = self.fba_tools.run_model_characterization({
                                                            'fbamodel_id': model_info[1],
                                                            'workspace': workspace_name,
                                                            'fbamodel_output_id': model_info[1]})
                    logging.warning('Generated new objects: {}'.format(ret))
                    # new_ref = ret.get('new_fbamodel_ref')

                    # idx = model_df.index.values.tolist().index(model_ref)
                    # model_df.index.values[idx] = '/'.join(new_ref.split('/')[:2])
                except Exception:
                    logging.warning('failed to run run_model_characterization')
                    logging.warning(traceback.format_exc())
                    logging.warning(sys.exc_info()[2])

            idx = model_df.index.values.tolist().index(model_ref)
            model_df.index.values[idx] = '/'.join(model_ref.split('/')[:2])

    def _build_model_comparison_data(self, model_df):

        logging.info('Start building overall and reaction statistics')

        model_refs = model_df.index.tolist()
        overall_stats = list()
        reaction_stats = list()

        for model_ref in model_refs:
            model_obj = self.dfu.get_objects({'object_refs': [model_ref]})['data'][0]
            (model_overall_stats,
             model_reaction_stats) = self._get_model_overall_stats(model_obj, model_df, model_ref)
            overall_stats.append(model_overall_stats)
            reaction_stats.append(model_reaction_stats)

        return overall_stats, reaction_stats

    def _normalize_data(self, df, normalization_type):

        logging.info('Start normalizing data')

        nor_df = df.copy(deep=True)
        nor_df.fillna(0, inplace=True)

        if normalization_type == 'zscore':
            X_scaled = preprocessing.scale(nor_df.values, axis=1)
            nor_df = pd.DataFrame(index=nor_df.index, columns=nor_df.columns, data=X_scaled)
        elif normalization_type == 'rownormalization':
            idx = nor_df.index.tolist()
            nor_values = list()
            for key in idx:
                line_values = nor_df.loc[key, :].tolist()
                max_value = max(line_values)

                try:
                    nor_values.append([i / float(max_value) for i in line_values])
                except Exception:
                    logging.warning('Cannot compute row normalize on\n{}: {}'.format(key,
                                                                                     line_values))
                    nor_values.append(line_values)

            nor_df = pd.DataFrame(index=nor_df.index, columns=nor_df.columns, data=nor_values)
        else:
            logging.warning('Unexpected normalization type: {}'.format(normalization_type))

        return nor_df

    def _get_pathway_heatmap_data(self, field_type, model_refs):

        logging.info('Start building pathway heatmap data for {}'.format(field_type))

        nor_type = field_type.split('_')[-1]
        if nor_type not in ['zscore', 'rownormalization', 'dividepathwaysize']:
            nor_type = None

        if nor_type:
            field_type = '_'.join(field_type.split('_')[:-1])

        first_model_ref = model_refs[0]
        model_obj = self.dfu.get_objects({'object_refs': [first_model_ref]})['data'][0]
        model_data = model_obj['data']
        model_info = model_obj['info']
        attributes = model_data.get('attributes', {})
        pathways = attributes.get('pathways', {})

        pathway_ids = list()
        pathway_names = list()
        pathway_class2_names = list()
        fetched_pathway_value = list()
        pathway_count = 0
        for pathway_id, pathway_data in pathways.items():
            pathway_ids.append(pathway_id)
            pathway_name = pathway_data.get('name') + ' [{}]'.format(pathway_count)
            classes = pathway_data.get('classes', [])
            pathway_class2 = classes[1] + ' [{}]'.format(pathway_count)
            pathway_names.append(pathway_name)
            pathway_class2_names.append(pathway_class2)
            pathway_value = pathway_data.get(field_type, 0)
            if nor_type == 'dividepathwaysize':
                pathway_size = pathway_data.get('pathway_size', 1)
                if pathway_size <= 0:
                    pathway_size = 1
                pathway_value = pathway_value / float(pathway_size)
            fetched_pathway_value.append(pathway_value)
            pathway_count += 1

        model_name = model_info[1] + ' [0]'
        pathway_df = pd.DataFrame({model_name: fetched_pathway_value}, index=pathway_names)

        for model_ref in model_refs[1:]:
            model_obj = self.dfu.get_objects({'object_refs': [model_ref]})['data'][0]
            model_data = model_obj['data']
            model_info = model_obj['info']
            attributes = model_data.get('attributes', {})
            pathways = attributes.get('pathways', {})
            # TODO: should check/add newly found pathways found in models
            fetched_pathway_value = list()
            for pathway_id in pathway_ids:
                pathway_data = pathways.get(pathway_id, {})
                pathway_value = pathway_data.get(field_type, 0)
                if nor_type == 'dividepathwaysize':
                    pathway_size = pathway_data.get('pathway_size', 1)
                    if pathway_size <= 0:
                        pathway_size = 1
                    pathway_value = pathway_value / float(pathway_size)
                fetched_pathway_value.append(pathway_value)

            model_name = model_info[1] + ' [{}]'.format(model_refs.index(model_ref))
            pathway_df[model_name] = fetched_pathway_value

        if nor_type:
            pathway_df = self._normalize_data(pathway_df, nor_type)

        col_ordered_label = self._compute_cluster_label_order(pathway_df.T.values.tolist(),
                                                              pathway_df.T.index.tolist())

        idx_ordered_label = self._compute_cluster_label_order(pathway_df.values.tolist(),
                                                              pathway_df.index.tolist())

        pathway_df = pathway_df.reindex(index=idx_ordered_label, columns=col_ordered_label)

        class_ordered_label = list()
        for label in pathway_df.index.tolist():
            class_ordered_label.append(pathway_class2_names[pathway_names.index(label)])

        pathway_info = {'name': pathway_df.index.tolist(),
                        'class': class_ordered_label}

        return pathway_df, pathway_info

    def _build_heatmap_meta(self, model_df):

        logging.info('Start building heatmap metadata')

        model_refs = model_df.index.tolist()

        meta_data_names = model_df.columns.tolist()

        model_names = list()
        for i, model_ref in enumerate(model_refs):
            model_obj = self.dfu.get_objects({'object_refs': [model_ref]})['data'][0]
            model_name = model_obj['info'][1]
            model_names.append(model_name + ' [{}]'.format(i))

        heatmap_meta = dict()
        for meta_data_name in meta_data_names:
            model_meta_mapping = dict()
            model_meta = model_df.loc[:, meta_data_name].tolist()

            model_meta_data = list()
            for i, meta_data in enumerate(model_meta):
                model_meta_data.append(meta_data + '[{}]'.format(i))

            for i, model_name in enumerate(model_names):
                model_meta_mapping.update({model_name: model_meta_data[i]})

            heatmap_meta.update({meta_data_name: model_meta_mapping})

        return heatmap_meta

    def _build_heatmap_data(self, model_df):

        logging.info('Start building heatmap data')

        model_refs = model_df.index.tolist()

        heatmap_data = dict()

        pathway_types = ['functional_rxn', 'functional_rxn_zscore',
                         'functional_rxn_rownormalization', 'functional_rxn_dividepathwaysize',
                         'gapfilled_rxn', 'gapfilled_rxn_zscore',
                         'gapfilled_rxn_rownormalization', 'gapfilled_rxn_dividepathwaysize',
                         'nonfunctional_rxn', 'nonfunctional_rxn_zscore',
                         'nonfunctional_rxn_rownormalization', 'nonfunctional_rxn_dividepathwaysize',
                         'gene_count', 'gene_count_zscore',
                         'gene_count_rownormalization', 'gene_count_dividepathwaysize',
                         'average_genes_per_reaction', 'stddev_genes_per_reaction',
                         'average_coverage_per_reaction', 'stddev_coverage_per_reaction']

        pathway_name_map = {
            'functional_rxn': 'Functional Reaction',
            'functional_rxn_zscore': 'Functional Reaction (Z-Score)',
            'functional_rxn_rownormalization': 'Functional Reaction (Divide Row Maximum Value)',
            'functional_rxn_dividepathwaysize': 'Functional Reaction (Divide Pathway Size)',
            'gapfilled_rxn': 'Gapfilled Reaction',
            'gapfilled_rxn_zscore': 'Gapfilled Reaction (Z-Score)',
            'gapfilled_rxn_rownormalization': 'Gapfilled Reaction (Divide Row Maximum Value)',
            'gapfilled_rxn_dividepathwaysize': 'Gapfilled Reaction (Divide Pathway Size)',
            'nonfunctional_rxn': 'Nonfunctional Reaction',
            'nonfunctional_rxn_zscore': 'Nonfunctional Reaction (Z-Score)',
            'nonfunctional_rxn_rownormalization': 'Nonfunctional Reaction (Divide Row Maximum Value)',
            'nonfunctional_rxn_dividepathwaysize': 'Nonfunctional Reaction (Divide Pathway Size)',
            'gene_count': 'Gene Count',
            'gene_count_zscore': 'Gene Count (Z-Score)',
            'gene_count_rownormalization': 'Gene Count (Divide Row Maximum Value)',
            'gene_count_dividepathwaysize': 'Gene Count (Divide Pathway Size)',
            'average_genes_per_reaction': 'Average Genes Per Reaction',
            'stddev_genes_per_reaction': 'Stddev Genes Per Reaction',
            'average_coverage_per_reaction': 'Average Coverage Per Reaction',
            'stddev_coverage_per_reaction': 'Stddev Coverage Per Reaction'}

        for pathway_type in pathway_types:
            pathway_df, pathway_info = self._get_pathway_heatmap_data(pathway_type, model_refs)
            pathway_name = pathway_name_map.get(pathway_type, pathway_type)
            heatmap_data.update({pathway_name: {'values': pathway_df.values.tolist(),
                                                'compound_names': pathway_df.columns.tolist(),
                                                'pathways': pathway_info}})

        return heatmap_data

    def _generate_heatmap_report(self, overall_stats, reaction_stats, heatmap_data, heatmap_meta):

        logging.info('Start building html report')

        html_report = list()

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        logging.info('Start building report files in dir: {}'.format(output_directory))
        self._mkdir_p(output_directory)
        model_comp_file_path = os.path.join(output_directory, 'model_comp.html')

        # build model_comp html (home page)
        shutil.copyfile(os.path.join(os.path.dirname(__file__),
                                     'templates', 'model_comp_template.html'),
                        model_comp_file_path)

        # build overall stats json and html
        logging.info('Start building overall tab')
        overall_stats_json = os.path.join(output_directory, 'overall_stats.json')
        total_records = len(overall_stats)
        overall_stats_data = {'draw': 1, 'recordsTotal': total_records,
                              'recordsFiltered': total_records,
                              'data': overall_stats}
        with open(overall_stats_json, 'w') as outfile:
            json.dump(overall_stats_data, outfile)

        metadata_info = ""
        for meta_name in heatmap_meta.keys():
            metadata_info += """\n<th class="metadata">{}</th>\n""".format(meta_name)

        overall_stats_html = os.path.join(output_directory, 'overall_table.html')
        with open(overall_stats_html, 'w') as overall_stats_html:
            with open(os.path.join(os.path.dirname(__file__),
                                   'templates', 'overall_table_template.html'),
                      'r') as overall_template_file:
                overall_template = overall_template_file.read()
                overall_template = overall_template.replace('deferLoading_count',
                                                            str(total_records))
                overall_template = overall_template.replace('<!-- metadata_info -->',
                                                            metadata_info)
                overall_stats_html.write(overall_template)

        # build reaction stats json and html
        logging.info('Start building reaction tab')
        reaction_stats_json = os.path.join(output_directory, 'reaction_stats.json')
        total_records = len(reaction_stats)
        reaction_stats_data = {'draw': 1, 'recordsTotal': total_records,
                               'recordsFiltered': total_records,
                               'data': reaction_stats}
        with open(reaction_stats_json, 'w') as outfile:
            json.dump(reaction_stats_data, outfile)

        reaction_stats_html = os.path.join(output_directory, 'reaction_table.html')
        with open(reaction_stats_html, 'w') as reaction_stats_html:
            with open(os.path.join(os.path.dirname(__file__),
                                   'templates', 'reaction_table_template.html'),
                      'r') as reaction_template_file:
                reaction_template = reaction_template_file.read()
                reaction_template = reaction_template.replace('deferLoading_count',
                                                              str(total_records))
                reaction_stats_html.write(reaction_template)

        # build heatmap stats json and html
        logging.info('Start building heatmap tab')
        heatmap_data_json = os.path.join(output_directory, 'heatmap_data.json')
        heatmap_meta_json = os.path.join(output_directory, 'heatmap_meta.json')
        with open(heatmap_data_json, 'w') as outfile:
            json.dump(heatmap_data, outfile)
        with open(heatmap_meta_json, 'w') as outfile:
            json.dump(heatmap_meta, outfile)

        metadata_info = ""
        for meta_name in heatmap_meta.keys():
            metadata_info += """\n<option value="{}">{}</option>\n""".format(meta_name, meta_name)

        data_info = ""
        for data_name in heatmap_data.keys():
            data_info += """\n<option value="{}">{}</option>\n""".format(data_name, data_name)

        pathway_info = ""
        for pathway_name in heatmap_data.get(list(heatmap_data.keys())[0])['pathways'].keys():
            pathway_info += """\n<option value="{}">{}</option>\n""".format(pathway_name, pathway_name)

        heatmap_html = os.path.join(output_directory, 'heatmap.html')
        with open(heatmap_html, 'w') as heatmap_html:
            with open(os.path.join(os.path.dirname(__file__),
                                   'templates', 'heatmap_template.html'),
                      'r') as heatmap_template_file:
                heatmap_template = heatmap_template_file.read()
                heatmap_template = heatmap_template.replace('<!-- metadata_info -->',
                                                            metadata_info)
                heatmap_template = heatmap_template.replace('<!-- data_info -->',
                                                            data_info)
                heatmap_template = heatmap_template.replace('<!-- pathway_info -->',
                                                            pathway_info)
                heatmap_html.write(heatmap_template)

        report_shock_id = self.dfu.file_to_shock({'file_path': output_directory,
                                                  'pack': 'zip'})['shock_id']

        html_report.append({'shock_id': report_shock_id,
                            'name': os.path.basename(model_comp_file_path),
                            'label': os.path.basename(model_comp_file_path),
                            'description': 'HTML summary report for Model Comparison App'
                            })

        return html_report

    def __init__(self, config):
        self.callback_url = config['SDK_CALLBACK_URL']
        self.endpoint = config['kbase-endpoint']
        self.scratch = config['scratch']
        self.token = config['KB_AUTH_TOKEN']
        self.dfu = DataFileUtil(self.callback_url)
        self.fba_tools = fba_tools(self.callback_url)
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)

    def run_model_heatmap_analysis(self, params):
        staging_file_path = params.get('staging_file_path')
        workspace_name = params.get('workspace_name')
        attri_mapping_ref = params.get('attri_mapping_ref')

        if staging_file_path:
            model_file_path = self.dfu.download_staging_file(
                            {'staging_file_subdir_path': staging_file_path}).get('copy_file_path')

            model_df = self._read_csv_file(model_file_path)
        elif attri_mapping_ref:
            attri_mapping_data = self.dfu.get_objects(
                                        {'object_refs': [attri_mapping_ref]})['data'][0]['data']
            attributes = pd.DataFrame(attri_mapping_data['attributes'])
            instances = pd.DataFrame(attri_mapping_data['instances'])
            am_df = attributes.join(instances)
            model_name_idx = am_df['attribute'].tolist().index('model_name')
            am_df.drop(columns=['source'], index=[model_name_idx], inplace=True)
            am_df = am_df.T
            am_df.columns = am_df.loc['attribute', :]
            am_df.drop(index=['attribute'], inplace=True)

            model_df = am_df
        else:
            raise ValueError("Please provide valide staging file or attribute mapping")

        self._check_model_obj_version(model_df, workspace_name)

        try:
            model_df.drop(columns=['model_name'], inplace=True)
        except KeyError:
            logging.info('model_name does not exist in excel')

        try:
            model_df.drop(columns=['model_index'], inplace=True)
        except KeyError:
            logging.info('model_index does not exist in excel')

        logging.info('start building stats on {}'.format(model_df.index.tolist()))

        (overall_stats, reaction_stats) = self._build_model_comparison_data(model_df)
        heatmap_data = self._build_heatmap_data(model_df)
        heatmap_meta = self._build_heatmap_meta(model_df)

        html_files = self._generate_heatmap_report(overall_stats, reaction_stats, heatmap_data,
                                                   heatmap_meta)

        report_params = {'message': '',
                         'workspace_name': workspace_name,
                         'html_links': html_files,
                         'direct_html_link_index': 0,
                         'html_window_height': 666,
                         'report_object_name': 'model_comparison_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url, token=self.token)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output
