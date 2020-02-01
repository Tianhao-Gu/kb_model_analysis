
import errno
import os
import uuid
import logging
import pandas as pd
from xlrd.biffh import XLRDError
import json
import shutil

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport


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
        dist_matrix = pdist(values)
        linkage_matrix = linkage(dist_matrix)

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

        genecount = 0
        for modelreaction in modelreactions:
            genecount += modelreaction.get('gene_count', 0)

        # build overall_stats
        overall_stats.append(model_name)  # model_name
        overall_stats.append(len(modelreactions))  # Total reactions
        overall_stats.append(genecount)  # Total genes
        overall_stats.append(len(modelcompounds))  # Total compounds
        overall_stats.append(attributes.get('core_gapfilling', 0))  # Core gapfilling
        overall_stats.append(attributes.get('baseline_gapfilling', 0))  # Rich media gapfilling
        overall_stats.append(attributes.get('base_atp', 0))  # ATP per mol glucose
        overall_stats.append(attributes.get('auxo_biomass', 0))  # Rich biomass yield

        attr_keys = attributes.keys()
        auxotrophy_keys = [s for s in attr_keys if ('auxotrophy_' in s)]
        auxotrophy_keys.remove('auxotrophy_gapfilling')
        auxotrophy_count = 0
        for auxotrophy_key in auxotrophy_keys:
            auxotrophy = attributes.get(auxotrophy_key)
            if auxotrophy.endswith('1'):
                auxotrophy_count += 1
        overall_stats.append(auxotrophy_count)  # Predicted auxotrophies

        meta_data = model_df.loc[model_ref, :].tolist()
        overall_stats.extend(meta_data)

        # build reaction_stats
        reaction_stats.append(model_name)  # model_name
        auxo_class_Negative = attributes.get('auxo_class_Negative', 0)
        auxo_class_Positive = attributes.get('auxo_class_Positive', 0)

        reaction_stats.append(auxo_class_Negative + auxo_class_Positive)  # Defined media essential
        reaction_stats.append(attributes.get('auxo_class_Negative variable', 0) +
                              attributes.get('auxo_class_Positive variable', 0) +
                              attributes.get('auxo_class_Variable', 0))  # Defined media functional
        reaction_stats.append(attributes.get('auxo_class_Blocked', 0))  # Defined media blocked

        complete_class_Negative = attributes.get('complete_class_Negative', 0)
        complete_class_Positive = attributes.get('complete_class_Positive', 0)

        reaction_stats.append(complete_class_Negative + complete_class_Positive)  # Rich media essential
        reaction_stats.append(attributes.get('complete_class_Negative variable', 0) +
                              attributes.get('complete_class_Positive variable', 0) +
                              attributes.get('complete_class_Variable', 0))  # Rich media functional
        reaction_stats.append(attributes.get('complete_class_Blocked', 0))  # Rich media blocked

        return overall_stats, reaction_stats

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

    def _get_pathway_heatmap_data(self, field_type, model_refs):

        logging.info('Start building pathway heatmap data for {}'.format(field_type))

        first_model_ref = model_refs[0]
        model_obj = self.dfu.get_objects({'object_refs': [first_model_ref]})['data'][0]
        model_data = model_obj['data']
        model_info = model_obj['info']
        attributes = model_data.get('attributes', {})
        attr_keys = attributes.keys()
        pathway_keys = [s for s in attr_keys if ('pathways_' in s)]

        pathway_names = list()
        fetched_pathway_value = list()
        for pathway_key in pathway_keys:
            pathway_data = attributes.get(pathway_key)
            pathway_names.append(pathway_key.split('pathways_')[-1])   # TODO: fetch pathway_name from pathway_key
            if field_type == 'gf':
                fetched_pathway_value.append(pathway_data.get('gf', 0))
            elif field_type == 'nonblocked':
                fetched_pathway_value.append(pathway_data.get('nonblocked', 0))
            elif field_type == 'rxn':
                fetched_pathway_value.append(pathway_data.get('rxn', 0))
            else:
                raise ValueError("Unexpected field type")
        model_name = model_info[1]
        pathway_df = pd.DataFrame({model_name: fetched_pathway_value}, index=pathway_names)

        for model_ref in model_refs[1:]:
            model_obj = self.dfu.get_objects({'object_refs': [model_ref]})['data'][0]
            model_data = model_obj['data']
            model_info = model_obj['info']
            attributes = model_data.get('attributes', {})

            fetched_pathway_value = list()
            for pathway_key in pathway_keys:
                pathway_data = attributes.get(pathway_key, {})

                if field_type == 'gf':
                    fetched_pathway_value.append(pathway_data.get('gf', 0))
                elif field_type == 'nonblocked':
                    fetched_pathway_value.append(pathway_data.get('nonblocked', 0))
                elif field_type == 'rxn':
                    fetched_pathway_value.append(pathway_data.get('rxn', 0))
                else:
                    raise ValueError("Unexpected field type")
            # TODO: should check/add newly found pathways found in models
            pathway_df[model_info[1]] = fetched_pathway_value

        ordered_label = self._compute_cluster_label_order(pathway_df.T.values.tolist(),
                                                          pathway_df.T.index.tolist())
        pathway_df = pathway_df.reindex(columns=ordered_label)

        return pathway_df

    def _build_heatmap_meta(self, model_df):

        logging.info('Start building heatmap metadata')

        model_refs = model_df.index.tolist()

        meta_data_names = model_df.columns.tolist()

        model_names = list()
        for model_ref in model_refs:
            model_obj = self.dfu.get_objects({'object_refs': [model_ref]})['data'][0]
            model_name = model_obj['info'][1]
            model_names.append(model_name)

        heatmap_meta = dict()
        for meta_data_name in meta_data_names:
            model_meta_mapping = dict()
            model_meta_data = model_df.loc[:, meta_data_name].tolist()

            for i, model_name in enumerate(model_names):
                model_meta_mapping.update({model_name: model_meta_data[i]})

            heatmap_meta.update({meta_data_name: model_meta_mapping})

        return heatmap_meta

    def _build_heatmap_data(self, model_df):

        logging.info('Start building heatmap data')

        model_refs = model_df.index.tolist()

        heatmap_data = dict()

        # build gf
        gf_df = self._get_pathway_heatmap_data('gf', model_refs)
        heatmap_data.update({'gf': {'values': gf_df.values.tolist(),
                                    'compound_names': gf_df.columns.tolist()}})

        # build nonblocked
        nonblocked_df = self._get_pathway_heatmap_data('nonblocked', model_refs)
        heatmap_data.update({'nonblocked': {'values': nonblocked_df.values.tolist(),
                                            'compound_names': nonblocked_df.columns.tolist()}})

        # build rxn
        rxn_df = self._get_pathway_heatmap_data('rxn', model_refs)
        heatmap_data.update({'rxn': {'values': rxn_df.values.tolist(),
                                     'compound_names': rxn_df.columns.tolist()}})

        pathway_names = gf_df.index.tolist()
        heatmap_data.update({'pathway_names': pathway_names})

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

        heatmap_html = os.path.join(output_directory, 'heatmap.html')
        with open(heatmap_html, 'w') as heatmap_html:
            with open(os.path.join(os.path.dirname(__file__),
                                   'templates', 'heatmap_template.html'),
                      'r') as heatmap_template_file:
                heatmap_template = heatmap_template_file.read()
                heatmap_template = heatmap_template.replace('<!-- metadata_info -->',
                                                            metadata_info)
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
        self.scratch = config['scratch']
        self.token = config['KB_AUTH_TOKEN']
        self.dfu = DataFileUtil(self.callback_url)
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)

    def run_model_heatmap_analysis(self, params):
        staging_file_path = params.get('staging_file_path')
        workspace_name = params.get('workspace_name')

        model_file_path = self.dfu.download_staging_file(
                        {'staging_file_subdir_path': staging_file_path}).get('copy_file_path')

        model_df = self._read_csv_file(model_file_path)

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
