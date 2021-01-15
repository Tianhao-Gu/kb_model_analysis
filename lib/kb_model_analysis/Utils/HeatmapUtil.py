
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
import copy
import math

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist
import plotly.figure_factory as ff
import plotly.graph_objects as go
from plotly.offline import plot

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WsLargeDataIOClient import WsLargeDataIO
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

    def _get_model_obj(self, model_ref):
        cached_obj = self.obj_cache.get(model_ref)

        if cached_obj:
            logging.info('Getting object from cache')
            return cached_obj
        else:
            logging.info('Getting object using DataFileUtil')
            model_obj = self.dfu.get_objects({'object_refs': [model_ref]})['data'][0]
            self.obj_cache[model_ref] = model_obj
            return model_obj

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
        overall_stats.append(attributes.get('baseline_gapfilling', 0))  # Rich media gapfilling
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

        reaction_stats.append(complete_class_Negative +
                              complete_class_Positive)  # Rich media essential
        reaction_stats.append(fbas_complete.get('NegativeVariable', 0) +
                              fbas_complete.get('PositiveVariable', 0) +
                              fbas_complete.get('Variable', 0))  # Rich media functional
        reaction_stats.append(fbas_complete.get('Blocked', 0))  # Rich media blocked

        return overall_stats, reaction_stats

    def _run_model_characterization(self, workspace_name, model_name):
        '''
        run fba_tools.run_model_characterization to generate pathway information

        '''
        try:
            logging.warning('Found empty attributes and pathways in {}'.format(model_name))
            logging.warning('Trying to run model characterization')

            ret = self.fba_tools.run_model_characterization({
                'fbamodel_id': model_name,
                'workspace': workspace_name,
                'fbamodel_output_id': model_name})
            logging.warning('Generated new objects: {}'.format(ret))
            new_model_ref = ret.get('new_fbamodel_ref')

        except Exception:
            new_model_ref = None
            logging.warning('failed to run run_model_characterization')
            logging.warning(traceback.format_exc())
            logging.warning(sys.exc_info()[2])

        return new_model_ref

    def _check_model_obj_version(self, model_df, workspace_name):
        model_refs = model_df.index.tolist()

        for model_ref in model_refs:
            # should not cache the object since object might change after running FBA run_model_characterization
            latest_model_ref = '/'.join(model_ref.split('/')[:2])
            model_obj = self.dfu.get_objects({'object_refs': [latest_model_ref]})['data'][0]
            model_info = model_obj['info']
            model_data = model_obj['data']
            model_name = model_info[1]
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
                new_model_ref = self._run_model_characterization(workspace_name, model_name)
                if new_model_ref:
                    logging.info('created new FBAModel object {} from {}'.format(new_model_ref,
                                                                                 model_ref))

            idx = model_df.index.values.tolist().index(model_ref)
            model_df.index.values[idx] = '/'.join(model_ref.split('/')[:2])

    def _build_model_comparison_data(self, model_df):

        logging.info('Start building overall and reaction statistics')

        model_refs = model_df.index.tolist()
        overall_stats = list()
        reaction_stats = list()

        for model_ref in model_refs:
            model_obj = self._get_model_obj(model_ref)
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
        elif normalization_type == 'dividepathwaysize':
            pass
        else:
            logging.warning('Unexpected normalization type: {}'.format(normalization_type))

        return nor_df

    def _get_fc_profile_heatmap_data(self, field_type, model_refs):
        logging.info('Start building functional profile heatmap data for {}'.format(field_type))

        nor_type = field_type.split('_')[-1]
        if nor_type not in ['zscore', 'rownormalization', 'dividepathwaysize']:
            nor_type = None

        if nor_type:
            field_type = '_'.join(field_type.split('_')[:-1])

        # calculate available pathway_ids for all models
        pathway_ids = set()
        pathway_id_name_map = dict()
        pathway_id_class2_map = dict()
        for model_ref in model_refs:
            model_obj = self._get_model_obj(model_ref)
            model_data = model_obj['data']
            attributes = model_data.get('attributes', {})
            pathways = attributes.get('pathways', {})
            model_pathway_ids = pathways.keys()
            pathway_ids = pathway_ids | set(model_pathway_ids)

            for pathway_id in model_pathway_ids:
                if not pathway_id_name_map.get(pathway_id):
                    pathway_data = pathways.get(pathway_id, {})
                    pathway_name = pathway_data.get('name', 'NO PATHWAY NAME')
                    pathway_id_name_map[pathway_id] = pathway_name

                if not pathway_id_class2_map.get(pathway_id):
                    pathway_data = pathways.get(pathway_id, {})
                    pathway_classes = pathway_data.get('classes', ['', 'NO PATHWAY CLASS NAME'])
                    try:
                        pathway_class2 = pathway_classes[1]
                    except Exception:
                        pathway_class2 = 'NO PATHWAY CLASS NAME'
                    pathway_id_class2_map[pathway_id] = pathway_class2

        # append index number to each value in case of duplicates
        for index, key in enumerate(pathway_id_name_map):
            pathway_name = pathway_id_name_map.get(key) + ' [{}]'.format(index)
            pathway_class2 = pathway_id_class2_map.get(key) + ' [{}]'.format(index)
            pathway_id_name_map[key] = pathway_name
            pathway_id_class2_map[key] = pathway_class2

        pathway_df = pd.DataFrame(index=pathway_ids)

        # fetch pathway values
        for model_ref in model_refs:
            model_obj = self._get_model_obj(model_ref)
            model_data = model_obj['data']
            attributes = model_data.get('attributes', {})
            pathways = attributes.get('pathways', {})

            fetched_pathway_value = list()
            for pathway_id in pathway_df.index.to_list():
                pathway_data = pathways.get(pathway_id, {})
                pathway_data['total_functional_coverage'] = pathway_data.get(
                    'average_coverage_per_reaction', 0) * pathway_data.get('functional_rxn', 0)
                pathway_value = pathway_data.get(field_type, 0)
                if nor_type == 'dividepathwaysize':
                    pathway_size = pathway_data.get('pathway_size', 1)
                    if pathway_size <= 0:
                        pathway_size = 1
                    pathway_value = pathway_value / float(pathway_size)
                fetched_pathway_value.append(pathway_value)

            pathway_df[model_ref] = fetched_pathway_value

        if nor_type:
            pathway_df = self._normalize_data(pathway_df, nor_type)

        # run cluster on pathway dataframe
        col_ordered_label = self._compute_cluster_label_order(pathway_df.T.values.tolist(),
                                                              pathway_df.T.index.tolist())

        idx_ordered_label = self._compute_cluster_label_order(pathway_df.values.tolist(),
                                                              pathway_df.index.tolist())

        pathway_df = pathway_df.reindex(index=idx_ordered_label, columns=col_ordered_label)

        pathway_col_mapping_info = {'Pathway Name': [pathway_id_name_map[key] for key in pathway_df.index],
                                    'Pathway Class': [pathway_id_class2_map[key] for key in pathway_df.index]}

        return pathway_df, pathway_col_mapping_info

    def _get_pathway_heatmap_data(self, field_type, model_refs):

        logging.info('Start building pathway heatmap data for {}'.format(field_type))

        ori_field_type = copy.deepcopy(field_type)

        nor_type = field_type.split('_')[-1]
        if nor_type not in ['zscore', 'rownormalization', 'dividepathwaysize']:
            nor_type = None

        if nor_type:
            field_type = '_'.join(field_type.split('_')[:-1])

        first_model_ref = model_refs[0]
        model_obj = self._get_model_obj(first_model_ref)
        model_data = model_obj['data']
        model_info = model_obj['info']
        attributes = model_data.get('attributes', {})
        pathways = attributes.get('pathways', {})

        pathway_ids = list()
        pathway_names = list()
        pathway_class2_names = list()
        fetched_pathway_value = list()
        pathway_name_id_map = dict()
        pathway_count = 0
        for pathway_id, pathway_data in pathways.items():
            if ori_field_type == 'functional_rxn':
                pathway_ids.append(pathway_id)
            elif pathway_id in self.functional_rxn_pathways:
                pathway_ids.append(pathway_id)
            else:
                continue
            pathway_data['total_functional_coverage'] = pathway_data.get(
                'average_coverage_per_reaction', 0) * pathway_data.get('functional_rxn', 0)
            pathway_name = pathway_data.get('name') + ' [{}]'.format(pathway_count)
            pathway_name_id_map[pathway_name] = pathway_id
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
            model_obj = self._get_model_obj(model_ref)
            model_data = model_obj['data']
            model_info = model_obj['info']
            attributes = model_data.get('attributes', {})
            pathways = attributes.get('pathways', {})
            # TODO: should check/add newly found pathways found in models
            fetched_pathway_value = list()
            for pathway_id in pathway_ids:
                pathway_data = pathways.get(pathway_id, {})
                pathway_data['total_functional_coverage'] = pathway_data.get(
                    'average_coverage_per_reaction', 0) * pathway_data.get('functional_rxn', 0)
                pathway_value = pathway_data.get(field_type, 0)
                if nor_type == 'dividepathwaysize':
                    pathway_size = pathway_data.get('pathway_size', 1)
                    if pathway_size <= 0:
                        pathway_size = 1
                    pathway_value = pathway_value / float(pathway_size)
                fetched_pathway_value.append(pathway_value)

            model_name = model_info[1] + ' [{}]'.format(model_refs.index(model_ref))
            pathway_df[model_name] = fetched_pathway_value

        if ori_field_type == 'functional_rxn':
            pathway_df = pathway_df.loc[(pathway_df != 0).any(1)]
            self.functional_rxn_pathways = [pathway_name_id_map[x]
                                            for x in pathway_df.index.tolist()]

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

    def _build_model_set_meta(self, model_df):
        logging.info('Start building heatmap metadata')

        model_refs = model_df.index.tolist()
        meta_data_names = model_df.columns.tolist()

        heatmap_meta = dict()
        for meta_data_name in meta_data_names:
            model_meta_mapping = dict()
            model_meta = model_df.loc[:, meta_data_name].tolist()

            for i, model_ref in enumerate(model_refs):
                model_meta_data = model_meta[i] + ' [{}]'.format(i)
                model_meta_mapping.update({model_ref: model_meta_data})

            heatmap_meta.update({meta_data_name: model_meta_mapping})

        model_name_mapping = dict()
        for i, model_ref in enumerate(model_refs):
            model_obj = self._get_model_obj(model_ref)
            model_name = model_obj['info'][1] + ' [{}]'.format(i)
            model_name_mapping.update({model_ref: model_name})

        heatmap_meta.update({'Model Name': model_name_mapping})

        return heatmap_meta

    def _build_heatmap_meta(self, model_df):

        logging.info('Start building heatmap metadata')

        model_refs = model_df.index.tolist()

        meta_data_names = model_df.columns.tolist()

        model_names = list()
        for i, model_ref in enumerate(model_refs):
            model_obj = self._get_model_obj(model_ref)
            model_name = model_obj['info'][1]
            model_names.append(model_name + ' [{}]'.format(i))

        heatmap_meta = dict()
        for meta_data_name in meta_data_names:
            model_meta_mapping = dict()
            model_meta = model_df.loc[:, meta_data_name].tolist()

            model_meta_data = list()
            for i, meta_data in enumerate(model_meta):
                model_meta_data.append(meta_data + ' [{}]'.format(i))

            for i, model_name in enumerate(model_names):
                model_meta_mapping.update({model_name: model_meta_data[i]})

            heatmap_meta.update({meta_data_name: model_meta_mapping})

        return heatmap_meta

    @staticmethod
    def _convert_size(size_bytes):
        if size_bytes == 0:
            return "0B"
        size_name = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
        i = int(math.floor(math.log(size_bytes, 1024)))
        p = math.pow(1024, i)
        s = round(size_bytes / p, 2)
        return "%s %s" % (s, size_name[i])

    def _calculate_object_size(self, func_profile_data):
        json_size = 0
        try:
            logging.info('start calculating object size')
            json_object = json.dumps(func_profile_data).encode("utf-8")
            json_size = len(json_object)
            size_str = self._convert_size(json_size)
            logging.info('serialized object JSON size: {}'.format(size_str))
        except Exception:
            logging.info('failed to calculate object size')

        return json_size

    def _create_func_profile(self, pathway_df, func_profile_obj_name, workspace_name,
                             attri_mapping_ref):
        logging.info('Start creating functional profile object: {}'.format(func_profile_obj_name))

        if not isinstance(workspace_name, int):
            workspace_id = self.dfu.ws_name_to_id(workspace_name)
        else:
            workspace_id = workspace_name

        func_profile_data = dict()

        func_profile_data['original_matrix_ref'] = attri_mapping_ref
        func_profile_data['profile_type'] = 'ModelSet'
        func_profile_data['profile_category'] = 'community'

        profile_data = {'row_ids': pathway_df.index.tolist(),
                        'col_ids': pathway_df.columns.tolist(),
                        'values': pathway_df.values.tolist()}

        func_profile_data['data'] = profile_data

        obj_size = self._calculate_object_size(func_profile_data)

        MB_200 = 200 * 1024 * 1024
        GB_1 = 1 * 1024 * 1024 * 1024

        if obj_size > GB_1:
            raise ValueError('Object is too large')
        elif obj_size <= MB_200:
            logging.info('Starting saving object via DataFileUtil')
            info = self.dfu.save_objects({
                "id": workspace_id,
                "objects": [{
                    "type": 'KBaseProfile.FunctionalProfile',
                    "data": func_profile_data,
                    "name": func_profile_obj_name
                }]
            })[0]
        else:
            logging.info('Starting saving object via WsLargeDataIO')
            data_path = os.path.join(self.scratch,
                                     func_profile_obj_name + "_" + str(uuid.uuid4()) + ".json")
            logging.info('Dumpping object data to file: {}'.format(data_path))
            json.dump(func_profile_data, open(data_path, 'w'))

            info = self.ws_large_data.save_objects({
                "id": workspace_id,
                "objects": [{
                    "type": 'KBaseProfile.FunctionalProfile',
                    "data_json_file": data_path,
                    "name": func_profile_obj_name
                }]
            })[0]

        obj_ref = "%s/%s/%s" % (info[6], info[0], info[4])

        return obj_ref

    def _build_func_profile_data(self, model_df, profile_types, attri_mapping_name,
                                 workspace_name, attri_mapping_ref):

        logging.info('Start building functional profile data')
        model_refs = model_df.index.tolist()

        fc_profile_data = dict()

        fc_profile_refs = list()
        for profile_type in profile_types:
            for suffix in ['', '_zscore', '_rownormalization', '_dividepathwaysize']:
                calculation_type = profile_type + suffix

                pathway_df, pathway_col_mapping_info = self._get_fc_profile_heatmap_data(
                    calculation_type,
                    model_refs)

                fc_profile_data.update({calculation_type: {'values': pathway_df.values.tolist(),
                                                           'compound_names': pathway_df.columns.tolist(),
                                                           'index': pathway_df.index.tolist(),
                                                           'pathways': pathway_col_mapping_info}})

                if suffix == '':
                    fcp_ref = self._create_func_profile(pathway_df,
                                                        attri_mapping_name +
                                                        '_{}'.format(profile_type),
                                                        workspace_name,
                                                        attri_mapping_ref)
                    fc_profile_refs.append(fcp_ref)

        return fc_profile_data, fc_profile_refs

    def _build_heatmap_data(self, model_df):

        logging.info('Start building heatmap data')

        model_refs = model_df.index.tolist()

        heatmap_data = dict()

        pathway_types = ['functional_rxn', 'functional_rxn_zscore',
                         'functional_rxn_rownormalization', 'functional_rxn_dividepathwaysize',
                         'total_functional_coverage', 'total_functional_coverage_zscore',
                         'total_functional_coverage_rownormalization',
                         'total_functional_coverage_dividepathwaysize',
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
            'total_functional_coverage': 'Total Functional Coverage',
            'total_functional_coverage_zscore': 'Total Functional Coverage (Z-Score)',
            'total_functional_coverage_rownormalization': 'Total Functional Coverage (Divide Row Maximum Value)',
            'total_functional_coverage_dividepathwaysize': 'Total Functional Coverage (Divide Pathway Size)',
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

    def _create_fc_profile_heatmap_drop_down(self, output_directory, profile_datas, heatmap_meta):
        logging.info('Start building heatmap tab')

        suffix = str(uuid.uuid4())

        profile_name_profile_page_map = dict()
        profile_name_profile_page_json_name = 'profile_name_profile_page_map_{}.json'.format(
            suffix)
        profile_name_profile_page_json = os.path.join(
            output_directory, profile_name_profile_page_json_name)

        for profile_name, profile_data in profile_datas.items():
            model_set_refs = profile_data['compound_names']
            heatmap_meta_copy = copy.deepcopy(heatmap_meta)

            model_names = [heatmap_meta_copy['Model Name'][model_set_ref]
                           for model_set_ref in model_set_refs]
            df = pd.DataFrame(profile_data['values'],
                              index=profile_data['pathways']['Pathway Name'],
                              columns=model_names)

            heatmap_name = '{}_heatmap_dropdown.html'.format(profile_name.replace(' ', '_'))
            heatmap_path = os.path.join(output_directory, heatmap_name)

            profile_name_profile_page_map[profile_name] = heatmap_name
            # Initialize figure by creating upper dendrogram
            fig = ff.create_dendrogram(df.T.values,
                                       orientation='bottom',
                                       labels=df.columns,
                                       colorscale=['lightsteelblue'],
                                       distfun=lambda x: pdist(x),
                                       linkagefun=lambda x: linkage(x, 'ward'))
            for i in range(len(fig['data'])):
                fig['data'][i]['yaxis'] = 'y2'

            # Create Side Dendrogram
            dendro_side = ff.create_dendrogram(df.values,
                                               orientation='left',
                                               labels=df.index,
                                               colorscale=['lightsteelblue'] * 8,
                                               distfun=lambda x: pdist(x),
                                               linkagefun=lambda x: linkage(x, 'ward'))
            for i in range(len(dendro_side['data'])):
                dendro_side['data'][i]['xaxis'] = 'x2'

            # Add Side Dendrogram Data to Figure
            for data in dendro_side['data']:
                fig.add_trace(data)

            # Create Heatmap
            idx_ordered_label = dendro_side['layout']['yaxis']['ticktext']
            col_ordered_label = fig['layout']['xaxis']['ticktext']

            colorscale = [
                [0, 'rgb(47, 15, 61)'],
                [0.2, 'rgb(107, 24, 93)'],
                [0.4, 'rgb(168, 40, 96)'],
                [0.6, 'rgb(215, 48, 31)'],
                [0.8, 'rgb(252,141,89)'],
                [1.0, 'rgb(255,247,236)']]

            df = df.reindex(index=idx_ordered_label, columns=col_ordered_label)

            heatmap = go.Heatmap(
                z=df.values,
                x=df.columns,
                y=df.index,
                hoverongaps=False,
                colorscale=colorscale)

            heatmap['x'] = fig['layout']['xaxis']['tickvals']
            heatmap['y'] = dendro_side['layout']['yaxis']['tickvals']

            # Add Heatmap Data to Figure
            fig.add_trace(heatmap)

            width = max(100 * df.columns.size, 1400)
            height = max(10 * df.index.size, 1000)
            y2_height = 100
            x2_width = 150
            y2_offset = y2_height / height
            x2_offset = x2_width / width

            fig.update_layout({'plot_bgcolor': 'rgba(0,0,0,0)',
                               'showlegend': False,
                               'width': width,
                               'height': height,
                               'autosize': True,
                               'hovermode': 'closest'})

            # Edit xaxis
            fig.update_layout(xaxis={'domain': [0, max(1-x2_offset, 0.825)],
                                     'mirror': False,
                                     'showgrid': False,
                                     'showline': False,
                                     'zeroline': False,
                                     'automargin': True,
                                     'tickangle': 45,
                                     'tickfont': dict(color='black', size=8),
                                     'ticks': ""})
            # Edit xaxis2
            fig.update_layout(xaxis2={'domain': [max(1-x2_offset, 0.825), 1],
                                      'mirror': False,
                                      'showgrid': False,
                                      'showline': False,
                                      'zeroline': False,
                                      'showticklabels': False,
                                      'ticks': ""})

            # Edit yaxis
            fig.update_layout(yaxis={'domain': [0, max(1-y2_offset, 0.85)],
                                     'mirror': False,
                                     'showgrid': False,
                                     'showline': False,
                                     'zeroline': False,
                                     'automargin': True,
                                     'tickfont': dict(color='black', size=8),
                                     'ticks': "",
                                     'ticktext': dendro_side['layout']['yaxis']['ticktext'],
                                     'tickvals': dendro_side['layout']['yaxis']['tickvals']
                                     })
            # Edit yaxis2
            fig.update_layout(yaxis2={'domain': [max(1-y2_offset, 0.85), 1],
                                      'mirror': False,
                                      'showgrid': False,
                                      'showline': False,
                                      'zeroline': False,
                                      'showticklabels': False,
                                      'ticks': ""})

            yaxis_copy = copy.deepcopy(fig['layout']['yaxis'])
            xaxis_copy = copy.deepcopy(fig['layout']['xaxis'])

            yaxis_copy['ticktext'] = profile_data['pathways']['Pathway Name']
            yaxis_pathway_name = copy.deepcopy(yaxis_copy)

            yaxis_copy['ticktext'] = profile_data['pathways']['Pathway Class']
            yaxis_pathway_class = copy.deepcopy(yaxis_copy)

            model_names = xaxis_copy['ticktext']

            update_x_buttons = list()
            model_ref_name_map = heatmap_meta_copy['Model Name']
            model_refs = list()

            for model_name in model_names:
                model_ref = list(model_ref_name_map.keys())[list(
                    model_ref_name_map.values()).index(model_name)]
                model_refs.append(model_ref)

            xaxis_model_name = copy.deepcopy(xaxis_copy)
            update_x_buttons.append(dict(label='Model Name',
                                         method="relayout",
                                         args=["xaxis", xaxis_model_name]))
            heatmap_meta_copy.pop('Model Name')
            for meta_name, model_ref_meta_map in heatmap_meta_copy.items():
                new_xaxis_names = [model_ref_meta_map[model_ref] for model_ref in model_refs]
                xaxis_copy['ticktext'] = new_xaxis_names
                xaxis_updated = copy.deepcopy(xaxis_copy)
                update_x_buttons.append(dict(label=meta_name,
                                             method="relayout",
                                             args=["xaxis", xaxis_updated]))

            # Add dropdowns
            button_layer_1_height = min(1 + y2_height/height, 1.06)
            fig.update_layout(
                updatemenus=[
                    dict(
                        buttons=list([
                            dict(label="Payway Name",
                                 method="relayout",
                                 args=["yaxis", yaxis_pathway_name]),
                            dict(label="Payway Class",
                                 method="relayout",
                                 args=["yaxis", yaxis_pathway_class]),
                        ]),
                        direction="down",
                        pad={"r": 10, "t": 10},
                        showactive=True,
                        x=0.06,
                        xanchor="left",
                        y=button_layer_1_height,
                        yanchor="top"
                    ),
                    dict(
                        buttons=update_x_buttons,
                        direction="down",
                        pad={"r": 10, "t": 10},
                        showactive=True,
                        x=0.27,
                        xanchor="left",
                        y=button_layer_1_height,
                        yanchor="top"
                    ),

                ]
            )

            fig.update_layout(
                annotations=[
                    dict(text="Pathway<br>Label", x=0, xref="paper",
                         y=button_layer_1_height - min(5/height, 0.01),
                         yref="paper", showarrow=False),
                    dict(text="Metadata<br>Label", x=0.2, xref="paper",
                         y=button_layer_1_height - min(5/height, 0.01),
                         yref="paper", align="left", showarrow=False),
                ])

            plot(fig, filename=heatmap_path)

        with open(profile_name_profile_page_json, 'w') as outfile:
            json.dump(profile_name_profile_page_map, outfile)

        data_info = ""
        for data_name in profile_datas.keys():
            data_info += """\n<option value="{}">{}</option>\n""".format(data_name, data_name)

        heatmap_html_name = 'heatmap_{}.html'.format(suffix)
        heatmap_html = os.path.join(output_directory, heatmap_html_name)

        with open(heatmap_html, 'w') as heatmap_html:
            with open(os.path.join(os.path.dirname(__file__),
                                   'templates', 'functional_profiles_dropdown_template.html'),
                      'r') as heatmap_template_file:
                heatmap_template = heatmap_template_file.read()
                heatmap_template = heatmap_template.replace('<!-- data_info -->',
                                                            data_info)
                heatmap_template = heatmap_template.replace('profile_name_profile_page_map.json',
                                                            profile_name_profile_page_json_name)
                heatmap_template = heatmap_template.replace(
                    'heatmap_dropdown.html',
                    profile_name_profile_page_map[list(profile_datas.keys())[0]])
                heatmap_html.write(heatmap_template)

        return heatmap_html_name

    def _create_fc_profile_heatmap(self, output_directory, profile_datas, heatmap_meta):

        logging.info('Start building heatmap tab')

        suffix = str(uuid.uuid4())

        heatmap_data_json_name = 'heatmap_data_{}.json'.format(suffix)
        heatmap_meta_json_name = 'heatmap_meta_{}.json'.format(suffix)
        heatmap_data_json = os.path.join(output_directory, heatmap_data_json_name)
        heatmap_meta_json = os.path.join(output_directory, heatmap_meta_json_name)
        with open(heatmap_data_json, 'w') as outfile:
            json.dump(profile_datas, outfile)
        with open(heatmap_meta_json, 'w') as outfile:
            json.dump(heatmap_meta, outfile)

        metadata_info = ""
        metadata_info += """\n<option value="{}">{}</option>\n""".format(
            'Model Name', 'Model Name')
        metadata_names = list(heatmap_meta.keys())
        metadata_names.remove('Model Name')
        for meta_name in metadata_names:
            metadata_info += """\n<option value="{}">{}</option>\n""".format(meta_name, meta_name)

        data_info = ""
        for data_name in profile_datas.keys():
            data_info += """\n<option value="{}">{}</option>\n""".format(data_name, data_name)

        pathway_info = ""
        for pathway_name in profile_datas.get(list(profile_datas.keys())[0])['pathways'].keys():
            pathway_info += """\n<option value="{}">{}</option>\n""".format(
                pathway_name, pathway_name)

        heatmap_html_name = 'heatmap_{}.html'.format(suffix)
        heatmap_html = os.path.join(output_directory, heatmap_html_name)

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
                heatmap_template = heatmap_template.replace('heatmap_data.json',
                                                            heatmap_data_json_name)
                heatmap_template = heatmap_template.replace('heatmap_meta.json',
                                                            heatmap_meta_json_name)
                heatmap_html.write(heatmap_template)

        return heatmap_html_name

    def _generate_fc_profile_report(self, fc_profile_data, heatmap_meta, profile_types):
        logging.info('Start building html report')

        pathway_name_map = {
            'functional_rxn': 'Functional Reaction',
            'functional_rxn_zscore': 'Functional Reaction (Z-Score)',
            'functional_rxn_rownormalization': 'Functional Reaction (Divide Row Maximum Value)',
            'functional_rxn_dividepathwaysize': 'Functional Reaction (Divide Pathway Size)',
            'total_functional_coverage': 'Total Functional Coverage',
            'total_functional_coverage_zscore': 'Total Functional Coverage (Z-Score)',
            'total_functional_coverage_rownormalization': 'Total Functional Coverage (Divide Row Maximum Value)',
            'total_functional_coverage_dividepathwaysize': 'Total Functional Coverage (Divide Pathway Size)',
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
            'gene_count_dividepathwaysize': 'Gene Count (Divide Pathway Size)'}

        html_report = list()

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        logging.info('Start building functional profiles report in dir: {}'.format(
            output_directory))
        self._mkdir_p(output_directory)
        model_set_file_path = os.path.join(output_directory, 'model_set_functional_profiles.html')

        tab_def_content = ''
        tab_content = ''

        # build first profile heatmap
        first_profile_type = profile_types[0]
        profile_datas = dict()
        for suffix in ['', '_zscore', '_rownormalization', '_dividepathwaysize']:
            profile_name = first_profile_type + suffix
            profile_datas.update({pathway_name_map[profile_name]: fc_profile_data[profile_name]})

            if suffix == '':
                html_tab_name = pathway_name_map[profile_name]

        # heatmap_page = self._create_fc_profile_heatmap(output_directory, profile_datas,
        #                                                heatmap_meta)
        heatmap_page = self._create_fc_profile_heatmap_drop_down(output_directory, profile_datas,
                                                                 heatmap_meta)

        tab_id = html_tab_name.replace(" ", "")
        tab_def_content += """
        <div class="tab">
            <button class="tablinks" onclick="openTab(event, '{}')" id="defaultOpen">{}</button>
        """.format(tab_id, html_tab_name)

        page_content = ''
        page_content += '\n<iframe height="1300px" width="100%" '
        page_content += 'src="{}" style="border:none;"></iframe>'.format(heatmap_page)

        tab_content += """\n<div id="{}" class="tabcontent">{}</div>""".format(tab_id,
                                                                               page_content)

        # build reset profiles heatmap
        for profile_type in profile_types[1:]:
            profile_datas = dict()
            for suffix in ['', '_zscore', '_rownormalization', '_dividepathwaysize']:
                profile_name = profile_type + suffix
                profile_datas.update(
                    {pathway_name_map[profile_name]: fc_profile_data[profile_name]})

                if suffix == '':
                    html_tab_name = pathway_name_map[profile_name]

            # heatmap_page = self._create_fc_profile_heatmap(output_directory, profile_datas,
            #                                                heatmap_meta)
            heatmap_page = self._create_fc_profile_heatmap_drop_down(output_directory,
                                                                     profile_datas,
                                                                     heatmap_meta)

            tab_id = html_tab_name.replace(" ", "")
            tab_def_content += """
                <button class="tablinks" onclick="openTab(event, '{}')">{}</button>
            """.format(tab_id, html_tab_name)

            page_content = ''
            page_content += '\n<iframe height="1300px" width="100%" '
            page_content += 'src="{}" style="border:none;"></iframe>'.format(heatmap_page)

            tab_content += """\n<div id="{}" class="tabcontent">{}</div>""".format(tab_id,
                                                                                   page_content)

        tab_def_content += """</div>"""
        visualization_content = tab_def_content + tab_content
        with open(model_set_file_path, 'w') as result_file:
            with open(os.path.join(os.path.dirname(__file__), 'templates', 'model_set_template.html'),
                      'r') as report_template_file:
                report_template = report_template_file.read()
                report_template = report_template.replace('<p>Visualization_Content</p>',
                                                          visualization_content)
                result_file.write(report_template)

        report_shock_id = self.dfu.file_to_shock({'file_path': output_directory,
                                                  'pack': 'zip'})['shock_id']

        html_report.append({'shock_id': report_shock_id,
                            'name': os.path.basename(model_set_file_path),
                            'label': os.path.basename(model_set_file_path),
                            'description': 'HTML summary report for Model Set Functional Profile'
                            })

        return html_report

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
            pathway_info += """\n<option value="{}">{}</option>\n""".format(
                pathway_name, pathway_name)

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
        self.ws_large_data = WsLargeDataIO(self.callback_url, service_ver='beta')
        self.fba_tools = fba_tools(self.callback_url)
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        self.obj_cache = dict()

        self.functional_rxn_pathways = list()

    def model_set_to_functional_profiles(self, params):
        workspace_name = params.get('workspace_name')
        attri_mapping_ref = params.get('attri_mapping_ref')
        profile_types = params.get('profile_types', {})

        profile_types = [item[0] for item in profile_types.items() if item[1]]
        if not profile_types:
            raise ValueError('Please choose at least one profile type to be created')

        attri_mapping_obj = self.dfu.get_objects(
            {'object_refs': [attri_mapping_ref]})['data'][0]
        attri_mapping_data = attri_mapping_obj['data']
        attri_mapping_name = attri_mapping_obj['info'][1]
        attributes = pd.DataFrame(attri_mapping_data['attributes'])
        instances = pd.DataFrame(attri_mapping_data['instances'])
        model_df = attributes.join(instances)
        model_name_idx = model_df['attribute'].tolist().index('model_name')
        model_df.drop(columns=['source'], index=[model_name_idx], inplace=True)
        model_df = model_df.T
        model_df.columns = model_df.loc['attribute', :]
        model_df.drop(index=['attribute'], inplace=True)

        self._check_model_obj_version(model_df, workspace_name)

        fc_profile_data, fc_profile_refs = self._build_func_profile_data(model_df, profile_types,
                                                                         attri_mapping_name,
                                                                         workspace_name,
                                                                         attri_mapping_ref)
        heatmap_meta = self._build_model_set_meta(model_df)

        html_files = self._generate_fc_profile_report(fc_profile_data, heatmap_meta, profile_types)

        objects_created = [{'ref': ref,
                            'description': 'Functional Porfile'} for ref in fc_profile_refs]
        report_params = {'message': '',
                         'workspace_name': workspace_name,
                         'objects_created': objects_created,
                         'html_links': html_files,
                         'direct_html_link_index': 0,
                         'html_window_height': 1350,
                         'report_object_name': 'model_set_to_func_profile_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url, token=self.token)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref'],
                         'functional_profile_refs': fc_profile_refs}

        return report_output

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
                         'html_window_height': 1250,
                         'report_object_name': 'model_comparison_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url, token=self.token)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output
