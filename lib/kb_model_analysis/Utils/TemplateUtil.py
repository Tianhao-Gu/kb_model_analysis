import errno
import os
import uuid
import logging
import pandas as pd
from xlrd.biffh import XLRDError
import json
import shutil

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WorkspaceClient import Workspace
from kb_model_analysis.authclient import KBaseAuth


class TemplateUtil:

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

    def _generate_report(self, obj_df, workspace_id, workspace_num):

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        logging.info('Start building report files in dir: {}'.format(output_directory))
        self._mkdir_p(output_directory)

        template_file = os.path.join(output_directory, 'model_comparison_template.xlsx')

        obj_df.to_excel(template_file)

        file_links = [{'path': template_file,
                       'name': os.path.basename(template_file),
                       'label': 'model comparison template file',
                       'description': 'use this file for Model Comparison app'}]

        report_params = {'message': 'Found {} objects in {} narratives'.format(
                                                            len(obj_df.index), workspace_num),
                         'workspace_id': workspace_id,
                         'file_links': file_links,
                         'report_object_name': 'model_comparison_template_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url, token=self.token)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

    def _build_ws_lookup_table(self, wsinfos):
        ''' builds a lookup table, skips anything without a 'narrative' metadata field set '''
        ws_lookup_table = {}
        for ws_info in wsinfos:
            if 'narrative' in ws_info[8]:
                if ws_info[8]['narrative'].isdigit() and int(ws_info[8]['narrative']) > 0:
                    ws_lookup_table[ws_info[0]] = ws_info
        return ws_lookup_table

    def _build_workspace_ids(self, wsinfos):
        workspace_ids = list()
        for ws_info in wsinfos:
            if 'narrative' in ws_info[8]:
                if ws_info[8]['narrative'].isdigit() and int(ws_info[8]['narrative']) > 0:
                    workspace_ids.append(ws_info[0])
        return workspace_ids

    def _get_workspace_ids(self, all_workspace, private_workspace, user_id):
        workspace_ids = list()

        if all_workspace:
            wsinfos = self.ws_client.list_workspace_info({'excludeGlobal': 1,
                                                          'showDeleted': 0})
            workspace_ids = self._build_workspace_ids(wsinfos)
        elif private_workspace:
            wsinfos = self.ws_client.list_workspace_info({'owners': [user_id],
                                                          'showDeleted': 0})
            workspace_ids = self._build_workspace_ids(wsinfos)

        return workspace_ids

    def __init__(self, config):
        self.callback_url = config['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.token = config['KB_AUTH_TOKEN']
        self.dfu = DataFileUtil(self.callback_url)
        self.ws_client = Workspace(config['workspace-url'])
        self.auth_client = KBaseAuth(config['auth-service-url'])
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)

    def create_heatmap_analysis_template(self, params):
        workspace_id = params.get('workspace_id')

        metadata_fields = params.get('metadata_fields')
        if metadata_fields is None:
            metadata_fields = ''
        metadata_fields = list(set([i.strip() for i in metadata_fields.split(',')]))
        try:
            metadata_fields.remove('')
        except Exception:
            pass
        object_types = params.get('object_type', ['KBaseFBA.FBAModel'])

        user_id = self.auth_client.get_user(self.token)

        workspace_scope = params.get('workspace_scope', {})

        all_workspace = workspace_scope.get('all_workspace', 0)
        private_workspace = workspace_scope.get('private_workspace', 0)

        workspace_ids = self._get_workspace_ids(all_workspace, private_workspace, user_id)

        if not workspace_ids:
            workspace_ids = [workspace_id]

        obj_infos = list()
        for object_type in object_types:
            objs = self.ws_client.list_objects({'ids': workspace_ids,
                                                'type': object_type,
                                                'showDeleted': 0,
                                                'showHidden': 0})
            obj_infos.extend(objs)

        if not obj_infos:
            raise ValueError("Cannot find any matching objects in any narrative")

        obj_mapping = dict()
        for obj_info in obj_infos:
            obj_ref = "%s/%s/%s" % (obj_info[6], obj_info[0], obj_info[4])
            obj_mapping.update({obj_ref: obj_info[1]})

        obj_df = pd.DataFrame.from_dict(obj_mapping, orient='index',
                                        columns=['model_name'])

        total_records = len(obj_mapping)
        for metadata_field in metadata_fields:
            obj_df[metadata_field] = [None] * total_records
        obj_df.index.name = 'model_ref'

        report_output = self._generate_report(obj_df, workspace_id, len(workspace_ids))

        return report_output
