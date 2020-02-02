# -*- coding: utf-8 -*-
import os
import time
import unittest
from configparser import ConfigParser
from mock import patch
import json
import random

from kb_model_analysis.kb_model_analysisImpl import kb_model_analysis
from kb_model_analysis.kb_model_analysisServer import MethodContext
from kb_model_analysis.authclient import KBaseAuth as _KBaseAuth

from installed_clients.WorkspaceClient import Workspace
from installed_clients.DataFileUtilClient import DataFileUtil


class kb_model_analysisTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_model_analysis'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'kb_model_analysis',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = kb_model_analysis(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_ContigFilter_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa
        cls.wsId = ret[0]

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def mock_download_staging_file(params):

        file_path = params.get('staging_file_subdir_path')

        return {'copy_file_path': file_path}

    def mock_get_objects(params):

        model_ref = params['object_refs'][0]

        random.seed(model_ref.split('/')[-1])

        data = {'data': []}

        example_model_file = os.path.join('data', 'example_model.json')

        with open(example_model_file) as json_file:
            fake_model_data = json.load(json_file)

        model_name = 'model_' + str(random.randint(0, 1024))

        attributes = fake_model_data['attributes']
        pathway_keys = [s for s in attributes.keys() if ('pathways_' in s)]

        for pathway_key in pathway_keys:
            pathway_data = attributes[pathway_key]
            pathway_data['gf'] = random.randint(0, 64)
            pathway_data['nonblocked'] = random.randint(0, 128)
            pathway_data['rxn'] = random.randint(0, 64)

        for i in range(500):
            pathway_key = 'pathways_rn0{}'.format(str(1055 + i*2))
            attributes[pathway_key] = {'gf': random.randint(0, 64),
                                       'nonblocked': random.randint(0, 128),
                                       'rxn': random.randint(0, 64)}

        data = {'data': [{'data': fake_model_data, 'info': [1, model_name]}]}

        return data

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    def test_bad_params(self):
        with self.assertRaises(ValueError) as context:
            self.serviceImpl.model_heatmap_analysis(self.ctx, {'workspace_name': self.wsName})
            self.assertIn("Required keys", str(context.exception.args))

        with self.assertRaises(ValueError) as context:
            self.serviceImpl.create_heatmap_analysis_template(self.ctx, {})
            self.assertIn("Required keys", str(context.exception.args))

    @patch.object(DataFileUtil, "download_staging_file", side_effect=mock_download_staging_file)
    @patch.object(DataFileUtil, "get_objects", side_effect=mock_get_objects)
    def test_model_heatmap_analysis_app(self, download_staging_file, get_objects):
        params = {'workspace_name': self.wsName,
                  'staging_file_path': os.path.join('data', 'model_compare_temp.xlsx')}
        returnVal = self.serviceImpl.model_heatmap_analysis(self.ctx, params)[0]

        self.assertIn('report_ref', returnVal)
        self.assertIn('report_name', returnVal)

    def test_create_heatmap_analysis_template(self):
        params = {'workspace_id': self.wsId,
                  'workspace_scope': 'all_workspace'}
        returnVal = self.serviceImpl.create_heatmap_analysis_template(self.ctx, params)[0]

        self.assertIn('report_ref', returnVal)
        self.assertIn('report_name', returnVal)
