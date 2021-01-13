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
from installed_clients.FakeObjectsForTestsClient import FakeObjectsForTests


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

        if model_ref.count('/') == 2 and model_ref.endswith('1'):

            # return fake ModelSet object
            attri_mapping_data = {
                "attributes": [
                    {
                        "attribute": "model_name",
                        "source": "upload"
                    },
                    {
                        "attribute": "metadata_3",
                        "source": "upload"
                    },
                    {
                        "attribute": "metadata_2",
                        "source": "upload"
                    },
                    {
                        "attribute": "metadata_1",
                        "source": "upload"
                    }
                ],
                "instances": {
                    '1/1/2': [
                        "Escherichia_coli.mdl.gapfilled-1",
                        "foo_1",
                        "foo_2",
                        "foo_3"
                    ],
                    '1/2/2': [
                        "Escherichia_coli.mdl.gapfilled-0",
                        "foo_1",
                        "bar_2",
                        ""
                    ],
                    '1/3/2': [
                        "Escherichia_coli.mdl.gapfilled",
                        "bar_1",
                        "",
                        "bar_3"
                    ]
                },
                "ontology_mapping_method": "User Curation - ISA format"
            }

            data = {'data': [{'data': attri_mapping_data,
                              'info': [1, 'FakeFBAModelSet', 'KBaseExperiments.FBAModelSet‑1.0']}]}
            return data

        # build fake FBAModel object
        random.seed(None)

        fake_model_data = {}

        fake_model_data['modelcompounds'] = ['modelcompounds'] * random.randint(0, 32)
        fake_model_data['modelreactions'] = ['modelreactions'] * random.randint(0, 32)

        attributes = {}

        attributes['base_gapfilling'] = random.randint(0, 32)
        attributes['core_gapfilling'] = random.randint(0, 32)
        attributes['gene_count'] = random.randint(0, 32)
        attributes['base_atp'] = random.randint(0, 32)

        fbas = dict()
        for fba_name in ['auxomedia', 'complete']:
            fbas[fba_name] = {'biomass': random.random(),
                              'Blocked': random.randint(0, 32),
                              'Negative': random.randint(0, 32),
                              'Positive': random.randint(0, 32),
                              'PositiveVariable': random.randint(0, 32),
                              'NegativeVariable': random.randint(0, 32),
                              'Variable': random.randint(0, 32)}

        attributes['fbas'] = fbas

        auxotrophy_size = random.randint(0, 9)
        auxotrophy = dict()
        for i in range(auxotrophy_size):
            auxotrophy['compound_{}'.format(i)] = {'is_auxotrophic': (random.randint(0, 32) < 16)}
        attributes['auxotrophy'] = auxotrophy

        pathway_size = 20
        pathways = dict()
        for i in range(pathway_size):
            pathway = dict()
            pathway['name'] = 'pathway_name_{}'.format(i)

            class_1 = ['Metabolism', 'class_1_B', 'class_1_C'][random.randint(0, 32) % 3]
            class_2 = ["Carbohydrate metabolism",
                       "Lipid metabolism",
                       "Nucleotide metabolism",
                       "Amino acid metabolism"][random.randint(0, 32) % 4]
            pathway['classes'] = [class_1, class_2]
            # pathway['class_1'] = ['class_1_A', 'class_1_B', 'class_1_C'][random.randint(0, 32) % 3]
            # pathway['class_2'] = ['class_2_A', 'class_2_B', 'class_2_C'][random.randint(0, 32) % 3]
            pathway['gapfilled_rxn'] = random.randint(0, 32)
            pathway['functional_rxn'] = random.randint(0, 32)
            pathway['nonfunctional_rxn'] = random.randint(0, 32)
            pathway['pathway_size'] = random.randint(0, 32)
            pathway['is_present'] = (random.randint(0, 32) < 16)
            pathway['gene_count'] = random.randint(0, 32)
            pathway['average_genes_per_reaction'] = random.random()
            pathway['stddev_genes_per_reaction'] = random.random()
            pathway['average_coverage_per_reaction'] = random.random()
            pathway['stddev_coverage_per_reaction'] = random.random()
            pathways['pathway_{}'.format(i)] = pathway
        attributes['pathways'] = pathways

        fake_model_data['attributes'] = attributes

        random.seed(model_ref.split('/')[-1])
        model_name = 'model_' + str(random.randint(0, 1024))
        data = {'data': [{'data': fake_model_data,
                          'info': [1, model_name, 'KBaseFBA.FBAModel-14.0']}]}

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

        foft = FakeObjectsForTests(self.callback_url)
        obj_names = ['test_reads.1', 'test_reads.2']
        foft.create_fake_reads({'ws_name': self.wsName, 'obj_names': obj_names})
        params = {'workspace_id': self.wsId,
                  'object_type': ['KBaseFile.SingleEndLibrary', 'KBaseFile.PairedEndLibrary'],
                  'workspace_scope': 'all_workspace',
                  }
        returnVal = self.serviceImpl.create_heatmap_analysis_template(self.ctx, params)[0]

        self.assertIn('report_ref', returnVal)
        self.assertIn('report_name', returnVal)

    @patch.object(DataFileUtil, "get_objects", side_effect=mock_get_objects)
    def test_model_set_to_functional_profiles_app(self, get_objects):
        profile_types = {"functional_rxn": 1,
                         "nonfunctional_rxn": 0,
                         "gapfilled_rxn": 1,
                         "total_functional_coverage": 0,
                         "gene_count": 1}

        foft = FakeObjectsForTests(self.callback_url)
        obj_name = 'test_obj.1'
        info = foft.create_any_objects({'ws_name': self.wsName, 'obj_names': [obj_name]})[0]
        attri_mapping_ref = "%s/%s/%s" % (info[6], info[0], info[4])

        params = {'workspace_name': self.wsName,
                  'profile_types': profile_types,
                  'attri_mapping_ref': attri_mapping_ref}
        returnVal = self.serviceImpl.model_set_to_functional_profiles(self.ctx, params)[0]

        self.assertIn('report_ref', returnVal)
        self.assertIn('report_name', returnVal)
        self.assertIn('functional_profile_refs', returnVal)
