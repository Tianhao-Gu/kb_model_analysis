# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os

from kb_model_analysis.Utils.HeatmapUtil import HeatmapUtil
from kb_model_analysis.Utils.TemplateUtil import TemplateUtil
#END_HEADER


class kb_model_analysis:
    '''
    Module Name:
    kb_model_analysis

    Module Description:
    A KBase module: kb_model_analysis
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/Tianhao-Gu/kb_model_analysis.git"
    GIT_COMMIT_HASH = "755733400ae81943052ee6e6be3a9e1c91414257"

    #BEGIN_CLASS_HEADER
    @staticmethod
    def validate_params(params, expected, opt_param=set()):
        """Validates that required parameters are present. Warns if unexpected parameters appear"""
        expected = set(expected)
        opt_param = set(opt_param)
        pkeys = set(params)
        if expected - pkeys:
            raise ValueError("Required keys {} not in supplied parameters"
                             .format(", ".join(expected - pkeys)))
        defined_param = expected | opt_param
        for param in params:
            if param not in defined_param:
                logging.warning("Unexpected parameter {} supplied".format(param))
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.config['SDK_CALLBACK_URL'] = os.environ['SDK_CALLBACK_URL']
        self.config['KB_AUTH_TOKEN'] = os.environ['KB_AUTH_TOKEN']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        self.heatmap_util = HeatmapUtil(config)
        self.template_util = TemplateUtil(config)
        #END_CONSTRUCTOR
        pass


    def model_heatmap_analysis(self, ctx, params):
        """
        :param params: instance of type "HeatmapAnalysisParams" -> structure:
           parameter "workspace_name" of String, parameter
           "staging_file_path" of String
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN model_heatmap_analysis
        self.validate_params(params, ['workspace_name', 'staging_file_path'])
        output = self.heatmap_util.run_model_heatmap_analysis(params)
        #END model_heatmap_analysis

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method model_heatmap_analysis return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def create_heatmap_analysis_template(self, ctx, params):
        """
        :param params: instance of type "HeatmapAnalysisTemplateParams" ->
           structure: parameter "workspace_id" of Long, parameter
           "object_types" of list of String, parameter "workspace_scope" of
           String, parameter "metadata_fields" of String
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN create_heatmap_analysis_template
        self.validate_params(params, ['workspace_id'],
                             opt_param=['object_types', 'metadata_fields', 'workspace_scope'])
        output = self.template_util.create_heatmap_analysis_template(params)
        #END create_heatmap_analysis_template

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method create_heatmap_analysis_template return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
