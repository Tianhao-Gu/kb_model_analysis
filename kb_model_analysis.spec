/*
A KBase module: kb_model_analysis
*/

module kb_model_analysis {
    /* A boolean - 0 for false, 1 for true.
    */
    typedef int boolean;

    /* An X/Y/Z style reference
    */
    typedef string obj_ref;

    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    typedef structure {
        string workspace_name;
        string staging_file_path;
        string attri_mapping_ref;
    } HeatmapAnalysisParams;

    funcdef model_heatmap_analysis(HeatmapAnalysisParams params) returns (ReportResults output) authentication required;

    /* workspace_scope: one of ['all_workspace', 'private_workspace', 'current_workspace'], default 'current_workspace'
       object_types: default ['KBaseFBA.FBAModel']
    */
    typedef structure {
        int workspace_id;
        list<string> object_types;
        string workspace_scope;
        string metadata_fields;
    } HeatmapAnalysisTemplateParams;

    funcdef create_heatmap_analysis_template(HeatmapAnalysisTemplateParams params) returns (ReportResults output) authentication required;

    /*
        required params:
        staging_file_subdir_path: subdirectory file path
        e.g.
          for file: /data/bulk/user_name/file_name
          staging_file_subdir_path is file_name
          for file: /data/bulk/user_name/subdir_1/subdir_2/file_name
          staging_file_subdir_path is subdir_1/subdir_2/file_name
        attribute_mapping_name: output ConditionSet object name
        workspace_id: workspace name/ID of the object
    */
    typedef structure {
        string staging_file_subdir_path;
        int workspace_id;
        string attribute_mapping_name;
    } FileToConditionSetParams;

    funcdef import_fbamodel_set_from_staging(FileToConditionSetParams params) returns (ReportResults output) authentication required;

};
