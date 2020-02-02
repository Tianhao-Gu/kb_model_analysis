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
    } HeatmapAnalysisParams;

    funcdef model_heatmap_analysis(HeatmapAnalysisParams params) returns (ReportResults output) authentication required;

    typedef structure {
        int workspace_id;
        list<string> object_types;
        string workspace_scope;
        string metadata_fields;
    } HeatmapAnalysisTemplateParams;

    funcdef create_heatmap_analysis_template(HeatmapAnalysisTemplateParams params) returns (ReportResults output) authentication required;

};
