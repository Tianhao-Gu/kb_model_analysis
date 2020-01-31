/*
A KBase module: kb_model_analysis
*/

module kb_model_analysis {
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    typedef structure {
        string workspace_name;
        string staging_file_path;
    } HeatmapAnalysisParams;

    funcdef model_heatmap_analysis(HeatmapAnalysisParams params) returns (ReportResults output) authentication required;

};
