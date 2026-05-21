================================================
  NLOC    CCN   token  PARAM  length  location  
------------------------------------------------
     238     36   2664      7     275 shiny_server@1-275@schardExplorer/R/shiny_server.R
       2      1     25      0       2 filename@275-276@schardExplorer/R/shiny_server.R
      12      1    130      1      15 content@276-290@schardExplorer/R/shiny_server.R
       2      1     25      0       2 filename@293-294@schardExplorer/R/shiny_server.R
      10      3    135      1      10 content@294-303@schardExplorer/R/shiny_server.R
       4      1     16      1       7 qc_report@16-22@schardExplorer/R/report.R
      51      8    691      3      64 qc_report.default@22-85@schardExplorer/R/report.R
       4      1     57      3       7 qc_report.character@85-91@schardExplorer/R/report.R
       4      1     57      3       7 qc_report.list@91-97@schardExplorer/R/report.R
       4      1     57      3       7 qc_report.SingleCellExperiment@97-103@schardExplorer/R/report.R
       4      1     57      3      12 qc_report.Seurat@103-114@schardExplorer/R/report.R
      96     12    574      3     102 write_qc_html@114-215@schardExplorer/R/report.R
       7      1     90      4       9 run_dashboard_app@10-18@schardExplorer/R/shiny_app.R
       8      1     55      3       8 server@18-25@schardExplorer/R/shiny_app.R
      44      2    328      3      48 shiny_ui@1-48@schardExplorer/R/shiny_ui.R
       4      1     16      1       7 launch_qc_dashboard@15-21@schardExplorer/R/launch.R
       8      1     29      1      11 launch_qc_dashboard.default@21-31@schardExplorer/R/launch.R
       5      1     55      4       8 launch_qc_dashboard.character@31-38@schardExplorer/R/launch.R
       5      1     55      4       8 launch_qc_dashboard.list@38-45@schardExplorer/R/launch.R
       5      1     55      4       8 launch_qc_dashboard.SingleCellExperiment@45-52@schardExplorer/R/launch.R
       4      1     54      4       4 launch_qc_dashboard.Seurat@52-55@schardExplorer/R/launch.R
      28      6    333      3      39 NULLplot_umap_qc@13-51@schardExplorer/R/plotting.R
      22      4    314      2      32 plot_qc_distributions@51-82@schardExplorer/R/plotting.R
      21      1    271      3      22 plot_impact_barplot@82-103@schardExplorer/R/plotting.R
      23      9    339      2      28 compute_qc_metrics@10-37@schardExplorer/R/qc_engine.R
       2      2     29      2      15 (anonymous)@37-51@schardExplorer/R/qc_engine.R
      38     16    484      2      53 belseaapply_filter_thresholds@51-103@schardExplorer/R/qc_engine.R
      28      3    312      4      31 impact_by_replicate@103-133@schardExplorer/R/qc_engine.R
      33      7    356      2      38 detect_qc_cols@8-45@schardExplorer/R/detect_qc_cols.R
      23      7    167      1      28 normalize_to_list@8-35@schardExplorer/R/utils.R
       3      1     25      1       3 format_pct@35-37@schardExplorer/R/utils.R
9 file analyzed.
==============================================================
NLOC    Avg.NLOC  AvgCCN  Avg.token  function_cnt    file
--------------------------------------------------------------
    263      52.8     8.4      595.8         5     schardExplorer/R/shiny_server.R
    161      23.9     3.6      215.6         7     schardExplorer/R/report.R
     15       7.5     1.0       72.5         2     schardExplorer/R/shiny_app.R
     44      44.0     2.0      328.0         1     schardExplorer/R/shiny_ui.R
     26       5.2     1.0       44.0         6     schardExplorer/R/launch.R
     70      23.7     3.7      306.0         3     schardExplorer/R/plotting.R
     88      22.8     7.5      291.0         4     schardExplorer/R/qc_engine.R
     33      33.0     7.0      356.0         1     schardExplorer/R/detect_qc_cols.R
     25      13.0     4.0       96.0         2     schardExplorer/R/utils.R

===========================================================================================================
!!!! Warnings (cyclomatic_complexity > 15 or length > 1000 or nloc > 1000000 or parameter_count > 100) !!!!
================================================
  NLOC    CCN   token  PARAM  length  location  
------------------------------------------------
     238     36   2664      7     275 shiny_server@1-275@schardExplorer/R/shiny_server.R
      38     16    484      2      53 belseaapply_filter_thresholds@51-103@schardExplorer/R/qc_engine.R
==========================================================================================
Total nloc   Avg.NLOC  AvgCCN  Avg.token   Fun Cnt  Warning cnt   Fun Rt   nloc Rt
------------------------------------------------------------------------------------------
       725      23.9     4.3      253.4       31            2      0.06    0.37
