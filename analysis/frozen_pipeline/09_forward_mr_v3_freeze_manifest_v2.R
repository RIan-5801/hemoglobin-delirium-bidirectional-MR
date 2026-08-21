args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") stop("Usage: Rscript 09_forward_mr_v3_freeze_manifest_v2.R --project-root <path>", call. = FALSE)
root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
qc_dir <- file.path(root,"results","qc"); log_dir <- file.path(root,"results","logs")
manifest_path <- file.path(qc_dir,"vuckovic_hb_finngen_r13_forward_mr_v3_freeze_manifest_v2.csv")
freeze_path <- file.path(qc_dir,"vuckovic_hb_finngen_r13_forward_mr_v3_freeze_v2.json")
log_path <- file.path(log_dir,"vuckovic_hb_finngen_r13_forward_mr_v3_freeze_v2.log")
stop_if <- function(x,msg) if(isTRUE(x)) stop(msg,call.=FALSE)
hash_file <- function(p) digest::digest(file=p,algo="sha256",serialize=FALSE)
log_line <- function(...) cat(paste0(...),"\n",file=log_path,append=TRUE)
numeric_equal <- function(x,y) is.finite(x)&&is.finite(y)&&abs(as.numeric(x)-as.numeric(y))<=1e-12
required_ivw_fields <- c("analysis_set","analysis_role","method","nsnp","beta","se","pval","ci_lower_beta","ci_upper_beta","OR","OR_lci","OR_uci")

normalize_single_row_object <- function(x, object_name) {
  if (is.data.frame(x)) {
    stop_if(nrow(x)!=1L,paste0(object_name," data.frame has ",nrow(x)," rows.")); out <- x
  } else if (is.list(x) && length(x)==1L && is.list(x[[1L]]) && is.null(names(x))) {
    return(normalize_single_row_object(x[[1L]],object_name))
  } else if (is.list(x) && all(required_ivw_fields %in% names(x))) {
    lens <- vapply(x[required_ivw_fields],length,integer(1)); stop_if(any(lens!=1L),paste0(object_name," has non-singleton required field."))
    out <- as.data.frame(lapply(x[required_ivw_fields],function(z) z[[1L]]),stringsAsFactors=FALSE,check.names=FALSE)
  } else stop(paste0("Cannot uniquely normalize ",object_name,"; class=",paste(class(x),collapse="|")," names=",paste(names(x),collapse="|")),call.=FALSE)
  stop_if(!all(required_ivw_fields %in% names(out)),paste0(object_name," lacks required named fields."))
  out <- out[,required_ivw_fields,drop=FALSE]
  out$analysis_set <- as.character(out$analysis_set); out$analysis_role <- as.character(out$analysis_role); out$method <- as.character(out$method)
  for (nm in setdiff(required_ivw_fields,c("analysis_set","analysis_role","method"))) out[[nm]] <- as.numeric(out[[nm]])
  stop_if(anyNA(out[,required_ivw_fields]),paste0(object_name," contains NA after normalization.")); out
}

synthetic_normalizer_test <- function() {
  base <- list(analysis_set="set",analysis_role="role",method="method",nsnp=3L,beta=0.1,se=0.2,pval=0.3,ci_lower_beta=-0.2,ci_upper_beta=0.4,OR=1.1,OR_lci=0.8,OR_uci=1.5)
  a <- normalize_single_row_object(base,"named_scalar_list")
  b <- normalize_single_row_object(lapply(base,function(z) c(z)),"length_one_vectors")
  c <- normalize_single_row_object(unname(list(as.data.frame(base,stringsAsFactors=FALSE))),"single_row_data_frame_like")
  identical(names(a),required_ivw_fields) && identical(names(b),required_ivw_fields) && identical(names(c),required_ivw_fields) && nrow(a)==1L && nrow(b)==1L && nrow(c)==1L && identical(a,b) && identical(a,c)
}

main <- function() {
  output_targets <- c(manifest_path,freeze_path,log_path)
  stop_if(any(file.exists(c(output_targets,paste0(output_targets,".partial")))),"Freeze V2 target or partial exists; refusing to overwrite.")
  files <- data.frame(file_role=c("forward_v3_script","forward_v3_decision_35","frozen_mr_input_copy_included","frozen_mr_input_copy_excluded","forward_v3_mr_estimates","forward_v3_heterogeneity","forward_v3_egger_intercept","forward_v3_instrument_strength","forward_v3_instrument_strength_summary","forward_v3_leave_one_out","forward_v3_leave_one_out_full_ivw","forward_v3_mr_presso","forward_v3_qc","forward_v3_log","upstream_harmonisation_v4_qc","upstream_loo_label_audit_v4_qc","authoritative_harmonised_input_included_v4","authoritative_harmonised_input_excluded_v4"),relative_path=c("R/09_forward_mr_vuckovic_hb_finngen_r13_primary_v3.R","docs/decisions/35_vuckovic_hb_finngen_r13_forward_mr_v3_v1.1.md","data_derived/mr_inputs/vuckovic_hb_finngen_r13_forward_primary_apoe_included_v3.parquet","data_derived/mr_inputs/vuckovic_hb_finngen_r13_forward_primary_apoe_excluded_v3.parquet","results/tables/vuckovic_hb_finngen_r13_forward_mr_estimates_v3.csv","results/tables/vuckovic_hb_finngen_r13_forward_heterogeneity_v3.csv","results/tables/vuckovic_hb_finngen_r13_forward_egger_intercept_v3.csv","results/tables/vuckovic_hb_finngen_r13_forward_instrument_strength_v3.csv","results/tables/vuckovic_hb_finngen_r13_forward_instrument_strength_summary_v3.csv","results/tables/vuckovic_hb_finngen_r13_forward_leave_one_out_v3.csv","results/tables/vuckovic_hb_finngen_r13_forward_leave_one_out_full_ivw_v3.csv","results/tables/vuckovic_hb_finngen_r13_forward_mr_presso_v3.csv","results/qc/vuckovic_hb_finngen_r13_forward_mr_v3.json","results/logs/vuckovic_hb_finngen_r13_forward_mr_v3.log","results/qc/vuckovic_hb_finngen_r13_primary_harmonisation_v4.json","results/qc/vuckovic_hb_finngen_r13_forward_loo_label_audit_v4.json","data_derived/harmonised/vuckovic_hb_finngen_r13_primary_apoe_included_v4.parquet","data_derived/harmonised/vuckovic_hb_finngen_r13_primary_apoe_excluded_v4.parquet"),stringsAsFactors=FALSE)
  paths <- file.path(root,files$relative_path); all_present <- all(file.exists(paths)); stop_if(!all_present,"Authoritative Forward V3 file is absent.")
  unit <- synthetic_normalizer_test(); stop_if(!unit,"single-row JSON normalizer synthetic unit test failed.")
  qc_path <- file.path(root,"results/qc/vuckovic_hb_finngen_r13_forward_mr_v3.json"); qc <- jsonlite::fromJSON(qc_path,simplifyVector=FALSE)
  for(nm in c("primary_ivw","apoe_excluded_ivw")){x<-qc[[nm]];log_line(nm,"_class=",paste(class(x),collapse="|"));log_line(nm,"_names=",paste(names(x),collapse="|"));capture.output(str(x),file=log_path,append=TRUE)}
  primary_qc <- normalize_single_row_object(qc$primary_ivw,"primary_ivw"); apoe_qc <- normalize_single_row_object(qc$apoe_excluded_ivw,"apoe_excluded_ivw")
  est <- read.csv(file.path(root,"results/tables/vuckovic_hb_finngen_r13_forward_mr_estimates_v3.csv"),check.names=FALSE,stringsAsFactors=FALSE)
  primary_csv <- est[est$analysis_set=="APOE_included" & est$method_id=="mr_ivw",,drop=FALSE]; apoe_csv <- est[est$analysis_set=="APOE_excluded" & est$method_id=="mr_ivw",,drop=FALSE]
  compare_ivw <- function(q,csv){nums<-c("nsnp","beta","se","pval","ci_lower_beta","ci_upper_beta","OR","OR_lci","OR_uci"); chars<-c("analysis_set","analysis_role","method"); diffs<-stats::setNames(vapply(nums,function(nm)abs(q[[nm]][[1L]]-as.numeric(csv[[nm]][[1L]])),numeric(1)),nums); list(pass=identical(as.character(q$analysis_set[[1L]]),as.character(csv$analysis_set[[1L]]))&&identical(as.character(q$analysis_role[[1L]]),as.character(csv$analysis_role[[1L]]))&&identical(as.character(q$method[[1L]]),as.character(csv$method[[1L]]))&&all(vapply(nums,function(nm)numeric_equal(q[[nm]][[1L]],csv[[nm]][[1L]]),logical(1))),absolute_difference=as.list(diffs))}
  primary_cmp <- if(nrow(primary_csv)==1L) compare_ivw(primary_qc,primary_csv) else list(pass=FALSE,absolute_difference=list()); apoe_cmp <- if(nrow(apoe_csv)==1L) compare_ivw(apoe_qc,apoe_csv) else list(pass=FALSE,absolute_difference=list())
  source_roles <- c("harmonised_included","harmonised_excluded","mr_input_included","mr_input_excluded","estimates","qc","script","log","decision_35")
  source_paths <- file.path(root,c("data_derived/harmonised/vuckovic_hb_finngen_r13_primary_apoe_included_v4.parquet","data_derived/harmonised/vuckovic_hb_finngen_r13_primary_apoe_excluded_v4.parquet","data_derived/mr_inputs/vuckovic_hb_finngen_r13_forward_primary_apoe_included_v3.parquet","data_derived/mr_inputs/vuckovic_hb_finngen_r13_forward_primary_apoe_excluded_v3.parquet","results/tables/vuckovic_hb_finngen_r13_forward_mr_estimates_v3.csv","results/qc/vuckovic_hb_finngen_r13_forward_mr_v3.json","R/09_forward_mr_vuckovic_hb_finngen_r13_primary_v3.R","results/logs/vuckovic_hb_finngen_r13_forward_mr_v3.log","docs/decisions/35_vuckovic_hb_finngen_r13_forward_mr_v3_v1.1.md")); names(source_paths)<-source_roles; source_before<-as.list(vapply(source_paths,hash_file,character(1)))
  hinc<-source_before$harmonised_included; hexc<-source_before$harmonised_excluded; qinc<-qc$harmonised_input_sha256$included_after; qexc<-qc$harmonised_input_sha256$excluded_after
  files$file_size_bytes<-as.numeric(file.info(paths)$size); files$sha256<-vapply(paths,hash_file,character(1))
  write.csv(files,paste0(manifest_path,".partial"),row.names=FALSE,na="")
  mr <- read.csv(paste0(manifest_path,".partial"),stringsAsFactors=FALSE,check.names=FALSE); manifest_complete<-nrow(mr)==18L && identical(sort(mr$relative_path),sort(files$relative_path)) && all(mr$file_size_bytes>0) && all(grepl("^[0-9a-f]{64}$",mr$sha256))
  manifest_sha<-hash_file(paste0(manifest_path,".partial"))
  source_after<-as.list(vapply(source_paths,hash_file,character(1))); no_source_mutation<-identical(source_before,source_after)
  checks<-list(all_authoritative_files_present=all_present,forward_v3_qc_passed=identical(qc$mr_status,"passed"),forward_v3_approved=isTRUE(qc$approved_for_forward_results_interpretation),no_forward_v3_hard_check_failures=length(qc$hard_check_failures)==0L,single_row_json_normalizer_unit_test=unit,primary_ivw_qc_normalized=nrow(primary_qc)==1L,apoe_excluded_ivw_qc_normalized=nrow(apoe_qc)==1L,primary_ivw_csv_unique=nrow(primary_csv)==1L,apoe_excluded_ivw_csv_unique=nrow(apoe_csv)==1L,qc_csv_primary_ivw_consistent=isTRUE(primary_cmp$pass),qc_csv_apoe_excluded_ivw_consistent=isTRUE(apoe_cmp$pass),included_harmonised_input_sha_matches_v3_qc=identical(tolower(hinc),tolower(qinc)),excluded_harmonised_input_sha_matches_v3_qc=identical(tolower(hexc),tolower(qexc)),mr_input_copy_included_present=file.exists(source_paths[["mr_input_included"]]),mr_input_copy_excluded_present=file.exists(source_paths[["mr_input_excluded"]]),mr_input_copy_included_sha_recorded=nzchar(source_before$mr_input_included),mr_input_copy_excluded_sha_recorded=nzchar(source_before$mr_input_excluded),manifest_complete=manifest_complete,manifest_paths_unique=!anyDuplicated(mr$relative_path),manifest_sha_complete=grepl("^[0-9a-f]{64}$",manifest_sha),no_new_partial_after_success=TRUE,no_source_mutation=no_source_mutation)
  failures<-names(checks)[!vapply(checks,isTRUE,logical(1))]; status<-if(length(failures)==0L)"passed" else "failed"
  freeze<-list(freeze_version="v2",authoritative_forward_mr_version="v3",analysis_direction="Hb_to_delirium",primary_analysis_set="APOE_included",sensitivity_analysis_set="APOE_excluded",forward_mr_status=qc$mr_status,approved_for_forward_results_interpretation=qc$approved_for_forward_results_interpretation,primary_ivw=primary_qc,apoe_excluded_ivw=apoe_qc,instrument_strength_summary=qc$instrument_strength_summary,heterogeneity_summary=qc$heterogeneity_summary,egger_intercept_summary=qc$egger_intercept_summary,mr_presso_summary=qc$mr_presso_summary,effect_interpretation="per genetically predicted 1-unit increase in standardized inverse-normal-transformed haemoglobin",steiger_status="deferred_to_directionality_sensitivity_stage",authoritative_harmonised_inputs=list(included_path="data_derived/harmonised/vuckovic_hb_finngen_r13_primary_apoe_included_v4.parquet",included_sha256=hinc,excluded_path="data_derived/harmonised/vuckovic_hb_finngen_r13_primary_apoe_excluded_v4.parquet",excluded_sha256=hexc),frozen_mr_input_copies=list(included_path="data_derived/mr_inputs/vuckovic_hb_finngen_r13_forward_primary_apoe_included_v3.parquet",included_sha256=source_before$mr_input_included,excluded_path="data_derived/mr_inputs/vuckovic_hb_finngen_r13_forward_primary_apoe_excluded_v3.parquet",excluded_sha256=source_before$mr_input_excluded),authoritative_files=files,manifest_sha256=manifest_sha,source_sha_before=source_before,source_sha_after=source_after,ivw_qc_csv_absolute_differences=list(primary=primary_cmp$absolute_difference,apoe_excluded=apoe_cmp$absolute_difference),freeze_status=status,hard_checks=checks,hard_check_failures=failures,informational_findings=list(freeze_v1="failed due to file-role SHA comparison and IVW JSON-shape handling errors; retained as historical evidence"))
  writeLines(jsonlite::toJSON(freeze,auto_unbox=TRUE,pretty=TRUE,digits=NA),paste0(freeze_path,".partial"))
  stop_if(status!="passed",paste0("Freeze V2 hard checks failed: ",paste(failures,collapse="; ")))
  stop_if(!file.rename(paste0(manifest_path,".partial"),manifest_path),"Manifest atomic rename failed."); stop_if(!identical(hash_file(manifest_path),manifest_sha),"Published manifest SHA differs from partial manifest SHA."); stop_if(!file.rename(paste0(freeze_path,".partial"),freeze_path),"Freeze JSON atomic rename failed.")
  stop_if(file.exists(paste0(manifest_path,".partial"))||file.exists(paste0(freeze_path,".partial")),"New V2 partial remains after publish.")
  log_line("freeze_status=passed");log_line("manifest_sha256=",manifest_sha);log_line("hard_check_failures=")
}
status<-tryCatch({main();0L},error=function(e){log_line("freeze_status=failed");log_line("error=",conditionMessage(e));1L});quit(status=status)
