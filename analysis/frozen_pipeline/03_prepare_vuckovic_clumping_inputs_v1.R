#!/usr/bin/env Rscript

# Creates audited Vuckovic Hb clumping inputs. This script never invokes PLINK.
suppressPackageStartupMessages({ library(DBI); library(duckdb); library(jsonlite); library(digest) })

root <- "E:/Research/hb_delirium_bidir_mr"
apoe_chr <- "19"; apoe_start <- 44000000L; apoe_end <- 46000000L
expected <- list(eligible=87991L, apoe=164L, conflicts=13L, coordinate=10L, allele=3L,
                 duplicate_records=78L, duplicate_groups=39L, resolved=4L, unresolved_groups=35L,
                 included=87982L, excluded=87818L)
inputs <- list(
  candidates=file.path(root,"data_derived/instruments/vuckovic_hb_p5e8_candidates_v3.parquet"),
  eligible=file.path(root,"data_derived/instruments/vuckovic_hb_clump_eligible_v1.parquet"),
  ineligible=file.path(root,"data_derived/instruments/vuckovic_hb_clump_ineligible_v1.parquet"),
  bim=file.path(root,"resources/ld/1kg_v3/EUR.bim"),
  bed=file.path(root,"resources/ld/1kg_v3/EUR.bed"),
  fam=file.path(root,"resources/ld/1kg_v3/EUR.fam"),
  protocol=file.path(root,"docs/protocol/analysis_plan_v1.1.md"),
  decision=file.path(root,"docs/decisions/05_ld_variant_identity_rule_v1.1.md")
)
outputs <- list(
  conflicts=file.path(root,"data_derived/instruments/reference_conflicts/vuckovic_hb_ld_identity_conflicts_v1.parquet"),
  audit=file.path(root,"data_derived/instruments/duplicates/vuckovic_hb_duplicate_rsid_audit_v1.parquet"),
  resolved=file.path(root,"data_derived/instruments/duplicates/vuckovic_hb_duplicate_rsid_resolved_v1.parquet"),
  unresolved=file.path(root,"data_derived/instruments/duplicates/vuckovic_hb_duplicate_rsid_unresolved_v1.parquet"),
  included=file.path(root,"data_derived/clumping_inputs/vuckovic_hb_clump_input_apoe_included_v1.tsv"),
  excluded=file.path(root,"data_derived/clumping_inputs/vuckovic_hb_clump_input_apoe_excluded_v1.tsv"),
  duplicate_json=file.path(root,"results/qc/vuckovic_hb_duplicate_resolution_v1.json"),
  conflict_json=file.path(root,"results/qc/vuckovic_hb_ld_identity_conflicts_v1.json"),
  counts_csv=file.path(root,"results/qc/vuckovic_hb_clumping_input_counts_v1.csv"),
  inputs_json=file.path(root,"results/qc/vuckovic_hb_clumping_inputs_v1.json"),
  log=file.path(root,"results/logs/vuckovic_hb_prepare_clumping_inputs_v1.log")
)

abort_if <- function(x, msg) if (isTRUE(x)) stop(msg, call.=FALSE)
for (p in inputs) abort_if(!file.exists(p), paste("Missing input:", p))
for (p in outputs) {
  abort_if(file.exists(p), paste("Refusing to overwrite output:", p))
  abort_if(file.exists(paste0(p,".partial")), paste("Refusing residual partial output:", paste0(p,".partial")))
}
for (p in outputs) dir.create(dirname(p), recursive=TRUE, showWarnings=FALSE)

qpath <- function(x) gsub("'", "''", gsub("\\\\", "/", x))
allele_set <- function(a, b) mapply(function(x,y) paste(sort(c(toupper(x),toupper(y))),collapse=":"), a,b,USE.NAMES=FALSE)
is_apoe <- function(x) as.character(x$chr)==apoe_chr & x$pos>=apoe_start & x$pos<=apoe_end
is_acgt <- function(x) grepl("^[ACGT]$",x$effect_allele) & grepl("^[ACGT]$",x$other_allele)
sha256 <- function(x) digest(x, algo="sha256", file=TRUE)
atomic_write <- function(path, fun) {
  partial <- paste0(path,".partial")
  abort_if(file.exists(path)||file.exists(partial), paste("Output occupied during write:",path))
  fun(partial)
  abort_if(!file.exists(partial), paste("Writer did not create:",partial))
  abort_if(!file.rename(partial,path), paste("Atomic promotion failed:",path))
}

con <- dbConnect(duckdb::duckdb(), dbdir=":memory:", config=list(shared_home=FALSE))
on.exit(dbDisconnect(con, shutdown=TRUE), add=TRUE)
dbExecute(con,sprintf("CREATE TEMP VIEW candidates AS SELECT * FROM read_parquet('%s')",qpath(inputs$candidates)))
dbExecute(con,sprintf("CREATE TEMP VIEW eligible AS SELECT * FROM read_parquet('%s')",qpath(inputs$eligible)))
dbExecute(con,sprintf("CREATE TEMP VIEW ineligible AS SELECT * FROM read_parquet('%s')",qpath(inputs$ineligible)))
dbExecute(con,sprintf("CREATE TEMP VIEW bim AS SELECT * FROM read_csv('%s',delim='\\t',header=false,columns={'chr':'VARCHAR','rsid':'VARCHAR','cm':'DOUBLE','pos':'BIGINT','a1':'VARCHAR','a2':'VARCHAR'})",qpath(inputs$bim)))

candidates <- dbGetQuery(con,"SELECT * FROM candidates")
eligible <- dbGetQuery(con,"SELECT * FROM eligible")
ineligible <- dbGetQuery(con,"SELECT * FROM ineligible")
dbWriteTable(con,"candidate_ids",data.frame(rsid=unique(candidates$rsid),stringsAsFactors=FALSE),overwrite=TRUE)
refs <- dbGetQuery(con,"SELECT b.chr,b.rsid,b.pos,b.a1,b.a2 FROM bim b INNER JOIN candidate_ids i USING(rsid)")
bim_duplicate_groups <- dbGetQuery(con,"SELECT count(*) AS n FROM (SELECT rsid FROM bim WHERE rsid <> '.' GROUP BY rsid HAVING count(*)>1)")$n[[1]]
abort_if(bim_duplicate_groups!=0,"EUR.bim contains duplicate rsIDs")

refs$reference_allele_set <- allele_set(refs$a1,refs$a2)
add_reference <- function(x) {
  x$apoe_region <- is_apoe(x)
  x$gwas_allele_set <- allele_set(x$effect_allele,x$other_allele)
  i <- match(x$rsid,refs$rsid)
  x$reference_chr <- refs$chr[i]; x$reference_pos <- refs$pos[i]
  x$reference_a1 <- refs$a1[i]; x$reference_a2 <- refs$a2[i]
  x$reference_allele_set <- refs$reference_allele_set[i]
  x$reference_identity_match <- !is.na(x$reference_chr) & x$chr==x$reference_chr & x$pos==x$reference_pos & x$gwas_allele_set==x$reference_allele_set
  x
}
candidates <- add_reference(candidates)
eligible <- add_reference(eligible)

abort_if(nrow(eligible)!=expected$eligible,"Eligible count differs from expectation")
abort_if(sum(eligible$apoe_region)!=expected$apoe,"Eligible APOE count differs from expectation")
abort_if(anyDuplicated(eligible$rsid),"Eligible input contains duplicate rsIDs")

conflicts <- eligible[!eligible$reference_identity_match,,drop=FALSE]
conflicts$identity_status <- "reference_identity_conflict"
conflicts$conflict_reason <- ifelse(conflicts$chr!=conflicts$reference_chr | conflicts$pos!=conflicts$reference_pos,"reference_coordinate_conflict","reference_allele_conflict")
conflicts$review_status <- "unresolved_excluded_from_current_clumping"
conflicts_out <- conflicts[,c("rsid","source_row_number","chr","pos","effect_allele","other_allele","gwas_allele_set","reference_chr","reference_pos","reference_a1","reference_a2","reference_allele_set","identity_status","conflict_reason","apoe_region","review_status")]
names(conflicts_out)[match(c("chr","pos","effect_allele","other_allele"),names(conflicts_out))] <- c("gwas_chr","gwas_pos","gwas_effect_allele","gwas_other_allele")
abort_if(nrow(conflicts_out)!=expected$conflicts,"Reference identity-conflict count differs from expectation")
abort_if(sum(conflicts_out$conflict_reason=="reference_coordinate_conflict")!=expected$coordinate,"Reference coordinate-conflict count differs from expectation")
abort_if(sum(conflicts_out$conflict_reason=="reference_allele_conflict")!=expected$allele,"Reference allele-conflict count differs from expectation")

dup_counts <- table(candidates$rsid)
duplicate_rsids <- names(dup_counts[dup_counts>1])
duplicates <- candidates[candidates$rsid %in% duplicate_rsids,,drop=FALSE]
abort_if(nrow(duplicates)!=expected$duplicate_records,"Duplicate-record count differs from expectation")
abort_if(length(duplicate_rsids)!=expected$duplicate_groups,"Duplicate-group count differs from expectation")
abort_if(sum(ineligible$rsid %in% duplicate_rsids)!=expected$duplicate_records,"Duplicate records are not all retained in ineligible")
duplicates$duplicate_group_size <- as.integer(dup_counts[duplicates$rsid])
duplicates$duplicate_group_class <- NA_character_
duplicates$duplicate_resolution_status <- NA_character_
duplicates$duplicate_resolution_reason <- NA_character_
duplicates$reference_identity_status <- ifelse(duplicates$reference_identity_match,"exact_reference_identity_match",ifelse(is.na(duplicates$reference_chr),"reference_rsid_absent","reference_identity_mismatch"))
selected <- logical(nrow(duplicates))
for (id in duplicate_rsids) {
  ix <- which(duplicates$rsid==id)
  matched <- ix[duplicates$reference_identity_match[ix]]
  if (length(matched)==1L) {
    duplicates$duplicate_group_class[ix] <- "unique_reference_match"
    duplicates$duplicate_resolution_status[ix] <- "unresolved_not_selected"
    duplicates$duplicate_resolution_reason[ix] <- "same_rsid_member_not_unique_reference_match"
    duplicates$duplicate_resolution_status[matched] <- "resolved_unique_reference_match"
    duplicates$duplicate_resolution_reason[matched] <- "only_group_member_matching_rsid_chr_pos_unordered_alleles"
    selected[matched] <- TRUE
  } else {
    duplicates$duplicate_group_class[ix] <- "coordinate_or_allele_conflict"
    duplicates$duplicate_resolution_status[ix] <- "unresolved_coordinate_or_allele_conflict"
    duplicates$duplicate_resolution_reason[ix] <- "no_unique_group_member_matching_rsid_chr_pos_unordered_alleles"
  }
}
abort_if(sum(selected)!=expected$resolved,"Resolved duplicate count differs from expectation")
abort_if(length(unique(duplicates$rsid[duplicates$duplicate_group_class=="coordinate_or_allele_conflict"]))!=expected$unresolved_groups,"Unresolved duplicate-group count differs from expectation")
audit_cols <- c("rsid","source_row_number","chr","pos","effect_allele","other_allele","allele_set_key","gwas_allele_set","beta","se","eaf","log10p","lp_raw","pval","F_stat","apoe_region","variant_id_source","reference_chr","reference_pos","reference_a1","reference_a2","reference_allele_set","reference_identity_status","duplicate_group_size","duplicate_group_class","duplicate_resolution_status","duplicate_resolution_reason")
duplicate_audit <- duplicates[,audit_cols]
duplicate_resolved <- duplicate_audit[selected,,drop=FALSE]
duplicate_unresolved <- duplicate_audit[!selected,,drop=FALSE]

base <- eligible[eligible$reference_identity_match,,drop=FALSE]
recovered <- candidates[candidates$source_row_number %in% duplicates$source_row_number[selected],,drop=FALSE]
recovered$eligibility_status <- "eligible_recovered_duplicate"; recovered$exclusion_reason <- NA_character_
abort_if(!all(recovered$reference_identity_match),"Recovered duplicate failed identity rule")
base$duplicate_resolution_status <- "not_duplicate"; base$reference_identity_status <- "exact_reference_identity_match"
recovered$duplicate_resolution_status <- "resolved_unique_reference_match"; recovered$reference_identity_status <- "exact_reference_identity_match"
final_source_columns <- c(names(eligible), "duplicate_resolution_status", "reference_identity_status")
final_source <- rbind(base[,final_source_columns],recovered[,final_source_columns])
final_source$SNP <- final_source$rsid; final_source$CHR <- final_source$chr; final_source$BP <- final_source$pos; final_source$LOG10P <- final_source$lp_raw
output_cols <- c("rsid","chr","pos","effect_allele","other_allele","beta","se","eaf","log10p","lp_raw","pval","F_stat","apoe_region","duplicate_resolution_status","reference_identity_status","source_row_number","variant_id_source","SNP","CHR","BP","LOG10P")
included <- final_source[,output_cols]
included <- included[order(-included$LOG10P,included$SNP),,drop=FALSE]
excluded <- included[!included$apoe_region,,drop=FALSE]

identity_matches <- function(x) {
  i <- match(x$rsid,refs$rsid)
  !is.na(refs$rsid[i]) & x$chr==refs$chr[i] & x$pos==refs$pos[i] & allele_set(x$effect_allele,x$other_allele)==refs$reference_allele_set[i]
}
validate <- function(x,n,apoe_n) list(
  records=nrow(x), unique_rsid=length(unique(x$rsid)), duplicate_rsid_records=sum(duplicated(x$rsid)|duplicated(x$rsid,fromLast=TRUE)),
  rsid_missing=sum(is.na(x$rsid)|x$rsid==""), eur_rsid_missing=sum(is.na(match(x$rsid,refs$rsid))), reference_identity_mismatch=sum(!identity_matches(x)),
  non_autosomal=sum(!grepl("^(?:[1-9]|1[0-9]|2[0-2])$",x$chr)), non_acgt_snp=sum(!is_acgt(x)), f_lt_10=sum(x$F_stat<10), p_not_lt_5e8=sum(is.na(x$pval)|x$pval>=5e-8),
  apoe_records=sum(x$apoe_region), expected_records=n, expected_apoe=apoe_n)
included_check <- validate(included,expected$included,expected$apoe)
excluded_check <- validate(excluded,expected$excluded,0L)
for (check in list(included_check,excluded_check)) {
  abort_if(check$records!=check$expected_records||check$unique_rsid!=check$records||check$duplicate_rsid_records!=0,"Final input count or rsID-uniqueness validation failed")
  abort_if(any(unlist(check[c("rsid_missing","eur_rsid_missing","reference_identity_mismatch","non_autosomal","non_acgt_snp","f_lt_10","p_not_lt_5e8")])!=0),"Final input eligibility validation failed")
  abort_if(check$apoe_records!=check$expected_apoe,"Final input APOE count differs from expectation")
}
non_apoe_included <- sort(included$rsid[!included$apoe_region])
abort_if(!identical(non_apoe_included,sort(excluded$rsid)),"Non-APOE rsID sets differ between inputs")
abort_if(!identical(sort(setdiff(included$rsid,excluded$rsid)),sort(included$rsid[included$apoe_region])),"Input difference is not exactly the APOE rsID set")
abort_if(any(conflicts_out$rsid %in% included$rsid)|any(conflicts_out$rsid %in% excluded$rsid),"Reference-conflict rsID entered a clumping input")
unresolved_rsids <- unique(duplicate_unresolved$rsid[duplicate_unresolved$duplicate_group_class=="coordinate_or_allele_conflict"])
abort_if(any(unresolved_rsids %in% included$rsid)|any(unresolved_rsids %in% excluded$rsid),"Unresolved duplicate rsID entered a clumping input")
abort_if(!all(table(included$rsid[included$rsid %in% recovered$rsid])==1L)|!all(table(excluded$rsid[excluded$rsid %in% recovered$rsid])==1L),"Recovered rsID occurrence validation failed")

write_parquet <- function(x,path) atomic_write(path,function(partial) { dbWriteTable(con,"write_table",x,overwrite=TRUE); dbExecute(con,sprintf("COPY write_table TO '%s' (FORMAT PARQUET)",qpath(partial))); dbExecute(con,"DROP TABLE write_table") })
write_tsv <- function(x,path) atomic_write(path,function(partial) write.table(x,partial,sep="\t",quote=FALSE,row.names=FALSE,na=""))
write_json_atomic <- function(x,path) atomic_write(path,function(partial) write_json(x,partial,pretty=TRUE,auto_unbox=TRUE,null="null"))
write_csv_atomic <- function(x,path) atomic_write(path,function(partial) write.csv(x,partial,row.names=FALSE,na=""))

write_parquet(conflicts_out,outputs$conflicts)
write_parquet(duplicate_audit,outputs$audit)
write_parquet(duplicate_resolved,outputs$resolved)
write_parquet(duplicate_unresolved,outputs$unresolved)
write_tsv(included,outputs$included)
write_tsv(excluded,outputs$excluded)
included_sha <- sha256(outputs$included); excluded_sha <- sha256(outputs$excluded)

duplicate_summary <- list(input_duplicate_records=nrow(duplicate_audit),input_duplicate_groups=length(duplicate_rsids),class_counts=as.list(table(duplicate_audit$duplicate_group_class)),resolved_records=nrow(duplicate_resolved),unresolved_records=nrow(duplicate_unresolved),resolved_rsids=as.list(duplicate_resolved$rsid),rule="unique rsID+GRCh37 chr+pos+unordered allele-set reference match only")
conflict_summary <- list(records=nrow(conflicts_out),coordinate_conflicts=sum(conflicts_out$conflict_reason=="reference_coordinate_conflict"),allele_conflicts=sum(conflicts_out$conflict_reason=="reference_allele_conflict"),review_status="unresolved_excluded_from_current_clumping",rsids=as.list(conflicts_out$rsid))
input_summary <- list(protocol="docs/protocol/analysis_plan_v1.1.md",identity_rule="rsID, GRCh37 chromosome, GRCh37 position, and unordered allele set all match EUR.bim",apoe_policy=list(primary="included",sensitivity="excluded before independent clumping",interval="GRCh37 chr19:44000000-46000000"),included=c(included_check,list(file=outputs$included,sha256=included_sha)),excluded=c(excluded_check,list(file=outputs$excluded,sha256=excluded_sha)),non_apoe_sets_identical=identical(non_apoe_included,sort(excluded$rsid)),included_minus_excluded=length(setdiff(included$rsid,excluded$rsid)),plink_executed=FALSE,ld_clumping_executed=FALSE)
counts <- data.frame(metric=c("eligible_input","identity_conflicts_excluded","duplicate_records_audited","duplicate_records_resolved","duplicate_records_unresolved","apoe_included_input","apoe_excluded_input","included_apoe_records","excluded_apoe_records"),value=c(nrow(eligible),nrow(conflicts_out),nrow(duplicate_audit),nrow(duplicate_resolved),nrow(duplicate_unresolved),nrow(included),nrow(excluded),sum(included$apoe_region),sum(excluded$apoe_region)))
write_json_atomic(duplicate_summary,outputs$duplicate_json)
write_json_atomic(conflict_summary,outputs$conflict_json)
write_csv_atomic(counts,outputs$counts_csv)
write_json_atomic(input_summary,outputs$inputs_json)
log_lines <- c(paste("Started:",format(Sys.time(),tz="Asia/Shanghai",usetz=TRUE)),"Protocol: docs/protocol/analysis_plan_v1.1.md",paste("Eligible input records:",nrow(eligible)),paste("Identity conflicts excluded:",nrow(conflicts_out),"(coordinate",sum(conflicts_out$conflict_reason=="reference_coordinate_conflict"),"; allele",sum(conflicts_out$conflict_reason=="reference_allele_conflict"),")"),paste("Duplicate groups/records:",length(duplicate_rsids),"/",nrow(duplicate_audit)),paste("Resolved duplicate records:",nrow(duplicate_resolved)),paste("APOE-included input:",nrow(included),"records; APOE",sum(included$apoe_region),"; SHA-256",included_sha),paste("APOE-excluded input:",nrow(excluded),"records; APOE",sum(excluded$apoe_region),"; SHA-256",excluded_sha),"All pre-clumping validations passed.","PLINK was not run.","LD clumping was not run.",paste("Completed:",format(Sys.time(),tz="Asia/Shanghai",usetz=TRUE)))
atomic_write(outputs$log,function(partial) writeLines(log_lines,partial,useBytes=TRUE))
message("Completed successfully. No PLINK or LD clumping was run.")
