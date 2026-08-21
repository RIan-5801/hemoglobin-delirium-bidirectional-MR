#!/usr/bin/env Rscript

# Validates two independent V2 PLINK clumping results; does not invoke PLINK.
suppressPackageStartupMessages({ library(DBI); library(duckdb); library(jsonlite); library(digest) })

argv <- commandArgs(trailingOnly=TRUE)
parse_cli <- function(x) {
  ans <- list(); i <- 1L
  while (i <= length(x)) {
    if (!startsWith(x[[i]],"-") || i == length(x)) stop("Invalid command-line arguments",call.=FALSE)
    ans[[sub("^-+","",x[[i]])]] <- x[[i+1L]]
    i <- i + 2L
  }
  ans
}
opt <- parse_cli(argv)
for (k in c("Root","PlinkPath","PlinkVersion")) if (is.null(opt[[k]])) stop("Missing argument: ",k,call.=FALSE)
root <- normalizePath(opt$Root,winslash="/",mustWork=TRUE)
plink_path <- normalizePath(opt$PlinkPath,winslash="/",mustWork=TRUE)
plink_version <- opt$PlinkVersion
included_input <- file.path(root,"data_derived/clumping_inputs/vuckovic_hb_clump_input_apoe_included_v1.tsv")
excluded_input <- file.path(root,"data_derived/clumping_inputs/vuckovic_hb_clump_input_apoe_excluded_v1.tsv")
included_prefix <- file.path(root,"data_derived/clumped/vuckovic_hb_apoe_included_v2")
excluded_prefix <- file.path(root,"data_derived/clumped/vuckovic_hb_apoe_excluded_v2")
pipeline_log <- file.path(root,"results/logs/vuckovic_hb_clumping_pipeline_v2.log")
outputs <- list(
  included_parquet=file.path(root,"data_derived/instruments/vuckovic_hb_instruments_apoe_included_v2.parquet"),
  included_tsv=file.path(root,"data_derived/instruments/vuckovic_hb_instruments_apoe_included_v2.tsv"),
  excluded_parquet=file.path(root,"data_derived/instruments/vuckovic_hb_instruments_apoe_excluded_v2.parquet"),
  excluded_tsv=file.path(root,"data_derived/instruments/vuckovic_hb_instruments_apoe_excluded_v2.tsv"),
  qc_json=file.path(root,"results/qc/vuckovic_hb_clumping_v2.json"),
  comparison_csv=file.path(root,"results/qc/vuckovic_hb_clumping_comparison_v2.csv"),
  by_chr_csv=file.path(root,"results/qc/vuckovic_hb_clumped_by_chr_v2.csv")
)
abort_if <- function(x,msg) if(isTRUE(x)) stop(msg,call.=FALSE)
for (p in c(included_input,excluded_input,paste0(included_prefix,".clumps"),paste0(excluded_prefix,".clumps"),paste0(included_prefix,".log"),paste0(excluded_prefix,".log"),pipeline_log)) abort_if(!file.exists(p),paste("Missing required file:",p))
for (p in outputs) { abort_if(file.exists(p),paste("Refusing to overwrite:",p)); abort_if(file.exists(paste0(p,".partial")),paste("Residual partial:",paste0(p,".partial"))); dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE) }
sha256 <- function(p) if(file.exists(p)) digest(p,algo="sha256",file=TRUE) else NA_character_
atomic_write <- function(p,fn) { q<-paste0(p,".partial"); fn(q); abort_if(!file.exists(q),paste("Writer failed:",q)); abort_if(!file.rename(q,p),paste("Promotion failed:",p)) }
allele_set <- function(a,b) mapply(function(x,y) paste(sort(c(toupper(x),toupper(y))),collapse=":"),a,b,USE.NAMES=FALSE)
read_input <- function(p) {
  x <- read.delim(p,sep="\t",check.names=FALSE,stringsAsFactors=FALSE)
  required <- c("rsid","chr","pos","effect_allele","other_allele","beta","se","eaf","lp_raw","pval","F_stat","apoe_region","source_row_number","variant_id_source")
  abort_if(!all(required %in% names(x)),paste("Input missing required columns:",p))
  x$apoe_region <- tolower(as.character(x$apoe_region))=="true"
  x
}
read_clumps <- function(p) {
  x <- read.table(p,header=TRUE,comment.char="",check.names=FALSE,stringsAsFactors=FALSE,fill=TRUE)
  abort_if(nrow(x)==0,paste("Zero index SNP rows:",p))
  normalized <- toupper(gsub("[^A-Za-z0-9]","",names(x)))
  hit <- match(c("ID","SNP","INDEXSNP"),normalized,nomatch=0L); hit<-hit[hit>0L]
  abort_if(length(hit)==0,paste("Cannot identify index-SNP column:",p))
  x$index_snp <- as.character(x[[hit[[1L]]]])
  abort_if(any(is.na(x$index_snp)|x$index_snp==""),paste("Missing index SNP:",p))
  x
}
count_missing <- function(p) if(!file.exists(p)) 0L else sum(nzchar(trimws(readLines(p,warn=FALSE))))
assert_log_clean <- function(p) abort_if(any(grepl("\\b(ERROR|FATAL)\\b",readLines(p,warn=FALSE),ignore.case=TRUE,perl=TRUE)),paste("PLINK log contains ERROR/FATAL:",p))
validate_result <- function(name,input_path,prefix,exclude_apoe) {
  input <- read_input(input_path); clumps <- read_clumps(paste0(prefix,".clumps")); assert_log_clean(paste0(prefix,".log"))
  abort_if(anyDuplicated(clumps$index_snp),paste("Duplicate index SNP:",name))
  ix <- match(clumps$index_snp,input$rsid); abort_if(any(is.na(ix)),paste("Index SNP not in original input:",name))
  inst <- input[ix,,drop=FALSE]
  abort_if(anyDuplicated(inst$rsid),paste("Duplicate joined instrument:",name))
  abort_if(any(is.na(inst$pval)|inst$pval>=5e-8),paste("Index SNP outside P threshold:",name))
  abort_if(any(is.na(inst$F_stat)|inst$F_stat<10),paste("Weak index SNP:",name))
  abort_if(any(!grepl("^(?:[1-9]|1[0-9]|2[0-2])$",as.character(inst$chr))),paste("Non-autosomal index SNP:",name))
  abort_if(any(!grepl("^[ACGT]$",inst$effect_allele)|!grepl("^[ACGT]$",inst$other_allele)),paste("Non-ACGT index SNP:",name))
  if(exclude_apoe) abort_if(any(inst$apoe_region),paste("APOE index SNP in excluded result:",name))
  inst$clumping_analysis <- name
  list(input=input,clumps=clumps,instruments=inst,missing_id=count_missing(paste0(prefix,".clumps.missing_id")),missing_allele=count_missing(paste0(prefix,".clumps.missing_allele")))
}
write_parquet <- function(x,p,con) atomic_write(p,function(q) { dbWriteTable(con,"instrument_write",x,overwrite=TRUE); dbExecute(con,sprintf("COPY instrument_write TO '%s' (FORMAT PARQUET)",gsub("'","''",gsub("\\\\","/",q)))); dbExecute(con,"DROP TABLE instrument_write") })
write_tsv <- function(x,p) atomic_write(p,function(q) write.table(x,q,sep="\t",quote=FALSE,row.names=FALSE,na=""))
write_csv <- function(x,p) atomic_write(p,function(q) write.csv(x,q,row.names=FALSE,na=""))
write_json <- function(x,p) atomic_write(p,function(q) jsonlite::write_json(x,q,pretty=TRUE,auto_unbox=TRUE,null="null"))

included <- validate_result("apoe_included",included_input,included_prefix,FALSE)
excluded <- validate_result("apoe_excluded",excluded_input,excluded_prefix,TRUE)
shared <- intersect(included$instruments$rsid,excluded$instruments$rsid)
included_only <- setdiff(included$instruments$rsid,excluded$instruments$rsid)
excluded_only <- setdiff(excluded$instruments$rsid,included$instruments$rsid)
instrument_fields <- c("rsid","chr","pos","effect_allele","other_allele","beta","se","eaf","lp_raw","pval","F_stat","apoe_region","source_row_number","variant_id_source","clumping_analysis")
con <- dbConnect(duckdb::duckdb(),dbdir=":memory:",config=list(shared_home=FALSE)); on.exit(dbDisconnect(con,shutdown=TRUE),add=TRUE)
write_parquet(included$instruments[,instrument_fields],outputs$included_parquet,con)
write_tsv(included$instruments[,instrument_fields],outputs$included_tsv)
write_parquet(excluded$instruments[,instrument_fields],outputs$excluded_parquet,con)
write_tsv(excluded$instruments[,instrument_fields],outputs$excluded_tsv)
by_chr <- rbind(data.frame(analysis="apoe_included",chr=names(table(included$instruments$chr)),index_snp_count=as.integer(table(included$instruments$chr))),data.frame(analysis="apoe_excluded",chr=names(table(excluded$instruments$chr)),index_snp_count=as.integer(table(excluded$instruments$chr))))
comparison <- data.frame(shared_index_snps=length(shared),included_only_index_snps=length(included_only),excluded_only_index_snps=length(excluded_only))
write_csv(by_chr,outputs$by_chr_csv); write_csv(comparison,outputs$comparison_csv)

roundtrip_ids <- function(parquet_path,tsv_path) {
  p <- dbGetQuery(con,sprintf("SELECT rsid FROM read_parquet('%s')",gsub("'","''",gsub("\\\\","/",parquet_path))))$rsid
  t <- read.delim(tsv_path,sep="\t",check.names=FALSE,stringsAsFactors=FALSE)$rsid
  identical(sort(p),sort(t))
}
abort_if(!roundtrip_ids(outputs$included_parquet,outputs$included_tsv),"Included Parquet/TSV rsID sets differ")
abort_if(!roundtrip_ids(outputs$excluded_parquet,outputs$excluded_tsv),"Excluded Parquet/TSV rsID sets differ")
stats <- function(x) list(index_snps=nrow(x),f_min=min(x$F_stat),f_median=median(x$F_stat),f_mean=mean(x$F_stat),f_max=max(x$F_stat),apoe_index_snps=sum(x$apoe_region))
qc <- list(
  protocol="docs/protocol/analysis_plan_v1.1.md",
  execution_model="two independent PLINK V2 runs; excluded results are not derived by deleting APOE from included results",
  plink_version=plink_version,
  parameters=list(r2=0.001,window_kb=10000,log10p1=7.301029995663981,log10p2=0,unphased=TRUE,threads=8,memory_mb=8000,clump_allow_overlap=FALSE),
  included=c(stats(included$instruments),list(missing_id=included$missing_id,missing_allele=included$missing_allele,index_snps=as.list(included$instruments$rsid))),
  excluded=c(stats(excluded$instruments),list(missing_id=excluded$missing_id,missing_allele=excluded$missing_allele,index_snps=as.list(excluded$instruments$rsid))),
  comparison=list(shared=length(shared),included_only=length(included_only),excluded_only=length(excluded_only),interpretation="Differences may extend beyond APOE because each input underwent independent greedy clumping."),
  parquet_tsv_sets_identical=list(included=TRUE,excluded=TRUE),
  hashes=list(plink_executable=sha256(plink_path),included_input=sha256(included_input),excluded_input=sha256(excluded_input),included_clumps=sha256(paste0(included_prefix,".clumps")),excluded_clumps=sha256(paste0(excluded_prefix,".clumps")),included_log=sha256(paste0(included_prefix,".log")),excluded_log=sha256(paste0(excluded_prefix,".log")),included_missing_id=sha256(paste0(included_prefix,".clumps.missing_id")),excluded_missing_id=sha256(paste0(excluded_prefix,".clumps.missing_id")),included_missing_allele=sha256(paste0(included_prefix,".clumps.missing_allele")),excluded_missing_allele=sha256(paste0(excluded_prefix,".clumps.missing_allele")),included_instruments_parquet=sha256(outputs$included_parquet),included_instruments_tsv=sha256(outputs$included_tsv),excluded_instruments_parquet=sha256(outputs$excluded_parquet),excluded_instruments_tsv=sha256(outputs$excluded_tsv))
)
write_json(qc,outputs$qc_json)
cat(paste0(format(Sys.time(),tz="Asia/Shanghai",usetz=TRUE)," V2 validation completed successfully.\n"),file=pipeline_log,append=TRUE)
message("V2 validation completed successfully.")

