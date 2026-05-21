library(stringr)
library(ape)
library(paleotree)
library(nloptr)
library(Matrix)
library(future.apply)
library(corpcor)
library(progressr)

suppressMessages(library(optparse))

source('/home/ec2-user/scripts/SCOUT//SCOUT_em_setup.R') 
source('/home/ec2-user/scripts/SCOUT/V2/SCOUT_EM_V4_TMP_RUNNER.R') ## SWITCHING TO TMP WHILE TROUBLESHOOTING 
source('/home/ec2-user/scripts/SCOUT/V2/SCOUT_EM_V4_utilities.R')

options(progressr.enable = TRUE)
handlers("progress")
#####################################################

option_list <- list(
  make_option(c("-s", "--samplesheet"),
              type = "character",
              help = "Model parameters."),

  make_option(c("-o", "--outdir"),
              type = "character",
              help = "path to output directory."),

  make_option(c("-l", "--logfile"),
              type = "character", default = NULL, 
              help = "path to logfile."),

  make_option(c("-c", "--cores"),
              type = "integer",
              help = "How many cores? Recommended max is 96.",
              default = 12),

  make_option(c("-r", "--regimes"),
              type = "character",
              default=NULL,
              help = "Comma-separated list of colnames with regime annotations. Defaults to c('BM1', 'OU1', 'OUM')."),

  make_option(c("-b", "--blacklist"),
              type = "character",
              default=NULL,
              help = "Comma-separated list of blacklist items to ignore when parsing metadata for gene columns.")
)

parser <- OptionParser(option_list = option_list)
options <- parse_args(parser)
#####################################################

cores <- options$cores
logfile <- options$logfile
ss <- read.csv(options$samplesheet,row.names=1)
wrkdir <- options$outdir
regimes <- if (is.null(options$regimes )) { c('BM1', 'OU1', 'OUM') } else {strsplit(options$regimes, ",")[[1]]}
blacklist <- if (is.null(options$blacklist )) { NULL } else {strsplit(options$blacklist, ",")[[1]]}
#####################################################

log_message(sprintf('Started logging @ %s\n', logfile, verbose = TRUE))
log_message(sprintf('Running sample sheet --> %s\n', options$samplesheet), logfile, verbose=TRUE)
log_message(sprintf('Data will be saved to --> %s\n', options$outdir), logfile, verbose=TRUE)
log_message('=========================================================\n', logfile, verbose=TRUE)
for (i in 1:nrow(ss)){

    counts <- read.csv(ss[i, 'data'],  row.names=1)
    tree <- ape::read.tree(ss[i, 'tree'])
    norm <- ss[i, 'normalize']
    scale <- ss[i, 'scale']

    idata <- formatSCOUT(tree_path = tree, 
                         metadata_path = counts, 
                         anc_infer = 'ape', 
                         outpath = wrkdir, 
                         regimes = regimes,
                         normalize = norm, 
                         scale = scale, 
                         smoothing = NULL, 
                         blacklist = blacklist, 
                         logfile=logfile)

    full.res <- runSCOUT.tester(idata, lambda1=ss[i, 'lambda1'], lambda2=ss[i, 'lambda2'], ss[i, 'fixed.root'], 
        M_only=ss[i, 'skipE'], runid=ss[i, 'testid'], cores = cores, logfile = logfile )

    final_list <- list(run_token = ss[i, 'testid'], test_settings = as.vector(ss[i,]), results = full.res )

    filename <- sprintf('%s/%s_%s.rds', wrkdir, code, format(Sys.Date(), "%Y%m%d"))
    saveRDS(final_list, filename)

    log_message(sprintf('Done with %s || %s / %s \n', ss[i, 'testid'], as.character(i), as.character(nrow(ss))), )
    log_message('=========================================================\n', logfile, verbose=TRUE)

}
