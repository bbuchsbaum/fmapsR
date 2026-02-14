#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
with_python <- "--with-python" %in% args

required <- c(
  "Matrix",
  "RSpectra",
  "testthat",
  "bench",
  "devtools",
  "roxygen2",
  "pkgload",
  "knitr",
  "rmarkdown"
)

installed <- rownames(installed.packages())
missing <- setdiff(required, installed)

if (length(missing) == 0) {
  message("All R development packages are already installed.")
} else {
  message("Installing missing R packages: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
}

if (!isTRUE(with_python)) {
  message("Skipping Python parity environment setup (pass --with-python to enable).")
  quit(status = 0)
}

find_python <- function() {
  candidates <- c(
    Sys.getenv("PYFM_BENCH_PYTHON", unset = ""),
    Sys.which("python3"),
    Sys.which("python")
  )
  candidates <- unique(candidates[nzchar(candidates)])

  for (cand in candidates) {
    if (nzchar(cand) && (file.exists(cand) || nzchar(Sys.which(cand)))) {
      return(cand)
    }
  }
  ""
}

run_cmd <- function(cmd, args) {
  out <- tryCatch(
    system2(cmd, args, stdout = TRUE, stderr = TRUE),
    error = function(e) paste0("system2_failure: ", conditionMessage(e))
  )
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) {
    stop(sprintf("Command failed: %s %s\n%s", cmd, paste(args, collapse = " "), paste(out, collapse = "\n")), call. = FALSE)
  }
  invisible(out)
}

base_py <- find_python()
if (!nzchar(base_py)) {
  stop("No Python interpreter found for parity setup. Install python3 or set PYFM_BENCH_PYTHON.", call. = FALSE)
}

venv_dir <- ".venv-bench"
venv_py <- file.path(venv_dir, "bin", "python")

if (!file.exists(venv_py)) {
  message("Creating virtual environment: ", venv_dir)
  run_cmd(base_py, c("-m", "venv", venv_dir))
}

message("Installing parity dependencies into ", venv_dir)
run_cmd(venv_py, c("-m", "pip", "install", "--upgrade", "pip"))
run_cmd(venv_py, c("-m", "pip", "install", "numpy", "scipy", "scikit-learn", "tqdm"))

message("Python parity environment ready: ", normalizePath(venv_py, winslash = "/", mustWork = FALSE))
message("Tip: export PYFM_BENCH_PYTHON=", normalizePath(venv_py, winslash = "/", mustWork = FALSE))
