shiny_write_docker = function(
    path = ".", appdir = "app", lockfile = "shiny_renv.lock",
    port = 3838, expose = TRUE, rspm = TRUE
) {
  rspm_env = ifelse(
    rspm,
    "ENV RENV_CONFIG_REPOS_OVERRIDE https://packagemanager.rstudio.com/cran/latest\n",
    ""
  )
  from_shiny_version = glue::glue("FROM rocker/shiny:{getRversion()}")
  renv::snapshot(
    project = path,
    lockfile = lockfile,
    prompt = FALSE,
    force = TRUE
  )
  pkgs = renv::dependencies(appdir)$Package
  sys_reqs = glue_sys_reqs(pkgs)
  copy_renv = glue::glue("COPY {lockfile} renv.lock")
  renv_install = 'RUN Rscript -e "install.packages(\'renv\')"'
  renv_restore  = 'RUN Rscript -e "renv::restore()"'
  
  copy_app = glue::glue("COPY {appdir} /srv/shiny-server/")
  expose = ifelse(expose, glue::glue("EXPOSE {port}"), "")
  cmd = 'CMD ["/usr/bin/shiny-server"]'
  
  ret = purrr::compact(list(
    from_shiny_version,
    rspm_env,
    sys_reqs,
    copy_renv,
    renv_install,
    renv_restore,
    copy_app,
    expose,
    cmd
  ))
  readr::write_lines(ret, file = file.path(path, "Dockerfile"))
}