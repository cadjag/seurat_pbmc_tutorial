FROM rocker/shiny:4.2.1
ENV RENV_CONFIG_REPOS_OVERRIDE https://packagemanager.rstudio.com/cran/latest

RUN apt-get update -qq && \ 
  apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libglpk-dev \
    libgmp3-dev \
    libicu-dev \
    libpng-dev \
    libssl-dev \
    libxml2-dev \
    make \
    pandoc \
    python3 \
    zlib1g-dev && \
  apt-get clean && \ 
  rm -rf /var/lib/apt/lists/*
COPY shiny_renv.lock renv.lock
RUN Rscript -e "install.packages('renv')"
RUN Rscript -e "renv::restore()"
COPY ./app/* /srv/shiny-server/
USER shiny
EXPOSE 3838
CMD ["/usr/bin/shiny-server"]
