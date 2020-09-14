FROM rocker/tidyverse:4.0.0
WORKDIR /tda
RUN apt-get update \
    && apt-get -y install libgsl0-dev \
    && apt -y install libxml2-dev
# renv and R packages
ENV RENV_VERSION 0.9.3
RUN echo "options(renv.consent = TRUE)" >> .Rprofile
COPY renv.lock .
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"
# restore from copied lockfile
RUN R -e "renv::restore(confirm = FALSE)"
