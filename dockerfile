FROM rocker/shiny:4.0.5

RUN apt-get update
RUN apt-get install -y libcurl4-gnutls-dev \
   libssl-dev \
   libxml2-dev \
   nano


RUN R -q -e "install.packages('plotly')" && \
    R -q -e "install.packages('deSolve',repos='http://cran.rstudio.com/')" && \
    R -q -e "install.packages('tidyverse')" && \
    R -q -e "install.packages('shinythemes')" && \
    R -q -e "install.packages('EpiModel',repos='http://cran.rstudio.com/')" && \
    R -q -e "install.packages('survival',repos='http://cran.rstudio.com/')" && \
    R -q -e "install.packages('EasyABC',repos='http://cran.rstudio.com/')" && \
    R -q -e "install.packages('viridis',repos='http://cran.rstudio.com/')" 
    
ADD . /srv/shiny-server/shiny

# run app
#CMD ["/usr/bin/shiny-server"]