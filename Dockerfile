FROM jozefhajnala/devr

# Clone the repo into RStudio's home for fast development
RUN git clone \
      https://github.com/SPan-18/spStack-dev \
      /home/rstudio/spStack-dev \
 && chown rstudio:rstudio /home/rstudio/spStack-dev -R \
 && R CMD INSTALL /home/rstudio/spStack-dev
