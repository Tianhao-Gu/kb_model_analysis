FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

# RUN apt-get update

RUN pip install --upgrade pip \
    && python --version

RUN pip install scipy==1.4.1 \
    && pip install pandas==1.0.0 \
    && pip install matplotlib==3.1.2 \
    && pip install openpyxl==3.0.3 \
    && pip install xlrd==1.2.0 \
    && pip install scikit-learn==0.22.1 \
    && pip install mock==3.0.5

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
