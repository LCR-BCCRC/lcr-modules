FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"
RUN mamba create --yes --name env \
        --channel conda-forge \
        --channel defaults \
        "python=3.12" \
        "numpy" \
        "requests" \
        "pyyaml" \
        "pip" \
    && mamba clean --all --yes \
    && /opt/conda/envs/env/bin/pip install \
        "vquest" \
        "lxml_html_clean"
ENV PATH=/opt/conda/envs/env/bin:$PATH
CMD ["/bin/bash"]
