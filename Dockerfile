# syntax=docker/dockerfile:1.7
FROM mambaorg/micromamba:2.3.2

LABEL org.opencontainers.image.title="PhageMine" \
      org.opencontainers.image.source="https://github.com/emmannaemeka/PhageMine" \
      org.opencontainers.image.version="1.3.0.dev0"

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
# The optional CA secret is for managed environments with a trusted proxy.
# Ordinary Docker builds use the base image's public CA store.
RUN --mount=type=secret,id=proxy_ca \
    if [ -f /run/secrets/proxy_ca ]; then export SSL_CERT_FILE=/run/secrets/proxy_ca PIP_CERT=/run/secrets/proxy_ca CURL_CA_BUNDLE=/run/secrets/proxy_ca CONDA_SSL_VERIFY=/run/secrets/proxy_ca; fi; \
    micromamba install --yes --name base --file /tmp/environment.yml && \
    micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/phagemine
WORKDIR /opt/phagemine
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN --mount=type=secret,id=proxy_ca \
    if [ -f /run/secrets/proxy_ca ]; then export PIP_CERT=/run/secrets/proxy_ca; fi; \
    python -m pip install --no-cache-dir . && \
    python -m pip uninstall --yes textwrap3 backports.tempfile && python -m pip check

ENV PHAGEMINE_REGISTRY_PATH=/data/resources.json
WORKDIR /data
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "phagemine"]
CMD ["--help"]
