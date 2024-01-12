# Specify parent image. Please select a fixed tag here.
ARG BASE_IMAGE=jupyter/minimal-notebook:python-3.9.13
FROM $BASE_IMAGE

# Fix: https://github.com/hadolint/hadolint/wiki/DL4006
# Fix: https://github.com/koalaman/shellcheck/wiki/SC3014
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

USER root

RUN apt-get update --yes && \
    apt-get install --yes --no-install-recommends \
    build-essential texlive-latex-recommended texlive-latex-extra texlive-fonts-recommended dvipng cm-super&& \
    apt-get clean && rm -rf /var/lib/apt/lists/*


# Add user for mybinder
ARG NB_USER=jovyan
ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}

RUN chown -R ${NB_USER} ${HOME}
USER ${NB_USER}

# Install Python 3 packages
COPY requirements.txt ${HOME}
RUN pip install -r "$HOME/requirements.txt" && rm "$HOME/requirements.txt"

# Copy ipynb + images
COPY DualQuaternionQuinticBlends ${HOME}/DualQuaternionQuinticBlends
COPY interpolator_example.ipynb ${HOME}
# ----

RUN pip cache purge 
