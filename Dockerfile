FROM continuumio/miniconda3:23.10.0-1

WORKDIR /app

RUN conda install -y python=3.8
RUN conda install -y gdal=3.6.0
#RUN conda install -y pandas=2.0.1
RUN conda install -y -c conda-forge julia=1.7 awscli
RUN python -m pip install awscli; pip install Pillow; pip install pystac==1.8.4; pip install spectral==0.23.1

# Get SpectralUnmixing repo
RUN git clone https://github.com/EnSpec/SpectralUnmixing.git -b v0.2.1-sister

# Get Unity-Py 
RUN pip install git+https://github.com/unity-sds/unity-py.git

# Get hytools
WORKDIR /app
RUN git clone https://github.com/EnSpec/hytools.git
WORKDIR /app/hytools
RUN pip install -e .

# Activate JULIA
WORKDIR /app/SpectralUnmixing
RUN export JULIA_SSL_CA_ROOTS_PATH=""
RUN julia -e 'using Pkg; Pkg.activate("."); Pkg.add(path="https://github.com/kmsquire/ArgParse2.jl"); Pkg.instantiate()'
RUN julia -e 'using Pkg; Pkg.activate("."); Pkg.precompile()'
ENV JULIA_PROJECT="/app/SpectralUnmixing"


COPY data data/.

COPY process.py .

#CMD ["python", "process.py", "catalog.json", "/app/SBG-unity-frac-cover/outputs/SBG-L2-FRAC_COVER/", "1", "1.0", "none", "001", "True", "SBG-L2-Fractional-Cover"]
