FROM docker.1ms.run/conda/miniconda3:latest

WORKDIR /app

COPY . .

RUN conda create -n longbow python=3.7 -y && \
conda init bash && \
. ~/.bashrc && \
conda activate longbow && \
pip install --no-cache-dir . && \
longbow --version && \
echo "LongBow test success" && \
echo "conda activate longbow" >> ~/.bashrc && \
. ~/.bashrc && \
echo "Auto activate written to bashrc"

WORKDIR /root



