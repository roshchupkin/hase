FROM python:3.13-slim@sha256:70729b46c69b4f1e97c4822c1af3df53a1476cf5ddc6c087c0c10bc3a5678c2f
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
COPY requirements.lock.txt /opt/pilot/requirements.lock.txt
RUN python -m pip install --disable-pip-version-check --no-cache-dir --only-binary=:all: -r /opt/pilot/requirements.lock.txt && python -m pip freeze > /opt/pilot/freeze.txt
WORKDIR /tmp
