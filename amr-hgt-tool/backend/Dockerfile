# Use the exact RGI image already referenced by AMRN's local workflow.
FROM finlaymaguire/rgi:latest

USER root
WORKDIR /app

COPY requirements.txt /app/requirements.txt
RUN python -m pip install --no-cache-dir -r /app/requirements.txt

COPY . /app

# wrapper.py intentionally remains unchanged. It still calls `docker run`
# exactly as it does locally. Inside Render there is no nested Docker daemon,
# so this small compatibility shim maps that existing invocation back to the
# RGI executable already bundled in this same image.
RUN install -m 0755 /app/docker-rgi-shim /usr/local/bin/docker

ENV PYTHONUNBUFFERED=1

CMD ["sh", "-c", "gunicorn --bind 0.0.0.0:${PORT:-10000} --timeout 300 app:app"]
