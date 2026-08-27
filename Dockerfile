FROM guomics2017/msdt-converter:v1.3

WORKDIR /app

COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r /app/requirements.txt
COPY . /app/

ENTRYPOINT ["python", "convert.py"]
