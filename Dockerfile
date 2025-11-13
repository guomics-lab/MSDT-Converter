FROM python:3.13.9

WORKDIR /app

COPY . /app/

RUN pip install -r /app/requirements.txt -i https://mirrors.aliyun.com/pypi/simple/

ENTRYPOINT ["python", "convert.py"]