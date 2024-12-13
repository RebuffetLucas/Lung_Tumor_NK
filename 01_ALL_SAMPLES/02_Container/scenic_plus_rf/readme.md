
# 2024_06_27_lucasrebuffet_jupyter_conda_scenicplus

Based on jupyter/minimal-notebook
Add:
 - mamba environment
 - scenicplus



## Build

docker build -t scenic_plus_rf .



## Save

docker save scenic_plus_rf | gzip > scenic_plus_rf.tar.gz



## Run Jupyter web server

```
docker run -u $(id -u ${USER}):$(id -g ${USER}) \
           -p 9999:8888 \
           -e TOKEN=myPass \
           -v /mnt:/mnt \
           -v ~/:/host_home \
           scenic_plus_rf
```

Then browse to the machine running docker on mapped port (9999).
For local hosting: http://127.0.0.1:9999

Use value specified in 'TOKEN' as password


