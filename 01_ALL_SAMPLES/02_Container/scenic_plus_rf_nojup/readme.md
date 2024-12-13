
# 2024_06_27_lucasrebuffet_jupyter_conda_scenicplus

Based on jupyter/minimal-notebook
Add:
 - mamba environment
 - scenicplus



## Build

docker build -t scenic_plus_rf_nojup .



## Save

docker save scenic_plus_rf_nojup | gzip > scenic_plus_rf_nojup.tar.gz



## Run Jupyter web server

```
docker run -u $(id -u ${USER}):$(id -g ${USER}) \
           -p 9996:8886 \
           -e TOKEN=myPass \
           -v /mnt:/mnt \
           -v ~/:/host_home \
           scenic_plus_rf_nojup
```

Then browse to the machine running docker on mapped port (9999).
For local hosting: http://127.0.0.1:9996

Use value specified in 'TOKEN' as password


