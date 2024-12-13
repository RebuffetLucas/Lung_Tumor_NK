
This image contains:

 - R 4.0.3
 - Rstudio server (installation requires the userconf.sh file)


# ######################
     COMPILE THE IMAGE
# ######################

docker build -t seuratv5_lung_cancer /mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/03_WORKSPACE/rebuffet/Lung_Basel/01_ALL_SAMPLES/02_Container/seuratv5_lung_cancer

# ######################
     RUN THE IMAGE
# ######################

docker run -d --name seuratv5_lung_cancer -p 8879:8787 -v /mnt:/mnt -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g) -e PASSWORD=coucou seuratv5_lung_cancer
