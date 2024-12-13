
This image contains:

 - R 4.0.3
 - Rstudio server (installation requires the userconf.sh file)


# ######################
     COMPILE THE IMAGE
# ######################

docker build -t seurat_lung_cancer /mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/03_WORKSPACE/rebuffet/Lung_Basel/01_ALL_SAMPLES/02_Container/seurat_lung_cancer

# ######################
     RUN THE IMAGE
# ######################

docker run -d --name seurat_lung_cancer -p 8881:8787 -v /mnt:/mnt -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g) -e PASSWORD=coucou seurat_lung_cancer
