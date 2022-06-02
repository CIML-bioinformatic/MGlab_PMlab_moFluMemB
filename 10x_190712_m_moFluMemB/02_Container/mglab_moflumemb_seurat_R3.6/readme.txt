This image contains:

 - R 3.5.3
 - Rstudio server (installation requires the userconf.sh file)
 - Packages for the Seurat analysis
  

# ######################
     COMPILE THE IMAGE
# ######################

docker build -t mglab_moflumemb_seurat_r3_6 /mnt/NAS7/Workspace/spinellil/ciml-bip/project/MGlab/moFluMemB/10x_190712_m_moFluMemB/02_Container/mglab_moflumemb_seurat_R3.6

# ######################
     RUN THE IMAGE
# ######################

docker run --name mglab_moflumemb_seurat_r3_6 -d -p 8787:8787 -v /mnt:/mnt -e PASSWORD= -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g)  mglab_moflumemb_seurat_r3_6
 
# ######################
     CONNECT TO RSTUDIO
# ######################
 
 In an Internet browser, type as url : http://127.0.0.1:8787 and use the login/password: <user>/rstudio
 
# ######################
	 NOTES
# ######################
 
 - To use knitr PDF compilation instead of Sweave, you have to go into Rstudio menu Tools->Global Options->Sweave->Weave Rnw files with.. and select "knitr".
 
