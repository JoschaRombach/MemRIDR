The app is still under development and not all functionalities currently work.
The data file Array_object.rds has not been added yet so app does not work as provided on github. 
If you want access email joscharom@outlook.com and I will share the file with you.

Frontpage:
![image](https://github.com/JoschaRombach/Proteome-wide_map_of_membrane_binding/assets/153042844/9088bc27-9fad-4d97-ae31-1dba8ea0d637)


Data page:
![image](https://github.com/JoschaRombach/Proteome-wide_map_of_membrane_binding/assets/153042844/f8dec209-dfd1-4dd8-94b7-59d360c226ef)

Instructions:
1. Download R (version > 4.2.1, https://www.r-project.org/) and Rstudio (https://posit.co/download/rstudio-desktop/)
2. Download all files as .zip, place it anywhere on your computer and unzip
3. Navigate to the file "app.R" and open it
4. Press the "Run App" button (green arrow) in the top right corner of the script panel
First time use only:
5. The app will prompt you to allow for download of packages from CRAN, Bioconductor and Github:

If you want to install packages seperately use "install.packages()" to download:
'tidyverse','berryFunctions','seqinr','stringr','data.table','viridis','bio3d',
'ggprism','ggExtra','gridExtra','zoo','RColorBrewer','shiny','DT','shiny.fluent',
'imola','cowplot','shiny.router','shiny.react','shinyWidgets','circlize'

Use "BiocManager::install()" to install: 'msa' and 'ggtree'

Use "devtools::install_github("swsoyee/r3dmol", upgrade = 'never')" to install 'r3dmol'
