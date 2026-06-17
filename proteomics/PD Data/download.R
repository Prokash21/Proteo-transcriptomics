## ----DownloadLabelfreeData----------------------------------------------------
url2 <- "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD019103/proteinGroups.txt"
options(timeout = 600)
download.file(url2, destfile = "./proteinGroups.txt", method = "auto")
#unzip("proteinGroups.txt")