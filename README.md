# covid19-severity
In this repository you can find code used to extract COVID-19 biomarkers and develop classification system 

code.R is the script that has been used to carry out the experiments

Data of the different GEO series can be found at -> https://drive.google.com/drive/folders/1FRNNoIVBP64MAd4cMIkEEwcOK-p8ypxh?usp=sharing

For each serie counts and severity information is available. Severity data for GSE156063  and  GSE152075 was obtained through their respective author.

# SHINI APP

app.R: code to contruct your local intuitive interface to carry out the experiments. To run the application download index.RData, www folder, app.R, labels.csv (severity labels of the samples) and matrix.csv (expression matrix of the three GEO series) from https://drive.google.com/drive/folders/1GGNyHkSocVvMstkHwgJ0BG_LX0MRAA33?usp=sharing and run it on R. (All files should be in the same folder)
