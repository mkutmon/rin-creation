# download relevant species file from http://mirtarbase.mbc.nctu.edu.tw/php/download.php
# create conversion.jar file - run ant in RINCreator  

java -cp RINCreator/dist/conversion.jar cytargetlinker.conversion.MirTarBase -i hsa_MTI.csv --bridgeDbFile Hs_Derby_20120602.bridge mirbase-v19.bridge -o mirtarbase-hsa-3.5.xgmml --organism "Homo sapiens" -v "3.5"

