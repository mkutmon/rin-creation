# download relevant species file from http://mirtarbase.mbc.nctu.edu.tw/php/download.php
# download identifier mapping database for gene products from http://bridgedb.org/data/gene_database/
# download identifier mapping database for miRNAs from https://github.com/mkutmon/rin-creation/blob/master/mirna-mapping/mirbase-v19.bridge
# create conversion.jar file - run ant in RINCreator  

java -cp RINCreator/dist/conversion.jar cytargetlinker.conversion.MirTarBase -i hsa_MTI.csv --bridgeDbFile Hs_Derby_20120602.bridge mirbase-v19.bridge -o mirtarbase-hsa-3.5.xgmml --organism "Homo sapiens" -v "3.5"

