# download relevant species file from http://mirecords.biolead.org/download.php
# download identifier mapping database for gene products from http://bridgedb.org/data/gene_database/
# download identifier mapping database for miRNAs from https://github.com/mkutmon/rin-creation/blob/master/mirna-mapping/mirbase-v19.bridge
# create conversion.jar file - run ant in RINCreator  

java -cp conversion.jar cytargetlinker.conversion.MiRecords -i miRecords_version4.csv --bridgeDbFile Hs_Derby_20120602.bridge mirbase-v19.bridge -o mirecords-hsa-2013-04-27.xgmml --organism "Homo sapiens" -v "2013-04-27"
