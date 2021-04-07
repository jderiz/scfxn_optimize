#!/bin/bash

# Merge sqlite3 databases into a single database.
# This can be used, for example with the separate_db_per_mpi_process flag to ReportToDB
#
# To merge the database parts 'features_<db_id>.db3_*' into the single database 'feature_<db_id>.db3'
#
#      bash merge.sh features_<db_id>.db3 feature_<db_id>.db3_*
#
#

### HAD TO SWITCH to /usr/bin/sqlite3 so the proper shared libraries were used libreadline.so.5
### HAD TO ADD LD_PRELOAD=/lib64/libsqlite3.so.0 to the from of each call to sqlite3 so header and source version mismatch don't occur

cache_size=4000000

target=$1
shift
inputs=$*
dumps=""
for x in $inputs
do
    echo "dumping $x ..."
    LD_PRELOAD=/lib64/libsqlite3.so.0 /usr/bin/sqlite3 $x '.dump' | sed -e 's/BEGIN TRANSACTION;//g' -e 's/COMMIT;//g' -e 's/CREATE TABLE/CREATE TABLE IF NOT EXISTS/g' -e 's/INSERT/INSERT OR IGNORE/g' > ${x}.dump
    rm $x
    dumps="$dumps ${x}.dump"
done

echo "merging dumps ..."

rm -rf $target
cat <(echo -e "PRAGMA foreign_keys=OFF;\nPRAGMA cache_size = ${cache_size};\nPRAGMA synchronous = OFF;\nPRAGMA journal_mode = OFF;\nPRAGMA locking_mode = EXCLUSIVE;\nPRAGMA count_changes = OFF;\nPRAGMA temp_store = MEMORY;\nPRAGMA auto_vacuum = NONE;\nBEGIN TRANSACTION;") $dumps <(echo "COMMIT;") | LD_PRELOAD=/lib64/libsqlite3.so.0 /usr/bin/sqlite3 $target

#cleanup dumps
rm -rf ${dumps}

echo "analyzing database ..."
LD_PRELOAD=/lib64/libsqlite3.so.0 /usr/bin/sqlite3 $target "analyze;" 
