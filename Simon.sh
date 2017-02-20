TRIBE_DIR="/scratch/cluster/monthly/ecabello/Simon/TRIBE/CODE/"
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A


CREATE USER 'username@'localhost' IDENTIFIED BY '';
GRANT ALL PRIVILEGES ON * . * TO 'username'@'localhost';
FLUSH PRIVILEGES;
#create mysql database
CREATE DATABASE mm9seq;
