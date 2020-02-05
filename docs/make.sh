DOCROOT=`pwd`

# user's manual
cd ${DOCROOT}/sphinx/ja
make html
cd ${DOCROOT}/sphinx/ja
make html

# api
cd $DOCROOT
sphinx-apidoc -f -e -o ./api ../abics
cd api
make html

