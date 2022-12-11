DOCROOT=`pwd`

# user's manual
cd ${DOCROOT}/sphinx/ja
make html
# make latexpdf
cd ${DOCROOT}/sphinx/en
make html
# make latexpdf

# api
cd $DOCROOT
sphinx-apidoc -f -e -o ./api ../abics
cd api
make html

