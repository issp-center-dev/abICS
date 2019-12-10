#!/bin/sh

rm -rf build dist
python3 ./setup.py bdist_wheel
wheel=$(ls -1 dist/*.whl)

echo ""
cat << EOF
To install abics, type
pip install --user $wheel

To update abics, type
pip install --user --no-deps --force-reinstall $wheel
EOF

