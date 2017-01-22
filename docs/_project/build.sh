#!/bin/bash

if ! type sphinx-build >/dev/null 2>&1; then
    echo "sphinx is not installed"
    exit 1
fi

#sphinx-apidoc -f -o ./apis ../../fmflow
sphinx-build -a -d ../_doctree ./ ../
sphinx-build -a -d ../_doctree ./ ../
