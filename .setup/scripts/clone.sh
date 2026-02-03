#!/bin/bash

git clone git@github.com:jeksterslab/metaDyn.git
rm -rf "$PWD.git"
mv metaDyn/.git "$PWD"
rm -rf metaDyn
