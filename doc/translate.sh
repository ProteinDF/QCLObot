#!/bin/bash

make gettext
cp -R _build/locale/*.pot local/pot
cp -R locale/pot/* locale/ja/LC_MESSAGES/
sphinx-intl update -p _build/local -l ja

