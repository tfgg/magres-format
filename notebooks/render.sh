echo "---
permalink: magres-python/tutorials/structural/
layout: page_magres
title:  \"MagresPython: Structural properties\"
date:   2014-10-03
category: magres-python
categoryname: \"MagresPython\"
navorder: 1
---" > ~/Dropbox/personal_website/jekyll/_pages/magres-python/tutorials/magres-tutorial-structural.markdown

ipython nbconvert --template=backticks.tpl  --to markdown Structural\ properties.ipynb --stdout >> ~/Dropbox/personal_website/jekyll/_pages/magres-python/tutorials/magres-tutorial-structural.markdown
