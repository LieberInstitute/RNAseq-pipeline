#!/bin/bash

echo "**** Cleaning job starts ****"
date

rm */.*.sh */.*.R */.paired_end
rm -fr */Coverage */logs */FastQC

echo "**** Cleaning job ends ****"
date
